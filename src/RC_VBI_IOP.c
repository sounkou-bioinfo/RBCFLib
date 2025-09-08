#include <R.h>
#include <Rinternals.h>
#include <R_ext/Print.h>
#include <R_ext/Error.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "htslib/hts.h"
#include "htslib/bgzf.h"
#include "htslib/vcf.h"
#include "htslib/hfile.h"
#include "vbi_index_capi.h"
#include "cgranges.h"
// Forward declaration for finalizer (after SEXP is defined)
SEXP RC_cgranges_destroy(SEXP cr_ptr);
// Forward declarations from vbi_index.c and vbi_index_capi.h
int do_index(const char *infile, const char *outfile, int n_threads);
SEXP RC_VBI_index_memory_usage(SEXP extPtr);
// Minimal cr_chrom implementation for cgranges_t (do not modify cgranges library files)
// Reconstruct contig index for interval i by searching ctg ranges
static inline const char *cr_chrom(const cgranges_t *cr, int64_t i) {
    if (!cr || i < 0 || i >= cr->n_r) return NULL;
    for (int32_t j = 0; j < cr->n_ctg; ++j) {
        int64_t off = cr->ctg[j].off;
        int64_t n = cr->ctg[j].n;
        if (i >= off && i < off + n) {
            return cr->ctg[j].name;
        }
    }
    return NULL;
}

// Helper: check file existence
static int file_exists(const char *path) {
    FILE *f = fopen(path, "rb");
    if (f) { fclose(f); return 1; }
    return 0;
}

//' Create a VBI index for a VCF/BCF file
//' @param vcf_path Path to VCF/BCF file
//' @param vbi_path Path to output VBI index file
//' @param threads Number of threads
//' @return Path to the created VBI index file
SEXP RC_VBI_index(SEXP vcf_path, SEXP vbi_path, SEXP threads) {
    const char *vcf = CHAR(STRING_ELT(vcf_path, 0));
    const char *vbi = CHAR(STRING_ELT(vbi_path, 0));
    int n_threads = asInteger(threads);
    int ret = do_index(vcf, vbi, n_threads);
    if (ret != 0) {
        Rf_error("Failed to create VBI index for %s", vcf);
    }
    if (!file_exists(vbi)) {
        Rf_error("VBI index file not found after creation: %s", vbi);
    }
    Rprintf("VBI index created: %s\n", vbi);
    return vbi_path;
}

// Helper: download remote file to local path (R's download.file)
// we will get back to this
// avoiding the R C API is advised 
static int download_file(const char *url, const char *dest) {
    SEXP call = PROTECT(lang4(install("download.file"), mkString(url), mkString(dest), mkString("auto")));
    SEXP res = PROTECT(R_tryEval(call, R_GlobalEnv, NULL));
    int status = asInteger(res);
    UNPROTECT(2);
    return status == 0;
}

//' Query VCF/BCF by region using VBI index (returns only records, no header)
// the threads are decompression threads only
SEXP RC_VBI_query_range(SEXP vcf_path, SEXP idx_ptr, SEXP region, SEXP threads) {
    const char *vcf = CHAR(STRING_ELT(vcf_path, 0));
    vbi_index_t *idx = (vbi_index_t*) R_ExternalPtrAddr(idx_ptr);
    if (!idx) Rf_error("[VBI] Index pointer is NULL");
    const char *reg = CHAR(STRING_ELT(region, 0));
    int n_threads = asInteger(threads); (void)n_threads; // reserved, not yet used

    int nfound = 0;
    int *indices = vbi_index_query_region(idx, reg, &nfound);
    if (!indices || nfound == 0) {
        // return 0-row data.frame
        SEXP df = PROTECT(allocVector(VECSXP, 0));
        SEXP names = PROTECT(allocVector(STRSXP, 0));
        SEXP rn = PROTECT(allocVector(INTSXP, 0));
        setAttrib(df, R_NamesSymbol, names);
        setAttrib(df, R_RowNamesSymbol, rn);
        setAttrib(df, R_ClassSymbol, mkString("data.frame"));
        UNPROTECT(3);
        if (indices) free(indices);
        return df;
    }

    htsFile *fp = hts_open(vcf, "r");
    if (!fp) {
        free(indices);
        Rf_error("Failed to open VCF/BCF: %s", vcf);
    }
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (!hdr) {
        hts_close(fp);
        free(indices);
        Rf_error("Failed to read VCF/BCF header: %s", vcf);
    }

    // Detect presence of CSQ / ANN INFO annotations (VEP/SnpEff)
    bcf_hrec_t *csq_hrec = bcf_hdr_get_hrec(hdr, BCF_HL_INFO, "ID", "CSQ", NULL);
    bcf_hrec_t *ann_hrec = bcf_hdr_get_hrec(hdr, BCF_HL_INFO, "ID", "ANN", NULL);

    // Pre-parse CSQ/ANN subfield names if available
    char **csq_fields = NULL; int n_csq_fields = 0;
    char **ann_fields = NULL; int n_ann_fields = 0;

    // helper lambda-like static functions not available; implement inline parsing
    if (csq_hrec) {
        const char *desc = NULL;
        for (int i=0;i<csq_hrec->nkeys;i++) if (strcmp(csq_hrec->keys[i],"Description")==0) { desc = csq_hrec->vals[i]; break; }
        if (desc) {
            const char *fmt = strstr(desc, "Format:");
            if (fmt) {
                fmt += 7; // skip 'Format:'
                while (*fmt==' ' || *fmt=='\t') fmt++;
                // copy until end
                char *tmp = strdup(fmt);
                // trim trailing quote / parenthesis if any
                size_t L = strlen(tmp);
                while (L>0 && (tmp[L-1]=='"' || tmp[L-1]==')' || tmp[L-1]==' ')) { tmp[--L]='\0'; }
                // split by '|'
                // first count
                for (size_t i=0;i<L;i++) if (tmp[i]=='|') n_csq_fields++;
                n_csq_fields++; // items = pipes + 1
                csq_fields = (char**)calloc(n_csq_fields, sizeof(char*));
                int fld=0; char *p=tmp; char *start=p;
                for (; *p; ++p) {
                    if (*p=='|') { *p='\0'; csq_fields[fld++] = strdup(start); start = p+1; }
                }
                if (fld < n_csq_fields) csq_fields[fld++] = strdup(start);
                free(tmp);
            }
        }
    }
    if (ann_hrec) {
        const char *desc = NULL;
        for (int i=0;i<ann_hrec->nkeys;i++) if (strcmp(ann_hrec->keys[i],"Description")==0) { desc = ann_hrec->vals[i]; break; }
        if (desc) {
            const char *fmt = strstr(desc, "'|"); // heuristics to find ANN format after first "'"
            const char *fmt2 = strstr(desc, "Format:");
            const char *use = fmt2?fmt2:fmt; // try Format: first
            if (use) {
                if (fmt2) { use += 7; while (*use==' '||*use=='\t') use++; }
                char *tmp = strdup(use);
                size_t L = strlen(tmp);
                while (L>0 && (tmp[L-1]=='"' || tmp[L-1]==')' || tmp[L-1]==' ')) { tmp[--L]='\0'; }
                for (size_t i=0;i<L;i++) if (tmp[i]=='|') n_ann_fields++;
                n_ann_fields++;
                ann_fields = (char**)calloc(n_ann_fields, sizeof(char*));
                int fld=0; char *p=tmp; char *start=p;
                for (; *p; ++p) {
                    if (*p=='|') { *p='\0'; ann_fields[fld++] = strdup(start); start = p+1; }
                }
                if (fld < n_ann_fields) ann_fields[fld++] = strdup(start);
                free(tmp);
            }
        }
    }

    // Allocate core columns
    SEXP chrom_col = PROTECT(allocVector(STRSXP, nfound));
    SEXP pos_col   = PROTECT(allocVector(INTSXP, nfound));
    SEXP id_col    = PROTECT(allocVector(STRSXP, nfound));
    SEXP ref_col   = PROTECT(allocVector(STRSXP, nfound));
    SEXP alt_col   = PROTECT(allocVector(STRSXP, nfound));
    SEXP qual_col  = PROTECT(allocVector(REALSXP, nfound));
    SEXP filter_col= PROTECT(allocVector(STRSXP, nfound));
    SEXP nallele_col=PROTECT(allocVector(INTSXP, nfound));
    SEXP index_col = PROTECT(allocVector(INTSXP, nfound));
    SEXP csq_col   = csq_hrec ? PROTECT(allocVector(VECSXP, nfound)) : R_NilValue; int prot_csx = csq_hrec?1:0;
    SEXP ann_col   = ann_hrec ? PROTECT(allocVector(VECSXP, nfound)) : R_NilValue; int prot_ann = ann_hrec?1:0;

    int nprotect = 9 + prot_csx + prot_ann;

    bcf1_t *rec = bcf_init();
    kstring_t kalt = {0,0,0};
    kstring_t kflt = {0,0,0};

    for (int i = 0; i < nfound; ++i) {
        int idx_var = indices[i];
        int seek_ok = 0;
        if (fp->format.compression == bgzf) {
            BGZF *bg = (BGZF *)fp->fp.bgzf;
            seek_ok = (bgzf_seek(bg, idx->offsets[idx_var], SEEK_SET) == 0);
        } else {
            hFILE *hf = (hFILE *)fp->fp.hfile;
            seek_ok = (hseek(hf, (off_t)idx->offsets[idx_var], SEEK_SET) == 0);
        }
        if (!seek_ok || bcf_read(fp, hdr, rec) < 0) {
            SET_STRING_ELT(chrom_col, i, NA_STRING);
            INTEGER(pos_col)[i] = NA_INTEGER;
            SET_STRING_ELT(id_col, i, NA_STRING);
            SET_STRING_ELT(ref_col, i, NA_STRING);
            SET_STRING_ELT(alt_col, i, NA_STRING);
            REAL(qual_col)[i] = NA_REAL;
            SET_STRING_ELT(filter_col, i, NA_STRING);
            INTEGER(nallele_col)[i] = NA_INTEGER;
            if (csq_hrec) SET_VECTOR_ELT(csq_col, i, R_NilValue);
            if (ann_hrec) SET_VECTOR_ELT(ann_col, i, R_NilValue);
            INTEGER(index_col)[i] = NA_INTEGER;
            continue;
        }
        bcf_unpack(rec, BCF_UN_STR|BCF_UN_INFO|BCF_UN_FLT);
        // chrom
        const char *chrom = bcf_hdr_id2name(hdr, rec->rid);
        SET_STRING_ELT(chrom_col, i, chrom?mkChar(chrom):NA_STRING);
        // pos (1-based)
        INTEGER(pos_col)[i] = rec->pos + 1;
        // id
        if (rec->d.id && strcmp(rec->d.id, ".")!=0) SET_STRING_ELT(id_col, i, mkChar(rec->d.id)); else SET_STRING_ELT(id_col, i, NA_STRING);
        // ref & alt
        if (rec->d.allele && rec->n_allele>0) {
            SET_STRING_ELT(ref_col, i, mkChar(rec->d.allele[0]));
            // build alt string
            kalt.l = 0; if (rec->n_allele>1) {
                for (int a=1;a<rec->n_allele;a++) {
                    if (a>1) kputc(',', &kalt);
                    kputs(rec->d.allele[a], &kalt);
                }
                SET_STRING_ELT(alt_col, i, mkChar(kalt.s));
            } else {
                SET_STRING_ELT(alt_col, i, mkChar("."));
            }
        } else {
            SET_STRING_ELT(ref_col, i, NA_STRING);
            SET_STRING_ELT(alt_col, i, NA_STRING);
        }
        INTEGER(nallele_col)[i] = rec->n_allele;
        INTEGER(index_col)[i] = idx_var + 1; // 1-based index of variant in VBI
        // qual
        REAL(qual_col)[i] = bcf_float_is_missing(rec->qual) ? NA_REAL : rec->qual;
        // filters
        kflt.l = 0;
        if (rec->d.n_flt == 0) {
            kputs("PASS", &kflt);
        } else {
            for (int f=0; f<rec->d.n_flt; ++f) {
                if (f) kputc(';', &kflt);
                const char *flt = bcf_hdr_int2id(hdr, BCF_DT_ID, rec->d.flt[f]);
                if (flt) kputs(flt, &kflt);
            }
        }
        SET_STRING_ELT(filter_col, i, kflt.s?mkChar(kflt.s):mkChar(""));

        // CSQ parsing
        if (csq_hrec) {
            bcf_info_t *info = bcf_get_info(hdr, rec, "CSQ");
            if (!info) {
                SET_VECTOR_ELT(csq_col, i, R_NilValue);
            } else {
                // get string value
                char *csq_str = NULL; int ns = 0;
                if (bcf_get_info_string(hdr, rec, "CSQ", &csq_str, &ns) > 0) {
                    // split by ',' entries
                    int ntrans=1; for (char *p=csq_str; *p; ++p) if (*p==',') ntrans++;
                    // allocate data.frame columns per field
                    SEXP cols = PROTECT(allocVector(VECSXP, n_csq_fields));
                    for (int c=0;c<n_csq_fields;c++) SET_VECTOR_ELT(cols, c, allocVector(STRSXP, ntrans));
                    // parse
                    int t=0; char *entry = csq_str; char *p=csq_str;
                    while (1) {
                        if (*p==',' || *p=='\0') {
                            char save = *p; *p='\0';
                            // parse one transcript entry by '|'
                            int fld=0; char *seg=entry; char *q=entry;
                            while (1) {
                                if (*q=='|' || *q=='\0') {
                                    char save2=*q; *q='\0';
                                    if (fld<n_csq_fields) {
                                        SEXP col = VECTOR_ELT(cols,fld);
                                        SET_STRING_ELT(col, t, seg && *seg? mkChar(seg): NA_STRING);
                                    }
                                    *q = save2;
                                    fld++; seg = q+1;
                                }
                                if (*q=='\0') break; q++;
                            }
                            *p = save;
                            t++;
                            if (save=='\0') break; // done
                            entry = p+1;
                        }
                        if (*p=='\0') break; else p++;
                    }
                    // set names
                    SEXP fn = PROTECT(allocVector(STRSXP, n_csq_fields));
                    for (int c=0;c<n_csq_fields;c++) SET_STRING_ELT(fn,c,mkChar(csq_fields[c]));
                    setAttrib(cols, R_NamesSymbol, fn);
                    // row.names
                    SEXP rn = PROTECT(allocVector(INTSXP, ntrans)); for (int r=0;r<ntrans;r++) INTEGER(rn)[r]=r+1; setAttrib(cols, R_RowNamesSymbol, rn);
                    setAttrib(cols, R_ClassSymbol, mkString("data.frame"));
                    SET_VECTOR_ELT(csq_col, i, cols);
                    UNPROTECT(3); // cols names rn
                    free(csq_str);
                } else {
                    if (csq_str) free(csq_str);
                    SET_VECTOR_ELT(csq_col, i, R_NilValue);
                }
            }
        }
        // ANN parsing (similar but simpler: store raw entries list)
        if (ann_hrec) {
            bcf_info_t *info = bcf_get_info(hdr, rec, "ANN");
            if (!info) {
                SET_VECTOR_ELT(ann_col, i, R_NilValue);
            } else {
                char *ann_str=NULL; int ns=0;
                if (bcf_get_info_string(hdr, rec, "ANN", &ann_str, &ns) > 0) {
                    int ntrans=1; for (char *p=ann_str; *p; ++p) if (*p==',') ntrans++;
                    SEXP entries = PROTECT(allocVector(STRSXP, ntrans));
                    int t=0; char *entry=ann_str; char *p=ann_str;
                    while (1) {
                        if (*p==',' || *p=='\0') {
                            char save=*p; *p='\0';
                            SET_STRING_ELT(entries, t, entry && *entry?mkChar(entry):NA_STRING);
                            *p=save; t++; if (save=='\0') break; entry=p+1;
                        }
                        if (*p=='\0') break; else p++;
                    }
                    SET_VECTOR_ELT(ann_col, i, entries);
                    UNPROTECT(1);
                    free(ann_str);
                } else {
                    if (ann_str) free(ann_str);
                    SET_VECTOR_ELT(ann_col, i, R_NilValue);
                }
            }
        }
    }

    bcf_destroy(rec);
    if (kalt.s) free(kalt.s);
    if (kflt.s) free(kflt.s);
    bcf_hdr_destroy(hdr);
    hts_close(fp);

    // Build final data.frame
    int ncols = 9 + (csq_hrec?1:0) + (ann_hrec?1:0);
    SEXP df = PROTECT(allocVector(VECSXP, ncols)); nprotect++;
    SEXP names = PROTECT(allocVector(STRSXP, ncols)); nprotect++;
    int col=0;
    SET_VECTOR_ELT(df, col, chrom_col); SET_STRING_ELT(names, col++, mkChar("chrom"));
    SET_VECTOR_ELT(df, col, pos_col);   SET_STRING_ELT(names, col++, mkChar("pos"));
    SET_VECTOR_ELT(df, col, id_col);    SET_STRING_ELT(names, col++, mkChar("id"));
    SET_VECTOR_ELT(df, col, ref_col);   SET_STRING_ELT(names, col++, mkChar("ref"));
    SET_VECTOR_ELT(df, col, alt_col);   SET_STRING_ELT(names, col++, mkChar("alt"));
    SET_VECTOR_ELT(df, col, qual_col);  SET_STRING_ELT(names, col++, mkChar("qual"));
    SET_VECTOR_ELT(df, col, filter_col);SET_STRING_ELT(names, col++, mkChar("filter"));
    SET_VECTOR_ELT(df, col, nallele_col);SET_STRING_ELT(names, col++, mkChar("n_allele"));
    SET_VECTOR_ELT(df, col, index_col); SET_STRING_ELT(names, col++, mkChar("index"));
    if (csq_hrec) { SET_VECTOR_ELT(df, col, csq_col); SET_STRING_ELT(names, col++, mkChar("CSQ")); }
    if (ann_hrec) { SET_VECTOR_ELT(df, col, ann_col); SET_STRING_ELT(names, col++, mkChar("ANN")); }
    setAttrib(df, R_NamesSymbol, names);
    SEXP rn = PROTECT(allocVector(INTSXP, nfound)); nprotect++;
    for (int i=0;i<nfound;i++) INTEGER(rn)[i]=i+1;
    setAttrib(df, R_RowNamesSymbol, rn);
    setAttrib(df, R_ClassSymbol, mkString("data.frame"));

    // free temps
    if (indices) free(indices);
    if (csq_fields) { for (int i=0;i<n_csq_fields;i++) free(csq_fields[i]); free(csq_fields);}    
    if (ann_fields) { for (int i=0;i<n_ann_fields;i++) free(ann_fields[i]); free(ann_fields);}    

    UNPROTECT(nprotect); // all protected objects
    return df;
}

//' Query VCF/BCF by contiguous variant index range (nth to kth, 1-based, inclusive, no header)
// we are reopening the vcf everytime
// this might be slow for many small queries
// an alternative would be to memory mapped the vcf
// issue with ths is that we file may be  remote
SEXP RC_VBI_query_by_indices(SEXP vcf_path, SEXP idx_ptr, SEXP start_idx, SEXP end_idx, SEXP threads) {
    const char *vcf = CHAR(STRING_ELT(vcf_path, 0));
    vbi_index_t *idx = (vbi_index_t*) R_ExternalPtrAddr(idx_ptr);
    if (!idx) Rf_error("[VBI] Index pointer is NULL");
    int n_threads = asInteger(threads); (void)n_threads;
    int start = asInteger(start_idx) - 1;
    int end = asInteger(end_idx) - 1;
    if (start < 0) start = 0;
    if (end >= idx->num_marker) end = idx->num_marker - 1;
    if (end < start || start >= idx->num_marker) {
        // return empty data.frame
        SEXP df = PROTECT(allocVector(VECSXP, 0));
        SEXP names = PROTECT(allocVector(STRSXP, 0));
        SEXP rn = PROTECT(allocVector(INTSXP, 0));
        setAttrib(df, R_NamesSymbol, names);
        setAttrib(df, R_RowNamesSymbol, rn);
        setAttrib(df, R_ClassSymbol, mkString("data.frame"));
        UNPROTECT(3);
        return df;
    }
    int nfound = end - start + 1;

    htsFile *fp = hts_open(vcf, "r");
    if (!fp) Rf_error("Failed to open VCF/BCF: %s", vcf);
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (!hdr) { hts_close(fp); Rf_error("Failed to read VCF/BCF header: %s", vcf);}    

    bcf_hrec_t *csq_hrec = bcf_hdr_get_hrec(hdr, BCF_HL_INFO, "ID", "CSQ", NULL);
    bcf_hrec_t *ann_hrec = bcf_hdr_get_hrec(hdr, BCF_HL_INFO, "ID", "ANN", NULL);

    char **csq_fields = NULL; int n_csq_fields = 0;
    char **ann_fields = NULL; int n_ann_fields = 0;
    if (csq_hrec) {
        const char *desc = NULL; for (int i=0;i<csq_hrec->nkeys;i++) if (strcmp(csq_hrec->keys[i],"Description")==0) { desc = csq_hrec->vals[i]; break; }
        if (desc) {
            const char *fmt = strstr(desc, "Format:");
            if (fmt) { fmt += 7; while (*fmt==' '||*fmt=='\t') fmt++; char *tmp=strdup(fmt); size_t L=strlen(tmp); while (L>0 && (tmp[L-1]=='"'||tmp[L-1]==')'||tmp[L-1]==' ')) tmp[--L]='\0'; for (size_t i=0;i<L;i++) if (tmp[i]=='|') n_csq_fields++; n_csq_fields++; csq_fields=(char**)calloc(n_csq_fields,sizeof(char*)); int fld=0; char *p=tmp; char *startp=p; for(;*p;++p){ if(*p=='|'){ *p='\0'; csq_fields[fld++]=strdup(startp); startp=p+1; }} if(fld<n_csq_fields) csq_fields[fld++]=strdup(startp); free(tmp);} }
    }
    if (ann_hrec) {
        const char *desc = NULL; for (int i=0;i<ann_hrec->nkeys;i++) if (strcmp(ann_hrec->keys[i],"Description")==0) { desc = ann_hrec->vals[i]; break; }
        if (desc) {
            const char *fmt2 = strstr(desc, "Format:");
            const char *use = fmt2?fmt2:NULL; if (use){ use +=7; while(*use==' '||*use=='\t') use++; char *tmp=strdup(use); size_t L=strlen(tmp); while (L>0 && (tmp[L-1]=='"'||tmp[L-1]==')'||tmp[L-1]==' ')) tmp[--L]='\0'; for (size_t i=0;i<L;i++) if(tmp[i]=='|') n_ann_fields++; n_ann_fields++; ann_fields=(char**)calloc(n_ann_fields,sizeof(char*)); int fld=0; char *p=tmp; char *startp=p; for(;*p;++p){ if(*p=='|'){ *p='\0'; ann_fields[fld++]=strdup(startp); startp=p+1; }} if(fld<n_ann_fields) ann_fields[fld++]=strdup(startp); free(tmp);} }
    }

    SEXP chrom_col = PROTECT(allocVector(STRSXP, nfound));
    SEXP pos_col   = PROTECT(allocVector(INTSXP, nfound));
    SEXP id_col    = PROTECT(allocVector(STRSXP, nfound));
    SEXP ref_col   = PROTECT(allocVector(STRSXP, nfound));
    SEXP alt_col   = PROTECT(allocVector(STRSXP, nfound));
    SEXP qual_col  = PROTECT(allocVector(REALSXP, nfound));
    SEXP filter_col= PROTECT(allocVector(STRSXP, nfound));
    SEXP nallele_col=PROTECT(allocVector(INTSXP, nfound));
    SEXP index_col = PROTECT(allocVector(INTSXP, nfound));
    SEXP csq_col   = csq_hrec ? PROTECT(allocVector(VECSXP, nfound)) : R_NilValue; int prot_csx = csq_hrec?1:0;
    SEXP ann_col   = ann_hrec ? PROTECT(allocVector(VECSXP, nfound)) : R_NilValue; int prot_ann = ann_hrec?1:0;
    int nprotect = 9 + prot_csx + prot_ann;

    bcf1_t *rec = bcf_init();
    // seek to first
    int seek_ok = 0;
    if (fp->format.compression == bgzf) {
        BGZF *bg = (BGZF *)fp->fp.bgzf;
        seek_ok = (bgzf_seek(bg, idx->offsets[start], SEEK_SET) == 0);
    } else {
        hFILE *hf = (hFILE *)fp->fp.hfile;
        seek_ok = (hseek(hf, (off_t)idx->offsets[start], SEEK_SET) == 0);
    }
    if (!seek_ok) {
        bcf_destroy(rec); bcf_hdr_destroy(hdr); hts_close(fp);
        Rf_error("Failed to seek to first record");
    }
    kstring_t kalt = {0,0,0};
    kstring_t kflt = {0,0,0};
    for (int i=0;i<nfound;i++) {
        if (bcf_read(fp, hdr, rec) < 0) {
            SET_STRING_ELT(chrom_col, i, NA_STRING);
            INTEGER(pos_col)[i] = NA_INTEGER;
            SET_STRING_ELT(id_col, i, NA_STRING);
            SET_STRING_ELT(ref_col, i, NA_STRING);
            SET_STRING_ELT(alt_col, i, NA_STRING);
            REAL(qual_col)[i] = NA_REAL; SET_STRING_ELT(filter_col, i, NA_STRING); INTEGER(nallele_col)[i]=NA_INTEGER;
            if (csq_hrec) SET_VECTOR_ELT(csq_col,i,R_NilValue);
            if (ann_hrec) SET_VECTOR_ELT(ann_col,i,R_NilValue);
            INTEGER(index_col)[i] = NA_INTEGER;
            continue;
        }
        bcf_unpack(rec, BCF_UN_STR|BCF_UN_INFO|BCF_UN_FLT);
        const char *chrom = bcf_hdr_id2name(hdr, rec->rid);
        SET_STRING_ELT(chrom_col, i, chrom?mkChar(chrom):NA_STRING);
        INTEGER(pos_col)[i] = rec->pos + 1;
        if (rec->d.id && strcmp(rec->d.id, ".")!=0) SET_STRING_ELT(id_col,i,mkChar(rec->d.id)); else SET_STRING_ELT(id_col,i,NA_STRING);
        if (rec->d.allele && rec->n_allele>0) {
            SET_STRING_ELT(ref_col,i,mkChar(rec->d.allele[0]));
            kalt.l=0; if(rec->n_allele>1){ for(int a=1;a<rec->n_allele;a++){ if(a>1) kputc(',',&kalt); kputs(rec->d.allele[a],&kalt);} SET_STRING_ELT(alt_col,i,mkChar(kalt.s)); } else SET_STRING_ELT(alt_col,i,mkChar("."));
        } else { SET_STRING_ELT(ref_col,i,NA_STRING); SET_STRING_ELT(alt_col,i,NA_STRING); }
        INTEGER(nallele_col)[i] = rec->n_allele;
        INTEGER(index_col)[i] = (start + i) + 1; // 1-based global index
        REAL(qual_col)[i] = bcf_float_is_missing(rec->qual) ? NA_REAL : rec->qual;
        kflt.l=0; if(rec->d.n_flt==0){ kputs("PASS",&kflt);} else { for(int f=0;f<rec->d.n_flt;f++){ if(f) kputc(';',&kflt); const char *flt=bcf_hdr_int2id(hdr, BCF_DT_ID, rec->d.flt[f]); if(flt) kputs(flt,&kflt);} }
        SET_STRING_ELT(filter_col,i, kflt.s?mkChar(kflt.s):mkChar(""));
        if (csq_hrec) { bcf_info_t *info = bcf_get_info(hdr, rec, "CSQ"); if(!info){ SET_VECTOR_ELT(csq_col,i,R_NilValue);} else { char *csq_str=NULL; int ns=0; if(bcf_get_info_string(hdr,rec,"CSQ",&csq_str,&ns)>0){ int ntrans=1; for(char *p=csq_str;*p;++p) if(*p==',') ntrans++; SEXP cols=PROTECT(allocVector(VECSXP,n_csq_fields)); for(int c=0;c<n_csq_fields;c++) SET_VECTOR_ELT(cols,c,allocVector(STRSXP,ntrans)); int t=0; char *entry=csq_str; char *p=csq_str; while(1){ if(*p==',' || *p=='\0'){ char save=*p; *p='\0'; int fld=0; char *seg=entry; char *q=entry; while(1){ if(*q=='|' || *q=='\0'){ char save2=*q; *q='\0'; if(fld<n_csq_fields){ SEXP colv=VECTOR_ELT(cols,fld); SET_STRING_ELT(colv,t, seg && *seg?mkChar(seg):NA_STRING);} *q=save2; fld++; seg=q+1;} if(*q=='\0') break; q++; } *p=save; t++; if(save=='\0') break; entry=p+1; } if(*p=='\0') break; else p++; } SEXP fn=PROTECT(allocVector(STRSXP,n_csq_fields)); for(int c=0;c<n_csq_fields;c++) SET_STRING_ELT(fn,c,mkChar(csq_fields[c])); setAttrib(cols,R_NamesSymbol,fn); SEXP rn=PROTECT(allocVector(INTSXP,ntrans)); for(int r=0;r<ntrans;r++) INTEGER(rn)[r]=r+1; setAttrib(cols,R_RowNamesSymbol,rn); setAttrib(cols,R_ClassSymbol,mkString("data.frame")); SET_VECTOR_ELT(csq_col,i,cols); UNPROTECT(3); free(csq_str);} else { if(csq_str) free(csq_str); SET_VECTOR_ELT(csq_col,i,R_NilValue);} }}
        if (ann_hrec) { bcf_info_t *info = bcf_get_info(hdr, rec, "ANN"); if(!info){ SET_VECTOR_ELT(ann_col,i,R_NilValue);} else { char *ann_str=NULL; int ns=0; if(bcf_get_info_string(hdr,rec,"ANN",&ann_str,&ns)>0){ int ntrans=1; for(char *p=ann_str;*p;++p) if(*p==',') ntrans++; SEXP entries=PROTECT(allocVector(STRSXP,ntrans)); int t=0; char *entry=ann_str; char *p=ann_str; while(1){ if(*p==',' || *p=='\0'){ char save=*p; *p='\0'; SET_STRING_ELT(entries,t, entry && *entry?mkChar(entry):NA_STRING); *p=save; t++; if(save=='\0') break; entry=p+1;} if(*p=='\0') break; else p++; } SET_VECTOR_ELT(ann_col,i,entries); UNPROTECT(1); free(ann_str);} else { if(ann_str) free(ann_str); SET_VECTOR_ELT(ann_col,i,R_NilValue);} }}
    }
    bcf_destroy(rec);
    if (kalt.s) free(kalt.s); if (kflt.s) free(kflt.s);
    bcf_hdr_destroy(hdr); hts_close(fp);

    int ncols = 9 + (csq_hrec?1:0) + (ann_hrec?1:0);
    SEXP df = PROTECT(allocVector(VECSXP, ncols)); nprotect++;
    SEXP names = PROTECT(allocVector(STRSXP, ncols)); nprotect++;
    int col=0; SET_VECTOR_ELT(df,col,chrom_col); SET_STRING_ELT(names,col++,mkChar("chrom"));
    SET_VECTOR_ELT(df,col,pos_col); SET_STRING_ELT(names,col++,mkChar("pos"));
    SET_VECTOR_ELT(df,col,id_col); SET_STRING_ELT(names,col++,mkChar("id"));
    SET_VECTOR_ELT(df,col,ref_col); SET_STRING_ELT(names,col++,mkChar("ref"));
    SET_VECTOR_ELT(df,col,alt_col); SET_STRING_ELT(names,col++,mkChar("alt"));
    SET_VECTOR_ELT(df,col,qual_col); SET_STRING_ELT(names,col++,mkChar("qual"));
    SET_VECTOR_ELT(df,col,filter_col); SET_STRING_ELT(names,col++,mkChar("filter"));
    SET_VECTOR_ELT(df,col,nallele_col); SET_STRING_ELT(names,col++,mkChar("n_allele"));
    SET_VECTOR_ELT(df,col,index_col); SET_STRING_ELT(names,col++,mkChar("index"));
    if (csq_hrec) { SET_VECTOR_ELT(df,col,csq_col); SET_STRING_ELT(names,col++,mkChar("CSQ")); }
    if (ann_hrec) { SET_VECTOR_ELT(df,col,ann_col); SET_STRING_ELT(names,col++,mkChar("ANN")); }
    setAttrib(df, R_NamesSymbol, names);
    SEXP rn = PROTECT(allocVector(INTSXP, nfound)); nprotect++;
    for (int i=0;i<nfound;i++) INTEGER(rn)[i]=i+1;
    setAttrib(df, R_RowNamesSymbol, rn);
    setAttrib(df, R_ClassSymbol, mkString("data.frame"));

    if (csq_fields) { for (int i=0;i<n_csq_fields;i++) free(csq_fields[i]); free(csq_fields);}    
    if (ann_fields) { for (int i=0;i<n_ann_fields;i++) free(ann_fields[i]); free(ann_fields);}    

    UNPROTECT(nprotect);
    return df;
}


SEXP RC_VBI_print_index(SEXP idx_ptr, SEXP n) {
    vbi_index_t *idx = (vbi_index_t*) R_ExternalPtrAddr(idx_ptr);
    if (!idx) Rf_error("[VBI] Index pointer is NULL");
    int nlines = asInteger(n);
    if (nlines > idx->num_marker) nlines = idx->num_marker;
    Rprintf("[VBI] num_marker: %ld\n", (long)idx->num_marker);
    for (int i = 0; i < nlines; ++i) {
        int chrom_id = idx->chrom_ids[i];
        const char *chrom_name = idx->chrom_names ? idx->chrom_names[chrom_id] : "?";
        Rprintf("[VBI] %d: %s\t%ld offset=%ld\n", i+1, chrom_name, (long)idx->positions[i], (long)idx->offsets[i]);
    }
    return R_NilValue;
}


// R-callable: load VBI index and return external pointer
SEXP RC_VBI_load_index(SEXP vbi_path) {
    const char *path = CHAR(STRING_ELT(vbi_path, 0));
    char local_path[4096];
    if (strstr(path, "://")) {
        snprintf(local_path, sizeof(local_path), "/tmp/vbi_%u.vbi", (unsigned)rand());
        if (!download_file(path, local_path)) {
            Rf_error("Failed to download VBI index: %s", path);
        }
        path = local_path;
        }
    vbi_index_t *idx = vbi_index_load(path);
    if (!idx) Rf_error("[VBI] Failed to load index: %s", path);
    SEXP ptr = PROTECT(R_MakeExternalPtr(idx, R_NilValue, R_NilValue));
    R_RegisterCFinalizerEx(ptr, vbi_index_finalizer, 1);
    UNPROTECT(1);
    return ptr;
}


// R-callable wrapper for cgranges region query
SEXP RC_VBI_query_region_cgranges(SEXP vcf_path, SEXP idx_ptr, SEXP region_str) {
    const char *vcf = CHAR(STRING_ELT(vcf_path, 0));
    vbi_index_t *idx = (vbi_index_t*) R_ExternalPtrAddr(idx_ptr);
    if (!idx) Rf_error("[VBI] Index pointer is NULL");
    const char *reg = CHAR(STRING_ELT(region_str, 0));
    int nfound = 0;
    int *hits = vbi_index_query_region_cgranges(idx, reg, &nfound);
    if (!hits || nfound == 0) {
        SEXP df = PROTECT(allocVector(VECSXP, 0));
        SEXP names = PROTECT(allocVector(STRSXP, 0));
        SEXP rn = PROTECT(allocVector(INTSXP, 0));
        setAttrib(df, R_NamesSymbol, names);
        setAttrib(df, R_RowNamesSymbol, rn);
        setAttrib(df, R_ClassSymbol, mkString("data.frame"));
        UNPROTECT(3);
        if (hits) free(hits);
        return df;
    }
    htsFile *fp = hts_open(vcf, "r");
    if (!fp) { free(hits); Rf_error("Failed to open VCF/BCF: %s", vcf);}    
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (!hdr) { hts_close(fp); free(hits); Rf_error("Failed to read VCF/BCF header: %s", vcf);}    
    bcf_hrec_t *csq_hrec = bcf_hdr_get_hrec(hdr, BCF_HL_INFO, "ID", "CSQ", NULL);
    bcf_hrec_t *ann_hrec = bcf_hdr_get_hrec(hdr, BCF_HL_INFO, "ID", "ANN", NULL);
    char **csq_fields=NULL; int n_csq_fields=0; char **ann_fields=NULL; int n_ann_fields=0; // consistent parsing (reuse code if needed later)
    if (csq_hrec) {
        const char *desc=NULL; for(int i=0;i<csq_hrec->nkeys;i++) if(strcmp(csq_hrec->keys[i],"Description")==0){ desc=csq_hrec->vals[i]; break; }
        if(desc){ const char *fmt=strstr(desc,"Format:"); if(fmt){ fmt+=7; while(*fmt==' '||*fmt=='\t') fmt++; char *tmp=strdup(fmt); size_t L=strlen(tmp); while(L>0 && (tmp[L-1]=='"'||tmp[L-1]==')'||tmp[L-1]==' ')) tmp[--L]='\0'; for(size_t i=0;i<L;i++) if(tmp[i]=='|') n_csq_fields++; n_csq_fields++; csq_fields=(char**)calloc(n_csq_fields,sizeof(char*)); int fld=0; char *p=tmp; char *startp=p; for(;*p;++p){ if(*p=='|'){ *p='\0'; csq_fields[fld++]=strdup(startp); startp=p+1; }} if(fld<n_csq_fields) csq_fields[fld++]=strdup(startp); free(tmp);} }
    }
    if (ann_hrec) {
        const char *desc=NULL; for(int i=0;i<ann_hrec->nkeys;i++) if(strcmp(ann_hrec->keys[i],"Description")==0){ desc=ann_hrec->vals[i]; break; }
        if(desc){ const char *fmt2=strstr(desc,"Format:"); const char *use=fmt2?fmt2:NULL; if(use){ use+=7; while(*use==' '||*use=='\t') use++; char *tmp=strdup(use); size_t L=strlen(tmp); while (L>0 && (tmp[L-1]=='"'||tmp[L-1]==')'||tmp[L-1]==' ')) tmp[--L]='\0'; for(size_t i=0;i<L;i++) if(tmp[i]=='|') n_ann_fields++; n_ann_fields++; ann_fields=(char**)calloc(n_ann_fields,sizeof(char*)); int fld=0; char *p=tmp; char *startp=p; for(;*p;++p){ if(*p=='|'){ *p='\0'; ann_fields[fld++]=strdup(startp); startp=p+1; }} if(fld<n_ann_fields) ann_fields[fld++]=strdup(startp); free(tmp);} }
    }
    SEXP chrom_col = PROTECT(allocVector(STRSXP, nfound));
    SEXP pos_col   = PROTECT(allocVector(INTSXP, nfound));
    SEXP id_col    = PROTECT(allocVector(STRSXP, nfound));
    SEXP ref_col   = PROTECT(allocVector(STRSXP, nfound));
    SEXP alt_col   = PROTECT(allocVector(STRSXP, nfound));
    SEXP qual_col  = PROTECT(allocVector(REALSXP, nfound));
    SEXP filter_col= PROTECT(allocVector(STRSXP, nfound));
    SEXP nallele_col=PROTECT(allocVector(INTSXP, nfound));
    SEXP index_col = PROTECT(allocVector(INTSXP, nfound));
    SEXP csq_col   = csq_hrec ? PROTECT(allocVector(VECSXP, nfound)) : R_NilValue; int prot_csx=csq_hrec?1:0;
    SEXP ann_col   = ann_hrec ? PROTECT(allocVector(VECSXP, nfound)) : R_NilValue; int prot_ann=ann_hrec?1:0;
    int nprotect = 9 + prot_csx + prot_ann;
    bcf1_t *rec = bcf_init();
    kstring_t kalt={0,0,0}; kstring_t kflt={0,0,0};
    for(int i=0;i<nfound;i++) {
        int idx_var = hits[i];
        int seek_ok=0; if (fp->format.compression==bgzf){ BGZF *bg=(BGZF*)fp->fp.bgzf; seek_ok=(bgzf_seek(bg, idx->offsets[idx_var], SEEK_SET)==0);} else { hFILE *hf=(hFILE*)fp->fp.hfile; seek_ok=(hseek(hf,(off_t)idx->offsets[idx_var],SEEK_SET)==0);} 
        if(!seek_ok || bcf_read(fp,hdr,rec)<0){ SET_STRING_ELT(chrom_col,i,NA_STRING); INTEGER(pos_col)[i]=NA_INTEGER; SET_STRING_ELT(id_col,i,NA_STRING); SET_STRING_ELT(ref_col,i,NA_STRING); SET_STRING_ELT(alt_col,i,NA_STRING); REAL(qual_col)[i]=NA_REAL; SET_STRING_ELT(filter_col,i,NA_STRING); INTEGER(nallele_col)[i]=NA_INTEGER; INTEGER(index_col)[i]=NA_INTEGER; if(csq_hrec) SET_VECTOR_ELT(csq_col,i,R_NilValue); if(ann_hrec) SET_VECTOR_ELT(ann_col,i,R_NilValue); continue; }
        bcf_unpack(rec, BCF_UN_STR|BCF_UN_INFO|BCF_UN_FLT);
        const char *chrom = bcf_hdr_id2name(hdr, rec->rid);
        SET_STRING_ELT(chrom_col, i, chrom?mkChar(chrom):NA_STRING);
        INTEGER(pos_col)[i] = rec->pos + 1;
        if (rec->d.id && strcmp(rec->d.id, ".")!=0) SET_STRING_ELT(id_col,i,mkChar(rec->d.id)); else SET_STRING_ELT(id_col,i,NA_STRING);
        if (rec->d.allele && rec->n_allele>0) {
            SET_STRING_ELT(ref_col,i,mkChar(rec->d.allele[0]));
            kalt.l=0; if(rec->n_allele>1){ for(int a=1;a<rec->n_allele;a++){ if(a>1) kputc(',',&kalt); kputs(rec->d.allele[a],&kalt);} SET_STRING_ELT(alt_col,i,mkChar(kalt.s)); } else SET_STRING_ELT(alt_col,i,mkChar("."));
        } else { SET_STRING_ELT(ref_col,i,NA_STRING); SET_STRING_ELT(alt_col,i,NA_STRING); }
        INTEGER(nallele_col)[i] = rec->n_allele; INTEGER(index_col)[i] = idx_var + 1;
        REAL(qual_col)[i] = bcf_float_is_missing(rec->qual) ? NA_REAL : rec->qual;
        kflt.l=0; if(rec->d.n_flt==0){ kputs("PASS",&kflt);} else { for(int f=0;f<rec->d.n_flt;f++){ if(f) kputc(';',&kflt); const char *flt=bcf_hdr_int2id(hdr, BCF_DT_ID, rec->d.flt[f]); if(flt) kputs(flt,&kflt);} }
        SET_STRING_ELT(filter_col,i, kflt.s?mkChar(kflt.s):mkChar(""));
        if (csq_hrec) { bcf_info_t *info=bcf_get_info(hdr, rec, "CSQ"); if(!info){ SET_VECTOR_ELT(csq_col,i,R_NilValue);} else { char *csq_str=NULL; int ns=0; if(bcf_get_info_string(hdr,rec,"CSQ",&csq_str,&ns)>0){ int ntrans=1; for(char *p=csq_str;*p;++p) if(*p==',') ntrans++; SEXP cols=PROTECT(allocVector(VECSXP,n_csq_fields)); for(int c=0;c<n_csq_fields;c++) SET_VECTOR_ELT(cols,c,allocVector(STRSXP,ntrans)); int t=0; char *entry=csq_str; char *p=csq_str; while(1){ if(*p==',' || *p=='\0'){ char save=*p; *p='\0'; int fld=0; char *seg=entry; char *q=entry; while(1){ if(*q=='|' || *q=='\0'){ char save2=*q; *q='\0'; if(fld<n_csq_fields){ SEXP colv=VECTOR_ELT(cols,fld); SET_STRING_ELT(colv,t, seg && *seg?mkChar(seg):NA_STRING);} *q=save2; fld++; seg=q+1;} if(*q=='\0') break; q++; } *p=save; t++; if(save=='\0') break; entry=p+1; } if(*p=='\0') break; else p++; } SEXP fn=PROTECT(allocVector(STRSXP,n_csq_fields)); for(int c=0;c<n_csq_fields;c++) SET_STRING_ELT(fn,c,mkChar(csq_fields[c])); setAttrib(cols,R_NamesSymbol,fn); SEXP rn=PROTECT(allocVector(INTSXP,ntrans)); for(int r=0;r<ntrans;r++) INTEGER(rn)[r]=r+1; setAttrib(cols,R_RowNamesSymbol,rn); setAttrib(cols,R_ClassSymbol,mkString("data.frame")); SET_VECTOR_ELT(csq_col,i,cols); UNPROTECT(3); free(csq_str);} else { if(csq_str) free(csq_str); SET_VECTOR_ELT(csq_col,i,R_NilValue);} }}
        if (ann_hrec) { bcf_info_t *info = bcf_get_info(hdr, rec, "ANN"); if(!info){ SET_VECTOR_ELT(ann_col,i,R_NilValue);} else { char *ann_str=NULL; int ns=0; if(bcf_get_info_string(hdr,rec,"ANN",&ann_str,&ns)>0){ int ntrans=1; for(char *p=ann_str;*p;++p) if(*p==',') ntrans++; SEXP entries=PROTECT(allocVector(STRSXP,ntrans)); int t=0; char *entry=ann_str; char *p=ann_str; while(1){ if(*p==',' || *p=='\0'){ char save=*p; *p='\0'; SET_STRING_ELT(entries,t, entry && *entry?mkChar(entry):NA_STRING); *p=save; t++; if(save=='\0') break; entry=p+1;} if(*p=='\0') break; else p++; } SET_VECTOR_ELT(ann_col,i,entries); UNPROTECT(1); free(ann_str);} else { if(ann_str) free(ann_str); SET_VECTOR_ELT(ann_col,i,R_NilValue);} }}
    }
    bcf_destroy(rec); if(kalt.s) free(kalt.s); if(kflt.s) free(kflt.s); bcf_hdr_destroy(hdr); hts_close(fp); free(hits);
    int ncols = 9 + (csq_hrec?1:0) + (ann_hrec?1:0);
    SEXP df = PROTECT(allocVector(VECSXP, ncols)); nprotect++;
    SEXP names = PROTECT(allocVector(STRSXP, ncols)); nprotect++;
    int col=0; SET_VECTOR_ELT(df,col,chrom_col); SET_STRING_ELT(names,col++,mkChar("chrom"));
    SET_VECTOR_ELT(df,col,pos_col); SET_STRING_ELT(names,col++,mkChar("pos"));
    SET_VECTOR_ELT(df,col,id_col); SET_STRING_ELT(names,col++,mkChar("id"));
    SET_VECTOR_ELT(df,col,ref_col); SET_STRING_ELT(names,col++,mkChar("ref"));
    SET_VECTOR_ELT(df,col,alt_col); SET_STRING_ELT(names,col++,mkChar("alt"));
    SET_VECTOR_ELT(df,col,qual_col); SET_STRING_ELT(names,col++,mkChar("qual"));
    SET_VECTOR_ELT(df,col,filter_col); SET_STRING_ELT(names,col++,mkChar("filter"));
    SET_VECTOR_ELT(df,col,nallele_col); SET_STRING_ELT(names,col++,mkChar("n_allele"));
    SET_VECTOR_ELT(df,col,index_col); SET_STRING_ELT(names,col++,mkChar("index"));
    if(csq_hrec){ SET_VECTOR_ELT(df,col,csq_col); SET_STRING_ELT(names,col++,mkChar("CSQ")); }
    if(ann_hrec){ SET_VECTOR_ELT(df,col,ann_col); SET_STRING_ELT(names,col++,mkChar("ANN")); }
    setAttrib(df,R_NamesSymbol,names); SEXP rn = PROTECT(allocVector(INTSXP,nfound)); nprotect++; for(int i=0;i<nfound;i++) INTEGER(rn)[i]=i+1; setAttrib(df,R_RowNamesSymbol,rn); setAttrib(df,R_ClassSymbol,mkString("data.frame"));
    if(csq_fields){ for(int i=0;i<n_csq_fields;i++) free(csq_fields[i]); free(csq_fields);} if(ann_fields){ for(int i=0;i<n_ann_fields;i++) free(ann_fields[i]); free(ann_fields);} 
    UNPROTECT(nprotect);
    return df;
}


// VBIExtractRanges: extract chrom/start/end/label from VBI index
SEXP RC_VBI_extract_ranges(SEXP idx_ptr, SEXP n) {
    vbi_index_t *idx = (vbi_index_t*) R_ExternalPtrAddr(idx_ptr);
    if (!idx) Rf_error("[VBI] Index pointer is NULL");
    int nvar = idx->num_marker;
    int nout = n == R_NilValue || INTEGER(n)[0] == NA_INTEGER ? nvar : INTEGER(n)[0];
    if (nout > nvar) nout = nvar;
    SEXP chroms = PROTECT(allocVector(STRSXP, nout));
    SEXP starts = PROTECT(allocVector(INTSXP, nout));
    SEXP ends   = PROTECT(allocVector(INTSXP, nout));
    SEXP labels = PROTECT(allocVector(INTSXP, nout));
    for (int i = 0; i < nout; ++i) {
        SET_STRING_ELT(chroms, i, mkChar(idx->chrom_names[idx->chrom_ids[i]]));
        INTEGER(starts)[i] = idx->positions[i];
        INTEGER(ends)[i]   = idx->positions[i];
        INTEGER(labels)[i] = i;
    }
    SEXP df = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(df, 0, chroms);
    SET_VECTOR_ELT(df, 1, starts);
    SET_VECTOR_ELT(df, 2, ends);
    SET_VECTOR_ELT(df, 3, labels);
    SEXP names = PROTECT(allocVector(STRSXP, 4));
    SET_STRING_ELT(names, 0, mkChar("chrom"));
    SET_STRING_ELT(names, 1, mkChar("start"));
    SET_STRING_ELT(names, 2, mkChar("end"));
    SET_STRING_ELT(names, 3, mkChar("label"));
    setAttrib(df, R_NamesSymbol, names);
    UNPROTECT(6);
    return df;
}

// Minimal cgranges R binding
typedef struct {
    cgranges_t *cr;
} cgranges_ptr_t;

SEXP RC_cgranges_create() {
    cgranges_ptr_t *ptr = (cgranges_ptr_t*) R_Calloc(1, cgranges_ptr_t);
    ptr->cr = cr_init();
    SEXP ext = PROTECT(R_MakeExternalPtr(ptr, R_NilValue, R_NilValue));
    R_RegisterCFinalizerEx(ext, (R_CFinalizer_t)RC_cgranges_destroy, 1);
    UNPROTECT(1);
    return ext;
}

SEXP RC_cgranges_add(SEXP cr_ptr, SEXP chrom, SEXP start, SEXP end, SEXP label) {
    cgranges_ptr_t *ptr = (cgranges_ptr_t*) R_ExternalPtrAddr(cr_ptr);
    if (!ptr || !ptr->cr) Rf_error("[cgranges] Null pointer");
    cr_add(ptr->cr, CHAR(STRING_ELT(chrom, 0)), INTEGER(start)[0], INTEGER(end)[0], INTEGER(label)[0]);
    return R_NilValue;
}

SEXP RC_cgranges_index(SEXP cr_ptr) {
    cgranges_ptr_t *ptr = (cgranges_ptr_t*) R_ExternalPtrAddr(cr_ptr);
    if (!ptr || !ptr->cr) Rf_error("[cgranges] Null pointer");
    cr_index(ptr->cr);
    return R_NilValue;
}

SEXP RC_cgranges_overlap(SEXP cr_ptr, SEXP chrom, SEXP start, SEXP end) {
    cgranges_ptr_t *ptr = (cgranges_ptr_t*) R_ExternalPtrAddr(cr_ptr);
    if (!ptr || !ptr->cr) Rf_error("[cgranges] Null pointer");
    int n = LENGTH(chrom);
    if (LENGTH(start) != n || LENGTH(end) != n) Rf_error("chrom, start, end must have same length");
    SEXP out = PROTECT(allocVector(VECSXP, n));
    for (int i = 0; i < n; ++i) {
        int64_t *b = 0, max_b = 0;
        int nhit = cr_overlap(ptr->cr, CHAR(STRING_ELT(chrom, i)), INTEGER(start)[i], INTEGER(end)[i], &b, &max_b);
        SEXP res = PROTECT(allocVector(INTSXP, nhit));
        for (int j = 0; j < nhit; ++j) INTEGER(res)[j] = (int)b[j] + 1; // 1-based for R
        if (b) free(b);
        SET_VECTOR_ELT(out, i, res);
        UNPROTECT(1);
    }
    UNPROTECT(1);
    return out;
}

SEXP RC_cgranges_destroy(SEXP cr_ptr) {
    cgranges_ptr_t *ptr = (cgranges_ptr_t*) R_ExternalPtrAddr(cr_ptr);
    if (ptr && ptr->cr) {
        cr_destroy(ptr->cr);
        ptr->cr = NULL;
        R_Free(ptr);
    }
    R_ClearExternalPtr(cr_ptr);
    return R_NilValue;
}

// Helper: extract interval info from cgranges by 0-based index
static void cr_extract_one(cgranges_t *cr, int idx, char **chrom, int *start, int *end, int *label) {
    if (!cr || !chrom || !start || !end || !label) return;
    if (idx < 0 || idx >= cr->n_r) {
        *chrom = NULL;
        *start = *end = *label = NA_INTEGER;
        return;
    }
    const char *cname = cr_chrom(cr, idx);
    *chrom = (char*)cname;
    *start = cr_start(cr, idx);
    *end = cr_end(cr, idx);
    *label = cr_label(cr, idx);
}

// Extract intervals by 1-based indices (R)
SEXP RC_cgranges_extract_by_index(SEXP cr_ptr, SEXP indices) {
    cgranges_ptr_t *ptr = (cgranges_ptr_t*) R_ExternalPtrAddr(cr_ptr);
    if (!ptr || !ptr->cr) Rf_error("[cgranges] Null pointer");
    int n = LENGTH(indices);
    SEXP chroms = PROTECT(allocVector(STRSXP, n));
    SEXP starts = PROTECT(allocVector(INTSXP, n));
    SEXP ends   = PROTECT(allocVector(INTSXP, n));
    SEXP labels = PROTECT(allocVector(INTSXP, n));
    for (int i = 0; i < n; ++i) {
        int idx = INTEGER(indices)[i] - 1; // 1-based to 0-based
        char *chrom = NULL; int start = NA_INTEGER, end = NA_INTEGER, label = NA_INTEGER;
        cr_extract_one(ptr->cr, idx, &chrom, &start, &end, &label);
        if (chrom)
            SET_STRING_ELT(chroms, i, mkChar(chrom));
        else
            SET_STRING_ELT(chroms, i, NA_STRING);
        INTEGER(starts)[i] = start;
        INTEGER(ends)[i] = end;
        INTEGER(labels)[i] = (label == NA_INTEGER) ? NA_INTEGER : label + 1; // 1-based for R
    }
    SEXP df = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(df, 0, chroms);
    SET_VECTOR_ELT(df, 1, starts);
    SET_VECTOR_ELT(df, 2, ends);
    SET_VECTOR_ELT(df, 3, labels);
    SEXP names = PROTECT(allocVector(STRSXP, 4));
    SET_STRING_ELT(names, 0, mkChar("chrom"));
    SET_STRING_ELT(names, 1, mkChar("start"));
    SET_STRING_ELT(names, 2, mkChar("end"));
    SET_STRING_ELT(names, 3, mkChar("label"));
    setAttrib(df, R_NamesSymbol, names);
    UNPROTECT(6);
    return df;
}


// Return the memory usage in bytes of a vbi_index_t structure
SEXP RC_VBI_index_memory_usage(SEXP extPtr) {
    vbi_index_t *idx = (vbi_index_t*) R_ExternalPtrAddr(extPtr);
    if (!idx) return ScalarReal(NA_REAL);
    size_t vbi_bytes = sizeof(vbi_index_t);
    if (idx->chrom_ids) vbi_bytes += idx->num_marker * sizeof(int32_t);
    if (idx->positions) vbi_bytes += idx->num_marker * sizeof(int64_t);
    if (idx->offsets) vbi_bytes += idx->num_marker * sizeof(int64_t);
    if (idx->chrom_names) {
        for (int i = 0; i < idx->n_chroms; ++i) {
            if (idx->chrom_names[i]) vbi_bytes += strlen(idx->chrom_names[i]) + 1;
        }
        vbi_bytes += idx->n_chroms * sizeof(char*);
    }
    size_t cr_bytes = 0;
    if (idx->cr) {
        cgranges_t *cr = (cgranges_t*) idx->cr;
        cr_bytes += sizeof(cgranges_t);
        if (cr->r) cr_bytes += cr->m_r * sizeof(cr_intv_t);
        if (cr->ctg) cr_bytes += cr->m_ctg * sizeof(cr_ctg_t);
        if (cr->hc) {/* skip, unknown size */}
        // Add up contig name strings
        if (cr->ctg) {
            for (int i = 0; i < cr->n_ctg; ++i) {
                if (cr->ctg[i].name) cr_bytes += strlen(cr->ctg[i].name) + 1;
            }
        }
    }
    SEXP out = PROTECT(allocVector(VECSXP, 2));
    SEXP nms = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(nms, 0, mkChar("vbi_index_t_bytes"));
    SET_STRING_ELT(nms, 1, mkChar("cgranges_t_bytes"));
    SET_VECTOR_ELT(out, 0, ScalarReal((double)vbi_bytes));
    SET_VECTOR_ELT(out, 1, ScalarReal((double)cr_bytes));
    setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(2);
    return out;
}
