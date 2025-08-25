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
#include "cgranges.h"
// Forward declaration for finalizer (after SEXP is defined)
SEXP RC_cgranges_destroy(SEXP cr_ptr);
// Forward declarations from vbi_index.c and vbi_index_capi.h
#include "vbi_index_capi.h"
int do_index(const char *infile, const char *outfile, int n_threads);

// Minimal cr_chrom implementation for cgranges_t (do not modify cgranges library files)
static inline const char *cr_chrom(const cgranges_t *cr, int64_t i) {
    int32_t ctg_id = (int32_t)(cr->r[i].x >> 32);
    return cr->ctg[ctg_id].name;
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
static int download_file(const char *url, const char *dest) {
    SEXP call = PROTECT(lang4(install("download.file"), mkString(url), mkString(dest), mkString("auto")));
    SEXP res = PROTECT(R_tryEval(call, R_GlobalEnv, NULL));
    int status = asInteger(res);
    UNPROTECT(2);
    return status == 0;
}

//' Query VCF/BCF by region using VBI index (returns only records, no header)
SEXP RC_VBI_query_range(SEXP vcf_path, SEXP idx_ptr, SEXP region, SEXP threads) {
    const char *vcf = CHAR(STRING_ELT(vcf_path, 0));
    vbi_index_t *idx = (vbi_index_t*) R_ExternalPtrAddr(idx_ptr);
    if (!idx) Rf_error("[VBI] Index pointer is NULL");
    const char *reg = CHAR(STRING_ELT(region, 0));
    int n_threads = asInteger(threads);
    SEXP lines = R_NilValue;
    PROTECT_INDEX idx_prot;
    PROTECT_WITH_INDEX(lines = R_NilValue, &idx_prot);
    int nfound = 0;
    int *indices = vbi_index_query_region(idx, reg, &nfound);
    if (!indices || nfound == 0) {
        lines = allocVector(STRSXP, 0);
        REPROTECT(lines, idx_prot);
        UNPROTECT(1);
        return lines;
    }
    htsFile *fp = hts_open(vcf, "r");
    if (!fp) {
        free(indices);
        UNPROTECT(1);
        Rf_error("Failed to open VCF/BCF: %s", vcf);
    }
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (!hdr) {
        hts_close(fp);
        free(indices);
        UNPROTECT(1);
        Rf_error("Failed to read VCF/BCF header: %s", vcf);
    }
    lines = allocVector(STRSXP, nfound);
    REPROTECT(lines, idx_prot);
    bcf1_t *rec = bcf_init();
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
        if (!seek_ok) {
            SET_STRING_ELT(lines, i, mkChar(""));
            continue;
        }
        int r = bcf_read(fp, hdr, rec);
        if (r < 0) {
            SET_STRING_ELT(lines, i, mkChar(""));
            continue;
        }
        bcf_unpack(rec, BCF_UN_STR);
        kstring_t str = {0,0,0};
        vcf_format1(hdr, rec, &str);
        SET_STRING_ELT(lines, i, mkChar(str.s));
        free(str.s);
    }
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    hts_close(fp);
    free(indices);
    UNPROTECT(1);
    return lines;
}

//' Query VCF/BCF by contiguous variant index range (nth to kth, 1-based, inclusive, no header)
SEXP RC_VBI_query_by_indices(SEXP vcf_path, SEXP idx_ptr, SEXP start_idx, SEXP end_idx, SEXP threads) {
    const char *vcf = CHAR(STRING_ELT(vcf_path, 0));
    vbi_index_t *idx = (vbi_index_t*) R_ExternalPtrAddr(idx_ptr);
    if (!idx) Rf_error("[VBI] Index pointer is NULL");
    int n_threads = asInteger(threads);
    int start = asInteger(start_idx) - 1;
    int end = asInteger(end_idx) - 1;
    SEXP lines = R_NilValue;
    PROTECT_INDEX idx_prot;
    PROTECT_WITH_INDEX(lines = R_NilValue, &idx_prot);
    int nfound = 0;
    nfound = end - start + 1;
    if (nfound <= 0) {
        lines = allocVector(STRSXP, 0);
        REPROTECT(lines, idx_prot);
        UNPROTECT(1);
        return lines;
    }
    htsFile *fp = hts_open(vcf, "r");
    if (!fp) {
        UNPROTECT(1);
        Rf_error("Failed to open VCF/BCF: %s", vcf);
    }
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (!hdr) {
        hts_close(fp);
        UNPROTECT(1);
        Rf_error("Failed to read VCF/BCF header: %s", vcf);
    }
    lines = allocVector(STRSXP, nfound);
    REPROTECT(lines, idx_prot);
    bcf1_t *rec = bcf_init();
    // Seek once to the first offset
    int seek_ok = 0;
    if (fp->format.compression == bgzf) {
        BGZF *bg = (BGZF *)fp->fp.bgzf;
        seek_ok = (bgzf_seek(bg, idx->offsets[start], SEEK_SET) == 0);
    } else {
        hFILE *hf = (hFILE *)fp->fp.hfile;
        seek_ok = (hseek(hf, (off_t)idx->offsets[start], SEEK_SET) == 0);
    }
    if (!seek_ok) {
        bcf_destroy(rec);
        bcf_hdr_destroy(hdr);
        hts_close(fp);
        UNPROTECT(1);
        Rf_error("Failed to seek to first record");
    }
    for (int i = 0; i < nfound; ++i) {
        int r = bcf_read(fp, hdr, rec);
        if (r < 0) {
            SET_STRING_ELT(lines, i, mkChar(""));
            continue;
        }
        bcf_unpack(rec, BCF_UN_STR);
        kstring_t str = {0,0,0};
        vcf_format1(hdr, rec, &str);
        SET_STRING_ELT(lines, i, mkChar(str.s));
        free(str.s);
    }
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    hts_close(fp);
    UNPROTECT(1);
    return lines;
}


SEXP RC_VBI_print_index(SEXP idx_ptr, SEXP n) {
    vbi_index_t *idx = (vbi_index_t*) R_ExternalPtrAddr(idx_ptr);
    if (!idx) Rf_error("[VBI] Index pointer is NULL");
    int nlines = asInteger(n);
    vbi_index_print(idx, nlines);
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
        SEXP out = PROTECT(allocVector(STRSXP, 0));
        UNPROTECT(1);
        return out;
    }
    // For uniformity, seek and extract records from VCF for the found indices
    htsFile *fp = hts_open(vcf, "r");
    if (!fp) {
        free(hits);
        Rf_error("Failed to open VCF/BCF: %s", vcf);
    }
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (!hdr) {
        hts_close(fp);
        free(hits);
        Rf_error("Failed to read VCF/BCF header: %s", vcf);
    }
    SEXP out = PROTECT(allocVector(STRSXP, nfound));
    bcf1_t *rec = bcf_init();
    for (int i = 0; i < nfound; ++i) {
        int idx_var = hits[i];
        int seek_ok = 0;
        if (fp->format.compression == bgzf) {
            BGZF *bg = (BGZF *)fp->fp.bgzf;
            seek_ok = (bgzf_seek(bg, idx->offsets[idx_var], SEEK_SET) == 0);
        } else {
            hFILE *hf = (hFILE *)fp->fp.hfile;
            seek_ok = (hseek(hf, (off_t)idx->offsets[idx_var], SEEK_SET) == 0);
        }
        if (!seek_ok) {
            SET_STRING_ELT(out, i, mkChar(""));
            continue;
        }
        int r = bcf_read(fp, hdr, rec);
        if (r < 0) {
            SET_STRING_ELT(out, i, mkChar(""));
            continue;
        }
        bcf_unpack(rec, BCF_UN_STR);
        kstring_t str = {0,0,0};
        vcf_format1(hdr, rec, &str);
        SET_STRING_ELT(out, i, mkChar(str.s));
        free(str.s);
    }
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    hts_close(fp);
    free(hits);
    UNPROTECT(1);
    return out;
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
        char *chrom; int start, end, label;
        cr_extract_one(ptr->cr, idx, &chrom, &start, &end, &label);
        SET_STRING_ELT(chroms, i, mkChar(chrom));
        INTEGER(starts)[i] = start;
        INTEGER(ends)[i] = end;
        INTEGER(labels)[i] = label + 1; // 1-based for R
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