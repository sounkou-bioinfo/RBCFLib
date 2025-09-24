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

// VCF/BCF header metadata structures
// Move these above VBIVcfContext

typedef enum {
    VCF_TYPE_INTEGER,
    VCF_TYPE_FLOAT,
    VCF_TYPE_FLAG,
    VCF_TYPE_CHARACTER,
    VCF_TYPE_STRING
} VcfFieldType;

typedef struct {
    char *id;
    char *number;
    VcfFieldType type;
    char *description;
    char *source;
    char *version;
} VcfInfoField;

typedef struct {
    char *id;
    char *number;
    VcfFieldType type;
    char *description;
} VcfFormatField;

typedef struct {
    char *id;
    char *description;
} VcfFilterField;

typedef struct {
    char *id;
    char *description;
} VcfAltField;

typedef struct {
    char *id;
    char *length;
    char *md5;
    char *url;
} VcfContigField;

typedef struct {
    VcfInfoField *info_fields;
    int n_info;
    VcfFormatField *format_fields;
    int n_format;
    VcfFilterField *filter_fields;
    int n_filter;
    VcfAltField *alt_fields;
    int n_alt;
    VcfContigField *contig_fields;
    int n_contig;
    char **sample_names;
    int n_samples;
    char *fileformat;
    char *header_line;
} VcfHeaderMetadata;

// VBI VCF Context structure definition
// Only define once
typedef struct vbi_vcf_context_t {
    htsFile *fp;
    bcf_hdr_t *hdr;
    vbi_index_t *vbi_idx;
    bcf1_t *tmp_ctx;  // temp variant context for reading
    kstring_t tmp_line;
    int query_failed;
    VcfHeaderMetadata *header_meta;
} VBIVcfContext, *VBIVcfContextPtr;

// CGRanges pointer type definition
// Move this above its first usage
typedef struct {
    cgranges_t *cr;
} cgranges_ptr_t;

// Forward declaration for cr_extract_one
static void cr_extract_one(cgranges_t *cr, int idx, char **chrom, int *start, int *end, int *label);

// Forward declarations
SEXP RC_cgranges_destroy(SEXP cr_ptr);
// Replace old forward declaration with new helper signature
static SEXP vbi_query_variants_basic(VBIVcfContextPtr ctx, int *hits, int nfound,
                                    int inc_info, int inc_format, int inc_genotypes);
// Forward declarations from vbi_index.c and vbi_index_capi.h
int do_index(const char *infile, const char *outfile, int n_threads);
SEXP RC_VBI_index_memory_usage(SEXP extPtr);
void RC_VBI_vcf_context_finalizer(SEXP extPtr);

// VBI indexing wrapper function
SEXP RC_VBI_index(SEXP vcf_path, SEXP vbi_path, SEXP threads) {
    const char *vcf = CHAR(STRING_ELT(vcf_path, 0));
    const char *vbi = CHAR(STRING_ELT(vbi_path, 0));
    int nthreads = asInteger(threads);
    
    int result = do_index(vcf, vbi, nthreads);
    if (result != 0) {
        Rf_error("[VBI] Indexing failed for %s", vcf);
    }
    
    return R_NilValue;
}

// Load VBI index from file
SEXP RC_VBI_load_index(SEXP vbi_path) {
    const char *path = CHAR(STRING_ELT(vbi_path, 0));
    vbi_index_t *idx = vbi_index_load(path);
    if (!idx) {
        Rf_error("[VBI] Failed to load index from %s", path);
    }
    
    SEXP extPtr = PROTECT(R_MakeExternalPtr(idx, R_NilValue, R_NilValue));
    R_RegisterCFinalizerEx(extPtr, (R_CFinalizer_t)vbi_index_finalizer, 1);
    UNPROTECT(1);
    return extPtr;
}


// Helper: Map VCF type string to enum
static VcfFieldType vcf_type_from_str(const char *type) {
    if (!type) return VCF_TYPE_STRING;
    if (strcmp(type, "Integer") == 0) return VCF_TYPE_INTEGER;
    if (strcmp(type, "Float") == 0) return VCF_TYPE_FLOAT;
    if (strcmp(type, "Flag") == 0) return VCF_TYPE_FLAG;
    if (strcmp(type, "Character") == 0) return VCF_TYPE_CHARACTER;
    if (strcmp(type, "String") == 0) return VCF_TYPE_STRING;
    return VCF_TYPE_STRING;
}

// Helper: Parse meta-information lines from bcf_hdr_t
static VcfHeaderMetadata *parse_vcf_header_metadata(bcf_hdr_t *hdr) {
    VcfHeaderMetadata *meta = calloc(1, sizeof(VcfHeaderMetadata));
    bcf_hrec_t *fileformat_hrec = bcf_hdr_get_hrec(hdr, BCF_HL_GEN, "fileformat", NULL, NULL);
    if (fileformat_hrec && fileformat_hrec->value)
        meta->fileformat = strdup(fileformat_hrec->value);
    meta->n_info = 0;
    for (int i = 0; i < hdr->nhrec; i++) {
        bcf_hrec_t *hrec = hdr->hrec[i];
        if (strcmp(hrec->key, "INFO") == 0) meta->n_info++;
    }
    meta->info_fields = calloc(meta->n_info, sizeof(VcfInfoField));
    int info_idx = 0;
    for (int i = 0; i < hdr->nhrec; i++) {
        bcf_hrec_t *hrec = hdr->hrec[i];
        if (strcmp(hrec->key, "INFO") == 0) {
            VcfInfoField *f = &meta->info_fields[info_idx++];
            int idx_id = bcf_hrec_find_key(hrec, "ID");
            int idx_num = bcf_hrec_find_key(hrec, "Number");
            int idx_type = bcf_hrec_find_key(hrec, "Type");
            int idx_desc = bcf_hrec_find_key(hrec, "Description");
            int idx_src = bcf_hrec_find_key(hrec, "Source");
            int idx_ver = bcf_hrec_find_key(hrec, "Version");
            f->id = (idx_id >= 0 && hrec->vals[idx_id]) ? strdup(hrec->vals[idx_id]) : NULL;
            f->number = (idx_num >= 0 && hrec->vals[idx_num]) ? strdup(hrec->vals[idx_num]) : NULL;
            f->type = (idx_type >= 0 && hrec->vals[idx_type]) ? vcf_type_from_str(hrec->vals[idx_type]) : VCF_TYPE_STRING;
            f->description = (idx_desc >= 0 && hrec->vals[idx_desc]) ? strdup(hrec->vals[idx_desc]) : NULL;
            f->source = (idx_src >= 0 && hrec->vals[idx_src]) ? strdup(hrec->vals[idx_src]) : NULL;
            f->version = (idx_ver >= 0 && hrec->vals[idx_ver]) ? strdup(hrec->vals[idx_ver]) : NULL;
        }
    }
    meta->n_format = 0;
    for (int i = 0; i < hdr->nhrec; i++) {
        bcf_hrec_t *hrec = hdr->hrec[i];
        if (strcmp(hrec->key, "FORMAT") == 0) meta->n_format++;
    }
    meta->format_fields = calloc(meta->n_format, sizeof(VcfFormatField));
    int fmt_idx = 0;
    for (int i = 0; i < hdr->nhrec; i++) {
        bcf_hrec_t *hrec = hdr->hrec[i];
        if (strcmp(hrec->key, "FORMAT") == 0) {
            VcfFormatField *f = &meta->format_fields[fmt_idx++];
            int idx_id = bcf_hrec_find_key(hrec, "ID");
            int idx_num = bcf_hrec_find_key(hrec, "Number");
            int idx_type = bcf_hrec_find_key(hrec, "Type");
            int idx_desc = bcf_hrec_find_key(hrec, "Description");
            f->id = (idx_id >= 0 && hrec->vals[idx_id]) ? strdup(hrec->vals[idx_id]) : NULL;
            f->number = (idx_num >= 0 && hrec->vals[idx_num]) ? strdup(hrec->vals[idx_num]) : NULL;
            f->type = (idx_type >= 0 && hrec->vals[idx_type]) ? vcf_type_from_str(hrec->vals[idx_type]) : VCF_TYPE_STRING;
            f->description = (idx_desc >= 0 && hrec->vals[idx_desc]) ? strdup(hrec->vals[idx_desc]) : NULL;
        }
    }
    meta->n_filter = 0;
    for (int i = 0; i < hdr->nhrec; i++) {
        bcf_hrec_t *hrec = hdr->hrec[i];
        if (strcmp(hrec->key, "FILTER") == 0) meta->n_filter++;
    }
    meta->filter_fields = calloc(meta->n_filter, sizeof(VcfFilterField));
    int filter_idx = 0;
    for (int i = 0; i < hdr->nhrec; i++) {
        bcf_hrec_t *hrec = hdr->hrec[i];
        if (strcmp(hrec->key, "FILTER") == 0) {
            VcfFilterField *f = &meta->filter_fields[filter_idx++];
            int idx_id = bcf_hrec_find_key(hrec, "ID");
            int idx_desc = bcf_hrec_find_key(hrec, "Description");
            f->id = (idx_id >= 0 && hrec->vals[idx_id]) ? strdup(hrec->vals[idx_id]) : NULL;
            f->description = (idx_desc >= 0 && hrec->vals[idx_desc]) ? strdup(hrec->vals[idx_desc]) : NULL;
        }
    }
    meta->n_alt = 0;
    for (int i = 0; i < hdr->nhrec; i++) {
        bcf_hrec_t *hrec = hdr->hrec[i];
        if (strcmp(hrec->key, "ALT") == 0) meta->n_alt++;
    }
    meta->alt_fields = calloc(meta->n_alt, sizeof(VcfAltField));
    int alt_idx = 0;
    for (int i = 0; i < hdr->nhrec; i++) {
        bcf_hrec_t *hrec = hdr->hrec[i];
        if (strcmp(hrec->key, "ALT") == 0) {
            VcfAltField *f = &meta->alt_fields[alt_idx++];
            int idx_id = bcf_hrec_find_key(hrec, "ID");
            int idx_desc = bcf_hrec_find_key(hrec, "Description");
            f->id = (idx_id >= 0 && hrec->vals[idx_id]) ? strdup(hrec->vals[idx_id]) : NULL;
            f->description = (idx_desc >= 0 && hrec->vals[idx_desc]) ? strdup(hrec->vals[idx_desc]) : NULL;
        }
    }
    meta->n_contig = 0;
    for (int i = 0; i < hdr->nhrec; i++) {
        bcf_hrec_t *hrec = hdr->hrec[i];
        if (strcmp(hrec->key, "contig") == 0) meta->n_contig++;
    }
    meta->contig_fields = calloc(meta->n_contig, sizeof(VcfContigField));
    int contig_idx = 0;
    for (int i = 0; i < hdr->nhrec; i++) {
        bcf_hrec_t *hrec = hdr->hrec[i];
        if (strcmp(hrec->key, "contig") == 0) {
            VcfContigField *f = &meta->contig_fields[contig_idx++];
            int idx_id = bcf_hrec_find_key(hrec, "ID");
            int idx_length = bcf_hrec_find_key(hrec, "length");
            int idx_md5 = bcf_hrec_find_key(hrec, "md5");
            int idx_url = bcf_hrec_find_key(hrec, "URL");
            f->id = (idx_id >= 0 && hrec->vals[idx_id]) ? strdup(hrec->vals[idx_id]) : NULL;
            f->length = (idx_length >= 0 && hrec->vals[idx_length]) ? strdup(hrec->vals[idx_length]) : NULL;
            f->md5 = (idx_md5 >= 0 && hrec->vals[idx_md5]) ? strdup(hrec->vals[idx_md5]) : NULL;
            f->url = (idx_url >= 0 && hrec->vals[idx_url]) ? strdup(hrec->vals[idx_url]) : NULL;
        }
    }
    meta->n_samples = bcf_hdr_nsamples(hdr);
    meta->sample_names = calloc(meta->n_samples, sizeof(char*));
    for (int i = 0; i < meta->n_samples; i++) {
        meta->sample_names[i] = strdup(bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, i));
    }
    meta->header_line = strdup(""); // hdr->samples is not a string, so just set to empty string or remove if not needed
    return meta;
}

static void free_vcf_header_metadata(VcfHeaderMetadata *meta) {
    if (!meta) return;
    for (int i = 0; i < meta->n_info; i++) {
        free(meta->info_fields[i].id);
        free(meta->info_fields[i].number);
        free(meta->info_fields[i].description);
        free(meta->info_fields[i].source);
        free(meta->info_fields[i].version);
    }
    free(meta->info_fields);
    for (int i = 0; i < meta->n_format; i++) {
        free(meta->format_fields[i].id);
        free(meta->format_fields[i].number);
        free(meta->format_fields[i].description);
    }
    free(meta->format_fields);
    for (int i = 0; i < meta->n_filter; i++) {
        free(meta->filter_fields[i].id);
        free(meta->filter_fields[i].description);
    }
    free(meta->filter_fields);
    for (int i = 0; i < meta->n_alt; i++) {
        free(meta->alt_fields[i].id);
        free(meta->alt_fields[i].description);
    }
    free(meta->alt_fields);
    for (int i = 0; i < meta->n_contig; i++) {
        free(meta->contig_fields[i].id);
        free(meta->contig_fields[i].length);
        free(meta->contig_fields[i].md5);
        free(meta->contig_fields[i].url);
    }
    free(meta->contig_fields);
    for (int i = 0; i < meta->n_samples; i++) {
        free(meta->sample_names[i]);
    }
    free(meta->sample_names);
    free(meta->fileformat);
    free(meta->header_line);
    free(meta);
}

// Forward declarations
SEXP RC_cgranges_destroy(SEXP cr_ptr);
// Replace old forward declaration with new helper signature
static SEXP vbi_query_variants_basic(VBIVcfContextPtr ctx, int *hits, int nfound,
                                    int inc_info, int inc_format, int inc_genotypes);
// Forward declarations from vbi_index.c and vbi_index_capi.h
int do_index(const char *infile, const char *outfile, int n_threads);
void RC_VBI_vcf_context_finalizer(SEXP extPtr);


// Print VBI index contents
SEXP RC_VBI_print_index(SEXP vbi_ptr, SEXP n) {
    vbi_index_t *idx = (vbi_index_t*) R_ExternalPtrAddr(vbi_ptr);
    if (!idx) {
        Rf_error("[VBI] Index pointer is NULL");
    }
    
    int num_lines = asInteger(n);
    vbi_index_print(idx, num_lines);
    
    return R_NilValue;
}

// Create VBI VCF context combining VCF file and VBI index
SEXP RC_VBI_vcf_load(SEXP vcf_path, SEXP vbi_path) {
    const char *vcf = CHAR(STRING_ELT(vcf_path, 0));
    const char *vbi = NULL;
    char *auto_vbi = NULL;  // Track if we allocated memory
    
    // Handle NULL vbi_path (auto-detect VBI file)
    if (vbi_path != R_NilValue && !isNull(vbi_path)) {
        vbi = CHAR(STRING_ELT(vbi_path, 0));
    } else {
        // Auto-detect VBI file by appending .vbi to VCF path
        auto_vbi = malloc(strlen(vcf) + 5);
        strcpy(auto_vbi, vcf);
        strcat(auto_vbi, ".vbi");
        vbi = auto_vbi;
    }
    
    VBIVcfContextPtr ctx = (VBIVcfContextPtr) R_Calloc(1, VBIVcfContext);
    ctx->query_failed = 0;
    memset(&ctx->tmp_line, 0, sizeof(kstring_t));
    
    // Load VBI index, create if it doesn't exist
    ctx->vbi_idx = vbi_index_load(vbi);
    if (!ctx->vbi_idx) {
        // Try to create the index if it doesn't exist
        Rprintf("[VBI] Index not found at %s, creating...\n", vbi);
        int result = do_index(vcf, vbi, 1); // Use 1 thread for indexing
        if (result != 0) {
            if (auto_vbi) free(auto_vbi);
            R_Free(ctx);
            Rf_error("[VBI] Failed to create index for %s", vcf);
        }
        
        // Try loading again after creation
        ctx->vbi_idx = vbi_index_load(vbi);
        if (!ctx->vbi_idx) {
            if (auto_vbi) free(auto_vbi);
            R_Free(ctx);
            Rf_error("[VBI] Failed to load newly created index from %s", vbi);
        }
    }
    
    // Open VCF file
    ctx->fp = hts_open(vcf, "r");
    if (!ctx->fp) {
        if (auto_vbi) free(auto_vbi);
        vbi_index_free(ctx->vbi_idx);
        R_Free(ctx);
        Rf_error("[VBI] Failed to open VCF file %s", vcf);
    }
    
    // Read header
    ctx->hdr = bcf_hdr_read(ctx->fp);
    if (!ctx->hdr) {
        if (auto_vbi) free(auto_vbi);
        hts_close(ctx->fp);
        vbi_index_free(ctx->vbi_idx);
        R_Free(ctx);
        Rf_error("[VBI] Failed to read VCF header from %s", vcf);
    }
    
    // Initialize temp variant context
    ctx->tmp_ctx = bcf_init();
    if (!ctx->tmp_ctx) {
        if (auto_vbi) free(auto_vbi);
        bcf_hdr_destroy(ctx->hdr);
        hts_close(ctx->fp);
        vbi_index_free(ctx->vbi_idx);
        R_Free(ctx);
        Rf_error("[VBI] Failed to initialize variant context");
    }
    
    // Parse header metadata
    ctx->header_meta = parse_vcf_header_metadata(ctx->hdr);
    
    // Clean up auto-allocated memory
    if (auto_vbi) free(auto_vbi);
    
    SEXP extPtr = PROTECT(R_MakeExternalPtr(ctx, R_NilValue, R_NilValue));
    R_RegisterCFinalizerEx(extPtr, (R_CFinalizer_t)RC_VBI_vcf_context_finalizer, 1);
    UNPROTECT(1);
    return extPtr;
}

// Finalizer for VBI VCF context
void RC_VBI_vcf_context_finalizer(SEXP extPtr) {
    VBIVcfContextPtr ctx = (VBIVcfContextPtr) R_ExternalPtrAddr(extPtr);
    if (ctx) {
        if (ctx->tmp_ctx) bcf_destroy(ctx->tmp_ctx);
        if (ctx->hdr) bcf_hdr_destroy(ctx->hdr);
        if (ctx->fp) hts_close(ctx->fp);
        if (ctx->vbi_idx) vbi_index_free(ctx->vbi_idx);
        if (ctx->tmp_line.s) free(ctx->tmp_line.s);
        if (ctx->header_meta) free_vcf_header_metadata(ctx->header_meta);
        R_Free(ctx);
        R_SetExternalPtrAddr(extPtr, NULL);
    }
}

// Get samples from VBI VCF context
SEXP RC_VBI_samples(SEXP vbi_vcf_ctx) {
    VBIVcfContextPtr ctx = (VBIVcfContextPtr) R_ExternalPtrAddr(vbi_vcf_ctx);
    if (!ctx || !ctx->hdr) Rf_error("[VBI] Invalid VCF context");
    
    int nsamples = bcf_hdr_nsamples(ctx->hdr);
    SEXP samples = PROTECT(allocVector(STRSXP, nsamples));
    
    for (int i = 0; i < nsamples; i++) {
        SET_STRING_ELT(samples, i, mkChar(bcf_hdr_int2id(ctx->hdr, BCF_DT_SAMPLE, i)));
    }
    
    UNPROTECT(1);
    return samples;
}

// Get number of samples from VBI VCF context
SEXP RC_VBI_nsamples(SEXP vbi_vcf_ctx) {
    VBIVcfContextPtr ctx = (VBIVcfContextPtr) R_ExternalPtrAddr(vbi_vcf_ctx);
    if (!ctx || !ctx->hdr) Rf_error("[VBI] Invalid VCF context");
    
    return ScalarInteger(bcf_hdr_nsamples(ctx->hdr));
}

// Get sample at specific index
SEXP RC_VBI_sample_at(SEXP vbi_vcf_ctx, SEXP index) {
    VBIVcfContextPtr ctx = (VBIVcfContextPtr) R_ExternalPtrAddr(vbi_vcf_ctx);
    if (!ctx || !ctx->hdr) Rf_error("[VBI] Invalid VCF context");
    
    int idx = asInteger(index) - 1; // Convert from 1-based to 0-based
    int nsamples = bcf_hdr_nsamples(ctx->hdr);
    
    if (idx < 0 || idx >= nsamples) {
        Rf_error("[VBI] Sample index %d out of range [1, %d]", idx + 1, nsamples);
    }
    
    const char *sample = bcf_hdr_int2id(ctx->hdr, BCF_DT_SAMPLE, idx);
    return ScalarString(mkChar(sample ? sample : ""));
}

// Convert sample name to index
SEXP RC_VBI_sample2index(SEXP vbi_vcf_ctx, SEXP sample_name) {
    VBIVcfContextPtr ctx = (VBIVcfContextPtr) R_ExternalPtrAddr(vbi_vcf_ctx);
    if (!ctx || !ctx->hdr) Rf_error("[VBI] Invalid VCF context");
    
    const char *name = CHAR(STRING_ELT(sample_name, 0));
    int idx = bcf_hdr_id2int(ctx->hdr, BCF_DT_SAMPLE, name);
    
    if (idx < 0) {
        return ScalarInteger(NA_INTEGER);
    }
    
    return ScalarInteger(idx + 1); // Convert from 0-based to 1-based
}

// Get INFO fields from VBI VCF context
SEXP RC_VBI_infos(SEXP vbi_vcf_ctx) {
    VBIVcfContextPtr ctx = (VBIVcfContextPtr) R_ExternalPtrAddr(vbi_vcf_ctx);
    if (!ctx || !ctx->hdr) Rf_error("[VBI] Invalid VCF context");
    
    // Count INFO entries
    int n_info = 0;
    for (int i = 0; i < ctx->hdr->n[BCF_DT_ID]; i++) {
        if (bcf_hdr_idinfo_exists(ctx->hdr, BCF_HL_INFO, i)) n_info++;
    }
    
    SEXP infos = PROTECT(allocVector(STRSXP, n_info));
    int idx = 0;
    
    for (int i = 0; i < ctx->hdr->n[BCF_DT_ID]; i++) {
        if (bcf_hdr_idinfo_exists(ctx->hdr, BCF_HL_INFO, i)) {
            const char *key = bcf_hdr_int2id(ctx->hdr, BCF_DT_ID, i);
            SET_STRING_ELT(infos, idx++, mkChar(key ? key : ""));
        }
    }
    
    UNPROTECT(1);
    return infos;
}

// Get FORMAT fields from VBI VCF context
SEXP RC_VBI_formats(SEXP vbi_vcf_ctx) {
    VBIVcfContextPtr ctx = (VBIVcfContextPtr) R_ExternalPtrAddr(vbi_vcf_ctx);
    if (!ctx || !ctx->hdr) Rf_error("[VBI] Invalid VCF context");
    
    // Count FORMAT entries
    int n_format = 0;
    for (int i = 0; i < ctx->hdr->n[BCF_DT_ID]; i++) {
        if (bcf_hdr_idinfo_exists(ctx->hdr, BCF_HL_FMT, i)) n_format++;
    }
    
    SEXP formats = PROTECT(allocVector(STRSXP, n_format));
    int idx = 0;
    
    for (int i = 0; i < ctx->hdr->n[BCF_DT_ID]; i++) {
        if (bcf_hdr_idinfo_exists(ctx->hdr, BCF_HL_FMT, i)) {
            const char *key = bcf_hdr_int2id(ctx->hdr, BCF_DT_ID, i);
            SET_STRING_ELT(formats, idx++, mkChar(key ? key : ""));
        }
    }
    
    UNPROTECT(1);
    return formats;
}

// Get FILTER fields from VBI VCF context
SEXP RC_VBI_filters(SEXP vbi_vcf_ctx) {
    VBIVcfContextPtr ctx = (VBIVcfContextPtr) R_ExternalPtrAddr(vbi_vcf_ctx);
    if (!ctx || !ctx->hdr) Rf_error("[VBI] Invalid VCF context");
    
    // Count FILTER entries
    int n_filter = 0;
    for (int i = 0; i < ctx->hdr->n[BCF_DT_ID]; i++) {
        if (bcf_hdr_idinfo_exists(ctx->hdr, BCF_HL_FLT, i)) n_filter++;
    }
    
    SEXP filters = PROTECT(allocVector(STRSXP, n_filter));
    int idx = 0;
    
    for (int i = 0; i < ctx->hdr->n[BCF_DT_ID]; i++) {
        if (bcf_hdr_idinfo_exists(ctx->hdr, BCF_HL_FLT, i)) {
            const char *key = bcf_hdr_int2id(ctx->hdr, BCF_DT_ID, i);
            SET_STRING_ELT(filters, idx++, mkChar(key ? key : ""));
        }
    }
    
    UNPROTECT(1);
    return filters;
}

// Extract ranges from VBI index
SEXP RC_VBI_extract_ranges(SEXP idx_ptr, SEXP n) {
    vbi_index_t *idx = (vbi_index_t*) R_ExternalPtrAddr(idx_ptr);
    if (!idx) Rf_error("[VBI] Index pointer is NULL");
    
    int num_ranges = asInteger(n);
    if (num_ranges <= 0 || num_ranges > idx->num_marker) {
        num_ranges = idx->num_marker;
    }
    
    SEXP chrom_col = PROTECT(allocVector(STRSXP, num_ranges));
    SEXP pos_col = PROTECT(allocVector(INTSXP, num_ranges));
    SEXP idx_col = PROTECT(allocVector(INTSXP, num_ranges));
    
    for (int i = 0; i < num_ranges; i++) {
        const char *chrom = idx->chrom_names[idx->chrom_ids[i]];
        SET_STRING_ELT(chrom_col, i, mkChar(chrom ? chrom : ""));
        INTEGER(pos_col)[i] = (int)idx->positions[i];
        INTEGER(idx_col)[i] = i + 1; // 1-based index for R
    }
    
    SEXP df = PROTECT(allocVector(VECSXP, 3));
    SET_VECTOR_ELT(df, 0, chrom_col);
    SET_VECTOR_ELT(df, 1, pos_col);
    SET_VECTOR_ELT(df, 2, idx_col);
    
    SEXP names = PROTECT(allocVector(STRSXP, 3));
    SET_STRING_ELT(names, 0, mkChar("chrom"));
    SET_STRING_ELT(names, 1, mkChar("pos"));
    SET_STRING_ELT(names, 2, mkChar("index"));
    
    setAttrib(df, R_NamesSymbol, names);
    setAttrib(df, R_ClassSymbol, mkString("data.frame"));
    
    SEXP rn = PROTECT(allocVector(INTSXP, num_ranges));
    for (int i = 0; i < num_ranges; i++) INTEGER(rn)[i] = i + 1;
    setAttrib(df, R_RowNamesSymbol, rn);
    
    UNPROTECT(6);
    return df;
}

// Shared helper to build a data.frame for a set of variant indices (hits)
// Columns: chrom,pos,id,ref,alt,qual,filter,n_allele,index,[CSQ],[ANN]
static SEXP vbi_query_variants_basic(VBIVcfContextPtr ctx, int *hits, int nfound,
                                    int inc_info, int inc_format, int inc_genotypes) {
    if (nfound <= 0) {
        SEXP df = PROTECT(allocVector(VECSXP, 0));
        SEXP names = PROTECT(allocVector(STRSXP, 0));
        SEXP rn = PROTECT(allocVector(INTSXP, 0));
        setAttrib(df, R_NamesSymbol, names);
        setAttrib(df, R_RowNamesSymbol, rn);
        setAttrib(df, R_ClassSymbol, mkString("data.frame"));
        UNPROTECT(3);
        return df;
    }
    // Detect CSQ / ANN only if info requested
    bcf_hrec_t *csq_hrec = inc_info ? bcf_hdr_get_hrec(ctx->hdr, BCF_HL_INFO, "ID", "CSQ", NULL) : NULL;
    bcf_hrec_t *ann_hrec = inc_info ? bcf_hdr_get_hrec(ctx->hdr, BCF_HL_INFO, "ID", "ANN", NULL) : NULL;
    int with_csq = csq_hrec != NULL;
    int with_ann = ann_hrec != NULL;

    // Base columns always
    int base_cols = 9; // chrom,pos,id,ref,alt,qual,filter,n_allele,index
    int extra_cols = (with_csq?1:0) + (with_ann?1:0);
    // Simple INFO/FORMAT/GT strategy: we expose aggregated INFO as a single JSON-ish string column when requested (placeholder)
    int info_cols = inc_info ? 1 : 0; // INFO
    int format_cols = inc_format ? 1 : 0; // FORMAT ids present (per record)
    int gt_cols = inc_genotypes ? 1 : 0; // GT (concatenated per sample)
    int ncols = base_cols + extra_cols + info_cols + format_cols + gt_cols;

    SEXP df = PROTECT(allocVector(VECSXP, ncols));
    SEXP names = PROTECT(allocVector(STRSXP, ncols));
    int col = 0;
    SEXP chrom_col = PROTECT(allocVector(STRSXP, nfound));
    SEXP pos_col   = PROTECT(allocVector(INTSXP, nfound));
    SEXP id_col    = PROTECT(allocVector(STRSXP, nfound));
    SEXP ref_col   = PROTECT(allocVector(STRSXP, nfound));
    SEXP alt_col   = PROTECT(allocVector(STRSXP, nfound));
    SEXP qual_col  = PROTECT(allocVector(REALSXP, nfound));
    SEXP filter_col= PROTECT(allocVector(STRSXP, nfound));
    SEXP nallele_col=PROTECT(allocVector(INTSXP, nfound));
    SEXP index_col = PROTECT(allocVector(INTSXP, nfound));
    int protect_count = 2 + 9; // df,names + 9 columns

    SEXP csq_col = R_NilValue, ann_col = R_NilValue, info_col = R_NilValue, fmt_col = R_NilValue, gt_col = R_NilValue;
    if (with_csq) { csq_col = PROTECT(allocVector(VECSXP, nfound)); protect_count++; }
    if (with_ann) { ann_col = PROTECT(allocVector(VECSXP, nfound)); protect_count++; }
    if (info_cols) { info_col = PROTECT(allocVector(STRSXP, nfound)); protect_count++; }
    if (format_cols){ fmt_col = PROTECT(allocVector(STRSXP, nfound)); protect_count++; }
    if (gt_cols) { gt_col = PROTECT(allocVector(STRSXP, nfound)); protect_count++; }

    bcf1_t *rec = bcf_init();
    kstring_t kalt = {0,0,0};
    kstring_t kflt = {0,0,0};
    kstring_t ktmp = {0,0,0};
    for (int i=0;i<nfound;i++) {
        int idx_var = hits[i];
        int seek_ok = 0;
        if (ctx->fp->format.compression == bgzf) {
            BGZF *bg = (BGZF *)ctx->fp->fp.bgzf;
            seek_ok = (bgzf_seek(bg, ctx->vbi_idx->offsets[idx_var], SEEK_SET) == 0);
        } else {
            hFILE *hf = (hFILE *)ctx->fp->fp.hfile;
            seek_ok = (hseek(hf, (off_t)ctx->vbi_idx->offsets[idx_var], SEEK_SET) == 0);
        }
        if (!seek_ok || bcf_read(ctx->fp, ctx->hdr, rec) < 0) {
            SET_STRING_ELT(chrom_col,i,NA_STRING); INTEGER(pos_col)[i]=NA_INTEGER; SET_STRING_ELT(id_col,i,NA_STRING);
            SET_STRING_ELT(ref_col,i,NA_STRING); SET_STRING_ELT(alt_col,i,NA_STRING); REAL(qual_col)[i]=NA_REAL;
            SET_STRING_ELT(filter_col,i,NA_STRING); INTEGER(nallele_col)[i]=NA_INTEGER; INTEGER(index_col)[i]=NA_INTEGER;
            if (with_csq) SET_VECTOR_ELT(csq_col,i,R_NilValue);
            if (with_ann) SET_VECTOR_ELT(ann_col,i,R_NilValue);
            if (info_cols) SET_STRING_ELT(info_col,i,NA_STRING);
            if (format_cols) SET_STRING_ELT(fmt_col,i,NA_STRING);
            if (gt_cols) SET_STRING_ELT(gt_col,i,NA_STRING);
            continue;
        }
        bcf_unpack(rec, BCF_UN_STR|BCF_UN_INFO| (inc_format||inc_genotypes?BCF_UN_FMT:0) | BCF_UN_FLT);
        const char *chrom = bcf_hdr_id2name(ctx->hdr, rec->rid);
        SET_STRING_ELT(chrom_col,i, chrom?mkChar(chrom):NA_STRING);
        INTEGER(pos_col)[i] = rec->pos + 1;
        SET_STRING_ELT(id_col,i, (rec->d.id && strcmp(rec->d.id, ".")!=0)? mkChar(rec->d.id): NA_STRING);
        if (rec->d.allele && rec->n_allele>0){
            SET_STRING_ELT(ref_col,i, mkChar(rec->d.allele[0]));
            kalt.l=0; if (rec->n_allele>1){ for(int a=1;a<rec->n_allele;a++){ if(a>1) kputc(',',&kalt); kputs(rec->d.allele[a],&kalt);} SET_STRING_ELT(alt_col,i,mkChar(kalt.s)); } else SET_STRING_ELT(alt_col,i,mkChar("."));
        } else { SET_STRING_ELT(ref_col,i,NA_STRING); SET_STRING_ELT(alt_col,i,NA_STRING);}        
        INTEGER(nallele_col)[i]=rec->n_allele; INTEGER(index_col)[i]=idx_var+1; REAL(qual_col)[i]= bcf_float_is_missing(rec->qual)? NA_REAL: rec->qual;
        kflt.l=0; if(rec->d.n_flt==0){ kputs("PASS",&kflt);} else { for(int f=0;f<rec->d.n_flt;f++){ if(f) kputc(';',&kflt); const char *flt=bcf_hdr_int2id(ctx->hdr, BCF_DT_ID, rec->d.flt[f]); if(flt) kputs(flt,&kflt);} }
        SET_STRING_ELT(filter_col,i, kflt.s?mkChar(kflt.s):mkChar(""));
        if (with_csq) SET_VECTOR_ELT(csq_col,i,R_NilValue); // placeholder
        if (with_ann) SET_VECTOR_ELT(ann_col,i,R_NilValue); // placeholder
        if (info_cols){
            // concatenate key=value for all INFO present (simple, may be large)
            ktmp.l=0; bcf_info_t *inf=NULL; for(int j=0;j<rec->n_info;j++){ inf=&rec->d.info[j]; const char *key = bcf_hdr_int2id(ctx->hdr,BCF_DT_ID,inf->key); if(!key) continue; if(j) kputc(';',&ktmp); kputs(key,&ktmp); if(inf->len>0 && inf->type!=BCF_BT_NULL){ kputc('=',&ktmp); if(inf->type==BCF_BT_INT32){ for(int k=0;k<inf->len;k++){ if(k) kputc(',',&ktmp); kputw(((int32_t*)inf->vptr)[k],&ktmp);} } else if(inf->type==BCF_BT_FLOAT){ for(int k=0;k<inf->len;k++){ if(k) kputc(',',&ktmp); ksprintf(&ktmp,"%g",((float*)inf->vptr)[k]); } } else if(inf->type==BCF_BT_CHAR){ kputsn(inf->vptr,inf->len,&ktmp);} }
            }
            SET_STRING_ELT(info_col,i, ktmp.s?mkChar(ktmp.s):NA_STRING);
        }
        if (format_cols){
            ktmp.l=0; for(int j=0;j<rec->n_fmt;j++){ if(j) kputc(';',&ktmp); const char *fid=bcf_hdr_int2id(ctx->hdr,BCF_DT_ID,rec->d.fmt[j].id); if(fid) kputs(fid,&ktmp);} SET_STRING_ELT(fmt_col,i, ktmp.s?mkChar(ktmp.s):NA_STRING);
        }
        if (gt_cols){
            ktmp.l=0; int ngt_arr=0; int32_t *gt=NULL; int ngt=bcf_get_genotypes(ctx->hdr, rec, &gt, &ngt_arr);
            int ns= bcf_hdr_nsamples(ctx->hdr);
            if(ngt>0 && ns>0){
                int ploidy = ngt / ns;
                for(int s=0;s<ns;s++){
                    if(s) kputc(';',&ktmp);
                    for(int p=0;p<ploidy;p++){
                        if(p) kputc('/', &ktmp);
                        int32_t g = gt[s*ploidy + p];
                        if (g==bcf_int32_vector_end) break;
                        if (g==bcf_gt_missing) { kputc('.',&ktmp); continue; }
                        ksprintf(&ktmp, "%d", bcf_gt_allele(g));
                    }
                }
            }
            if(gt) free(gt);
            SET_STRING_ELT(gt_col,i, ktmp.s?mkChar(ktmp.s):NA_STRING);
        }
    }
    bcf_destroy(rec);
    if(kalt.s) free(kalt.s); if(kflt.s) free(kflt.s); if(ktmp.s) free(ktmp.s);
    // assemble df columns
    col=0; // reset column counter (fix redefinition)
    SET_VECTOR_ELT(df,col,chrom_col); SET_STRING_ELT(names,col++,mkChar("chrom"));
    SET_VECTOR_ELT(df,col,pos_col); SET_STRING_ELT(names,col++,mkChar("pos"));
    SET_VECTOR_ELT(df,col,id_col); SET_STRING_ELT(names,col++,mkChar("id"));
    SET_VECTOR_ELT(df,col,ref_col); SET_STRING_ELT(names,col++,mkChar("ref"));
    SET_VECTOR_ELT(df,col,alt_col); SET_STRING_ELT(names,col++,mkChar("alt"));
    SET_VECTOR_ELT(df,col,qual_col); SET_STRING_ELT(names,col++,mkChar("qual"));
    SET_VECTOR_ELT(df,col,filter_col); SET_STRING_ELT(names,col++,mkChar("filter"));
    SET_VECTOR_ELT(df,col,nallele_col); SET_STRING_ELT(names,col++,mkChar("n_allele"));
    SET_VECTOR_ELT(df,col,index_col); SET_STRING_ELT(names,col++,mkChar("index"));
    if (with_csq){ SET_VECTOR_ELT(df,col,csq_col); SET_STRING_ELT(names,col++,mkChar("CSQ")); }
    if (with_ann){ SET_VECTOR_ELT(df,col,ann_col); SET_STRING_ELT(names,col++,mkChar("ANN")); }
    if (info_cols){ SET_VECTOR_ELT(df,col,info_col); SET_STRING_ELT(names,col++,mkChar("INFO")); }
    if (format_cols){ SET_VECTOR_ELT(df,col,fmt_col); SET_STRING_ELT(names,col++,mkChar("FORMAT_IDS")); }
    if (gt_cols){ SET_VECTOR_ELT(df,col,gt_col); SET_STRING_ELT(names,col++,mkChar("GT")); }
    setAttrib(df, R_NamesSymbol, names);
    SEXP rn = PROTECT(allocVector(INTSXP, nfound)); protect_count++;
    for(int i=0;i<nfound;i++) INTEGER(rn)[i]=i+1;
    setAttrib(df,R_RowNamesSymbol,rn);
    setAttrib(df,R_ClassSymbol,mkString("data.frame"));
    UNPROTECT(protect_count);
    return df;
}

// Print header metadata from VBI VCF context
SEXP RC_VBI_print_header_metadata(SEXP vbi_vcf_ctx_ext) {
    VBIVcfContextPtr ctx = (VBIVcfContextPtr) R_ExternalPtrAddr(vbi_vcf_ctx_ext);
    if (!ctx || !ctx->header_meta) Rf_error("[VBI] Invalid VCF context or missing header metadata");
    VcfHeaderMetadata *meta = ctx->header_meta;
    SEXP out = PROTECT(allocVector(VECSXP, 7));
    SEXP names = PROTECT(allocVector(STRSXP, 7));
    SET_STRING_ELT(names, 0, mkChar("fileformat"));
    SET_STRING_ELT(names, 1, mkChar("info_fields"));
    SET_STRING_ELT(names, 2, mkChar("format_fields"));
    SET_STRING_ELT(names, 3, mkChar("filter_fields"));
    SET_STRING_ELT(names, 4, mkChar("alt_fields"));
    SET_STRING_ELT(names, 5, mkChar("contig_fields"));
    SET_STRING_ELT(names, 6, mkChar("sample_names"));
    SET_VECTOR_ELT(out, 0, ScalarString(mkChar(meta->fileformat ? meta->fileformat : "")));
    // info_fields
    SEXP info_df = PROTECT(allocVector(VECSXP, 6));
    SEXP info_names = PROTECT(allocVector(STRSXP, 6));
    SET_STRING_ELT(info_names, 0, mkChar("id"));
    SET_STRING_ELT(info_names, 1, mkChar("number"));
    SET_STRING_ELT(info_names, 2, mkChar("type"));
    SET_STRING_ELT(info_names, 3, mkChar("description"));
    SET_STRING_ELT(info_names, 4, mkChar("source"));
    SET_STRING_ELT(info_names, 5, mkChar("version"));
    SEXP id_col = PROTECT(allocVector(STRSXP, meta->n_info));
    SEXP num_col = PROTECT(allocVector(STRSXP, meta->n_info));
    SEXP type_col = PROTECT(allocVector(STRSXP, meta->n_info));
    SEXP desc_col = PROTECT(allocVector(STRSXP, meta->n_info));
    SEXP src_col = PROTECT(allocVector(STRSXP, meta->n_info));
    SEXP ver_col = PROTECT(allocVector(STRSXP, meta->n_info));
    for (int i = 0; i < meta->n_info; i++) {
        SET_STRING_ELT(id_col, i, mkChar(meta->info_fields[i].id ? meta->info_fields[i].id : ""));
        SET_STRING_ELT(num_col, i, mkChar(meta->info_fields[i].number ? meta->info_fields[i].number : ""));
        SET_STRING_ELT(type_col, i, mkChar(
            meta->info_fields[i].type == VCF_TYPE_INTEGER ? "Integer" :
            meta->info_fields[i].type == VCF_TYPE_FLOAT ? "Float" :
            meta->info_fields[i].type == VCF_TYPE_FLAG ? "Flag" :
            meta->info_fields[i].type == VCF_TYPE_CHARACTER ? "Character" : "String"
        ));
        SET_STRING_ELT(desc_col, i, mkChar(meta->info_fields[i].description ? meta->info_fields[i].description : ""));
        SET_STRING_ELT(src_col, i, mkChar(meta->info_fields[i].source ? meta->info_fields[i].source : ""));
        SET_STRING_ELT(ver_col, i, mkChar(meta->info_fields[i].version ? meta->info_fields[i].version : ""));
    }
    SET_VECTOR_ELT(info_df, 0, id_col);
    SET_VECTOR_ELT(info_df, 1, num_col);
    SET_VECTOR_ELT(info_df, 2, type_col);
    SET_VECTOR_ELT(info_df, 3, desc_col);
    SET_VECTOR_ELT(info_df, 4, src_col);
    SET_VECTOR_ELT(info_df, 5, ver_col);
    setAttrib(info_df, R_NamesSymbol, info_names);
    SET_VECTOR_ELT(out, 1, info_df);
    // format_fields
    SEXP fmt_df = PROTECT(allocVector(VECSXP, 4));
    SEXP fmt_names = PROTECT(allocVector(STRSXP, 4));
    SET_STRING_ELT(fmt_names, 0, mkChar("id"));
    SET_STRING_ELT(fmt_names, 1, mkChar("number"));
    SET_STRING_ELT(fmt_names, 2, mkChar("type"));
    SET_STRING_ELT(fmt_names, 3, mkChar("description"));
    SEXP fid_col = PROTECT(allocVector(STRSXP, meta->n_format));
    SEXP fnum_col = PROTECT(allocVector(STRSXP, meta->n_format));
    SEXP ftype_col = PROTECT(allocVector(STRSXP, meta->n_format));
    SEXP fdesc_col = PROTECT(allocVector(STRSXP, meta->n_format));
    for (int i = 0; i < meta->n_format; i++) {
        SET_STRING_ELT(fid_col, i, mkChar(meta->format_fields[i].id ? meta->format_fields[i].id : ""));
        SET_STRING_ELT(fnum_col, i, mkChar(meta->format_fields[i].number ? meta->format_fields[i].number : ""));
        SET_STRING_ELT(ftype_col, i, mkChar(
            meta->format_fields[i].type == VCF_TYPE_INTEGER ? "Integer" :
            meta->format_fields[i].type == VCF_TYPE_FLOAT ? "Float" :
            meta->format_fields[i].type == VCF_TYPE_FLAG ? "Flag" :
            meta->format_fields[i].type == VCF_TYPE_CHARACTER ? "Character" : "String"
        ));
        SET_STRING_ELT(fdesc_col, i, mkChar(meta->format_fields[i].description ? meta->format_fields[i].description : ""));
    }
    SET_VECTOR_ELT(fmt_df, 0, fid_col);
    SET_VECTOR_ELT(fmt_df, 1, fnum_col);
    SET_VECTOR_ELT(fmt_df, 2, ftype_col);
    SET_VECTOR_ELT(fmt_df, 3, fdesc_col);
    setAttrib(fmt_df, R_NamesSymbol, fmt_names);
    SET_VECTOR_ELT(out, 2, fmt_df);
    // filter_fields
    SEXP filter_df = PROTECT(allocVector(VECSXP, 2));
    SEXP filter_names = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(filter_names, 0, mkChar("id"));
    SET_STRING_ELT(filter_names, 1, mkChar("description"));
    SEXP fid_col2 = PROTECT(allocVector(STRSXP, meta->n_filter));
    SEXP fdesc_col2 = PROTECT(allocVector(STRSXP, meta->n_filter));
    for (int i = 0; i < meta->n_filter; i++) {
        SET_STRING_ELT(fid_col2, i, mkChar(meta->filter_fields[i].id ? meta->filter_fields[i].id : ""));
        SET_STRING_ELT(fdesc_col2, i, mkChar(meta->filter_fields[i].description ? meta->filter_fields[i].description : ""));
    }
    SET_VECTOR_ELT(filter_df, 0, fid_col2);
    SET_VECTOR_ELT(filter_df, 1, fdesc_col2);
    setAttrib(filter_df, R_NamesSymbol, filter_names);
    SET_VECTOR_ELT(out, 3, filter_df);
    // alt_fields
    SEXP alt_df = PROTECT(allocVector(VECSXP, 2));
    SEXP alt_names = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(alt_names, 0, mkChar("id"));
    SET_STRING_ELT(alt_names, 1, mkChar("description"));
    SEXP aid_col = PROTECT(allocVector(STRSXP, meta->n_alt));
    SEXP adesc_col = PROTECT(allocVector(STRSXP, meta->n_alt));
    for (int i = 0; i < meta->n_alt; i++) {
        SET_STRING_ELT(aid_col, i, mkChar(meta->alt_fields[i].id ? meta->alt_fields[i].id : ""));
        SET_STRING_ELT(adesc_col, i, mkChar(meta->alt_fields[i].description ? meta->alt_fields[i].description : ""));
    }
    SET_VECTOR_ELT(alt_df, 0, aid_col);
    SET_VECTOR_ELT(alt_df, 1, adesc_col);
    setAttrib(alt_df, R_NamesSymbol, alt_names);
    SET_VECTOR_ELT(out, 4, alt_df);
    // contig_fields
    SEXP contig_df = PROTECT(allocVector(VECSXP, 4));
    SEXP contig_names = PROTECT(allocVector(STRSXP, 4));
    SET_STRING_ELT(contig_names, 0, mkChar("id"));
    SET_STRING_ELT(contig_names, 1, mkChar("length"));
    SET_STRING_ELT(contig_names, 2, mkChar("md5"));
    SET_STRING_ELT(contig_names, 3, mkChar("url"));
    SEXP cid_col = PROTECT(allocVector(STRSXP, meta->n_contig));
    SEXP clen_col = PROTECT(allocVector(STRSXP, meta->n_contig));
    SEXP cmd5_col = PROTECT(allocVector(STRSXP, meta->n_contig));
    SEXP curl_col = PROTECT(allocVector(STRSXP, meta->n_contig));
    for (int i = 0; i < meta->n_contig; i++) {
        SET_STRING_ELT(cid_col, i, mkChar(meta->contig_fields[i].id ? meta->contig_fields[i].id : ""));
        SET_STRING_ELT(clen_col, i, mkChar(meta->contig_fields[i].length ? meta->contig_fields[i].length : ""));
        SET_STRING_ELT(cmd5_col, i, mkChar(meta->contig_fields[i].md5 ? meta->contig_fields[i].md5 : ""));
        SET_STRING_ELT(curl_col, i, mkChar(meta->contig_fields[i].url ? meta->contig_fields[i].url : ""));
    }
    SET_VECTOR_ELT(contig_df, 0, cid_col);
    SET_VECTOR_ELT(contig_df, 1, clen_col);
    SET_VECTOR_ELT(contig_df, 2, cmd5_col);
    SET_VECTOR_ELT(contig_df, 3, curl_col);
    setAttrib(contig_df, R_NamesSymbol, contig_names);
    SET_VECTOR_ELT(out, 5, contig_df);
    // sample_names
    SEXP samples = PROTECT(allocVector(STRSXP, meta->n_samples));
    for (int i = 0; i < meta->n_samples; i++)
        SET_STRING_ELT(samples, i, mkChar(meta->sample_names[i] ? meta->sample_names[i] : ""));
    SET_VECTOR_ELT(out, 6, samples);
    setAttrib(out, R_NamesSymbol, names);
    UNPROTECT(12);
    return out;
}

// Backwards-compatible by-indices query reopening VCF each time (no info/format/genotypes)
SEXP RC_VBI_query_by_indices(SEXP vcf_path, SEXP idx_ptr, SEXP start_idx, SEXP end_idx, SEXP threads) {
    const char *vcf = CHAR(STRING_ELT(vcf_path, 0));
    vbi_index_t *idx = (vbi_index_t*) R_ExternalPtrAddr(idx_ptr);
    if (!idx) Rf_error("[VBI] Index pointer is NULL");
    int start = asInteger(start_idx) - 1;
    int end = asInteger(end_idx) - 1;
    if (start < 0) start = 0;
    if (end >= idx->num_marker) end = idx->num_marker - 1;
    if (end < start) {
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
    int *hits = (int*) malloc(sizeof(int)*nfound); if(!hits) Rf_error("[VBI] OOM");
    for(int i=0;i<nfound;i++) hits[i]= start + i;
    VBIVcfContextPtr ctx = (VBIVcfContextPtr) R_Calloc(1, VBIVcfContext);
    ctx->query_failed=0; memset(&ctx->tmp_line,0,sizeof(kstring_t));
    ctx->fp = hts_open(vcf, "r"); if(!ctx->fp){ R_Free(ctx); free(hits); Rf_error("Failed to open %s", vcf);}    
    ctx->hdr = bcf_hdr_read(ctx->fp); if(!ctx->hdr){ hts_close(ctx->fp); R_Free(ctx); free(hits); Rf_error("Failed header %s", vcf);}    
    ctx->vbi_idx = idx; ctx->tmp_ctx = bcf_init(); if(!ctx->tmp_ctx){ bcf_hdr_destroy(ctx->hdr); hts_close(ctx->fp); R_Free(ctx); free(hits); Rf_error("[VBI] OOM variant ctx"); }
    SEXP out = vbi_query_variants_basic(ctx, hits, nfound, 0,0,0);
    if(ctx->tmp_ctx) bcf_destroy(ctx->tmp_ctx); if(ctx->hdr) bcf_hdr_destroy(ctx->hdr); if(ctx->fp) hts_close(ctx->fp); if(ctx->tmp_line.s) free(ctx->tmp_line.s); R_Free(ctx); free(hits);
    return out;
}

// Linear scan region query using VBI index
SEXP RC_VBI_query_region(SEXP vbi_vcf_ctx, SEXP region_str, SEXP include_info, SEXP include_format, SEXP include_genotypes) {
    VBIVcfContextPtr ctx = (VBIVcfContextPtr) R_ExternalPtrAddr(vbi_vcf_ctx);
    if (!ctx) Rf_error("[VBI] VCF context pointer is NULL");
    if (!ctx->vbi_idx) Rf_error("[VBI] No VBI index available in context");
    int inc_info = asLogical(include_info);
    int inc_format = asLogical(include_format);
    int inc_genotypes = asLogical(include_genotypes);
    const char *reg = CHAR(STRING_ELT(region_str, 0));
    int nfound = 0;
    int *indices = vbi_index_query_region(ctx->vbi_idx, reg, &nfound);
    if (!indices || nfound == 0) {
        if (indices) free(indices);
        SEXP df = PROTECT(allocVector(VECSXP, 0));
        SEXP names = PROTECT(allocVector(STRSXP, 0));
        SEXP rn = PROTECT(allocVector(INTSXP, 0));
        setAttrib(df, R_NamesSymbol, names);
        setAttrib(df, R_RowNamesSymbol, rn);
        setAttrib(df, R_ClassSymbol, mkString("data.frame"));
        UNPROTECT(3);
        return df;
    }
    SEXP out = vbi_query_variants_basic(ctx, indices, nfound, inc_info, inc_format, inc_genotypes);
    free(indices);
    return out;
}

// Linear scan into VBI index for matches in a range
SEXP RC_VBI_query_range(SEXP vbi_vcf_ctx, SEXP chrom, SEXP start, SEXP end, SEXP include_info, SEXP include_format, SEXP include_genotypes) {
    VBIVcfContextPtr ctx = (VBIVcfContextPtr) R_ExternalPtrAddr(vbi_vcf_ctx);
    if (!ctx) Rf_error("[VBI] VCF context pointer is NULL");
    if (!ctx->vbi_idx) Rf_error("[VBI] No VBI index available in context");
    const char *chr = CHAR(STRING_ELT(chrom, 0));
    int st = asInteger(start);
    int en = asInteger(end);
    int inc_info = asLogical(include_info);
    int inc_format = asLogical(include_format);
    int inc_geno = asLogical(include_genotypes);
    int nfound = 0;
    // Linear scan: collect all indices matching chrom and position in [st, en]
    vbi_index_t *idx = ctx->vbi_idx;
    int *hits = (int*) malloc(sizeof(int) * idx->num_marker);
    if (!hits) Rf_error("[VBI] OOM");
    for (int i = 0; i < idx->num_marker; i++) {
        const char *c = idx->chrom_names[idx->chrom_ids[i]];
        int64_t pos = idx->positions[i];
        if (strcmp(c, chr) == 0 && pos >= st && pos <= en) {
            hits[nfound++] = i;
        }
    }
    SEXP out;
    if (nfound > 0) {
        out = vbi_query_variants_basic(ctx, hits, nfound, inc_info, inc_format, inc_geno);
    } else {
        SEXP df = PROTECT(allocVector(VECSXP, 0));
        SEXP names = PROTECT(allocVector(STRSXP, 0));
        SEXP rn = PROTECT(allocVector(INTSXP, 0));
        setAttrib(df, R_NamesSymbol, names);
        setAttrib(df, R_RowNamesSymbol, rn);
        setAttrib(df, R_ClassSymbol, mkString("data.frame"));
        UNPROTECT(3);
        out = df;
    }
    free(hits);
    return out;
}

// Minimal cgranges R binding
// 
//
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
    
    // Find which contig this interval belongs to by checking the offset ranges
    int32_t ctg_id = -1;
    for (int i = 0; i < cr->n_ctg; i++) {
        if (idx >= cr->ctg[i].off && idx < cr->ctg[i].off + cr->ctg[i].n) {
            ctg_id = i;
            break;
        }
    }
    
    // Extract contig name
    if (ctg_id >= 0 && ctg_id < cr->n_ctg && cr->ctg[ctg_id].name) {
        *chrom = cr->ctg[ctg_id].name;
    } else {
        *chrom = NULL;
    }
    
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

// Header info for VCF/BCF file (no VBI)
SEXP RC_VCF_header_info(SEXP vcf_path) {
    const char *vcf = CHAR(STRING_ELT(vcf_path, 0));
    htsFile *fp = hts_open(vcf, "r");
    if (!fp) Rf_error("[VCF] Failed to open file %s", vcf);
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (!hdr) {
        hts_close(fp);
        Rf_error("[VCF] Failed to read header from %s", vcf);
    }
    VcfHeaderMetadata *meta = parse_vcf_header_metadata(hdr);
    bcf_hdr_destroy(hdr);
    hts_close(fp);
    // Build R list
    SEXP out = PROTECT(allocVector(VECSXP, 7));
    SEXP names = PROTECT(allocVector(STRSXP, 7));
    SET_STRING_ELT(names, 0, mkChar("fileformat"));
    SET_STRING_ELT(names, 1, mkChar("info_fields"));
    SET_STRING_ELT(names, 2, mkChar("format_fields"));
    SET_STRING_ELT(names, 3, mkChar("filter_fields"));
    SET_STRING_ELT(names, 4, mkChar("alt_fields"));
    SET_STRING_ELT(names, 5, mkChar("contig_fields"));
    SET_STRING_ELT(names, 6, mkChar("sample_names"));
    SET_VECTOR_ELT(out, 0, ScalarString(mkChar(meta->fileformat ? meta->fileformat : "")));
    // info_fields
    SEXP info_df = PROTECT(allocVector(VECSXP, 6));
    SEXP info_names = PROTECT(allocVector(STRSXP, 6));
    SET_STRING_ELT(info_names, 0, mkChar("id"));
    SET_STRING_ELT(info_names, 1, mkChar("number"));
    SET_STRING_ELT(info_names, 2, mkChar("type"));
    SET_STRING_ELT(info_names, 3, mkChar("description"));
    SET_STRING_ELT(info_names, 4, mkChar("source"));
    SET_STRING_ELT(info_names, 5, mkChar("version"));
    SEXP id_col = PROTECT(allocVector(STRSXP, meta->n_info));
    SEXP num_col = PROTECT(allocVector(STRSXP, meta->n_info));
    SEXP type_col = PROTECT(allocVector(STRSXP, meta->n_info));
    SEXP desc_col = PROTECT(allocVector(STRSXP, meta->n_info));
    SEXP src_col = PROTECT(allocVector(STRSXP, meta->n_info));
    SEXP ver_col = PROTECT(allocVector(STRSXP, meta->n_info));
    for (int i = 0; i < meta->n_info; i++) {
        SET_STRING_ELT(id_col, i, mkChar(meta->info_fields[i].id ? meta->info_fields[i].id : ""));
        SET_STRING_ELT(num_col, i, mkChar(meta->info_fields[i].number ? meta->info_fields[i].number : ""));
        SET_STRING_ELT(type_col, i, mkChar(
            meta->info_fields[i].type == VCF_TYPE_INTEGER ? "Integer" :
            meta->info_fields[i].type == VCF_TYPE_FLOAT ? "Float" :
            meta->info_fields[i].type == VCF_TYPE_FLAG ? "Flag" :
            meta->info_fields[i].type == VCF_TYPE_CHARACTER ? "Character" : "String"
        ));
        SET_STRING_ELT(desc_col, i, mkChar(meta->info_fields[i].description ? meta->info_fields[i].description : ""));
        SET_STRING_ELT(src_col, i, mkChar(meta->info_fields[i].source ? meta->info_fields[i].source : ""));
        SET_STRING_ELT(ver_col, i, mkChar(meta->info_fields[i].version ? meta->info_fields[i].version : ""));
    }
    SET_VECTOR_ELT(info_df, 0, id_col);
    SET_VECTOR_ELT(info_df, 1, num_col);
    SET_VECTOR_ELT(info_df, 2, type_col);
    SET_VECTOR_ELT(info_df, 3, desc_col);
    SET_VECTOR_ELT(info_df, 4, src_col);
    SET_VECTOR_ELT(info_df, 5, ver_col);
    setAttrib(info_df, R_NamesSymbol, info_names);
    SET_VECTOR_ELT(out, 1, info_df);
    // format_fields
    SEXP fmt_df = PROTECT(allocVector(VECSXP, 4));
    SEXP fmt_names = PROTECT(allocVector(STRSXP, 4));
    SET_STRING_ELT(fmt_names, 0, mkChar("id"));
    SET_STRING_ELT(fmt_names, 1, mkChar("number"));
    SET_STRING_ELT(fmt_names, 2, mkChar("type"));
    SET_STRING_ELT(fmt_names, 3, mkChar("description"));
    SEXP fid_col = PROTECT(allocVector(STRSXP, meta->n_format));
    SEXP fnum_col = PROTECT(allocVector(STRSXP, meta->n_format));
    SEXP ftype_col = PROTECT(allocVector(STRSXP, meta->n_format));
    SEXP fdesc_col = PROTECT(allocVector(STRSXP, meta->n_format));
    for (int i = 0; i < meta->n_format; i++) {
        SET_STRING_ELT(fid_col, i, mkChar(meta->format_fields[i].id ? meta->format_fields[i].id : ""));
        SET_STRING_ELT(fnum_col, i, mkChar(meta->format_fields[i].number ? meta->format_fields[i].number : ""));
        SET_STRING_ELT(ftype_col, i, mkChar(
            meta->format_fields[i].type == VCF_TYPE_INTEGER ? "Integer" :
            meta->format_fields[i].type == VCF_TYPE_FLOAT ? "Float" :
            meta->format_fields[i].type == VCF_TYPE_FLAG ? "Flag" :
            meta->format_fields[i].type == VCF_TYPE_CHARACTER ? "Character" : "String"
        ));
        SET_STRING_ELT(fdesc_col, i, mkChar(meta->format_fields[i].description ? meta->format_fields[i].description : ""));
    }
    SET_VECTOR_ELT(fmt_df, 0, fid_col);
    SET_VECTOR_ELT(fmt_df, 1, fnum_col);
    SET_VECTOR_ELT(fmt_df, 2, ftype_col);
    SET_VECTOR_ELT(fmt_df, 3, fdesc_col);
    setAttrib(fmt_df, R_NamesSymbol, fmt_names);
    SET_VECTOR_ELT(out, 2, fmt_df);
    // filter_fields
    SEXP filter_df = PROTECT(allocVector(VECSXP, 2));
    SEXP filter_names = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(filter_names, 0, mkChar("id"));
    SET_STRING_ELT(filter_names, 1, mkChar("description"));
    SEXP fid_col2 = PROTECT(allocVector(STRSXP, meta->n_filter));
    SEXP fdesc_col2 = PROTECT(allocVector(STRSXP, meta->n_filter));
    for (int i = 0; i < meta->n_filter; i++) {
        SET_STRING_ELT(fid_col2, i, mkChar(meta->filter_fields[i].id ? meta->filter_fields[i].id : ""));
        SET_STRING_ELT(fdesc_col2, i, mkChar(meta->filter_fields[i].description ? meta->filter_fields[i].description : ""));
    }
    SET_VECTOR_ELT(filter_df, 0, fid_col2);
    SET_VECTOR_ELT(filter_df, 1, fdesc_col2);
    setAttrib(filter_df, R_NamesSymbol, filter_names);
    SET_VECTOR_ELT(out, 3, filter_df);
    // alt_fields
    SEXP alt_df = PROTECT(allocVector(VECSXP, 2));
    SEXP alt_names = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(alt_names, 0, mkChar("id"));
    SET_STRING_ELT(alt_names, 1, mkChar("description"));
    SEXP aid_col = PROTECT(allocVector(STRSXP, meta->n_alt));
    SEXP adesc_col = PROTECT(allocVector(STRSXP, meta->n_alt));
    for (int i = 0; i < meta->n_alt; i++) {
        SET_STRING_ELT(aid_col, i, mkChar(meta->alt_fields[i].id ? meta->alt_fields[i].id : ""));
        SET_STRING_ELT(adesc_col, i, mkChar(meta->alt_fields[i].description ? meta->alt_fields[i].description : ""));
    }
    SET_VECTOR_ELT(alt_df, 0, aid_col);
    SET_VECTOR_ELT(alt_df, 1, adesc_col);
    setAttrib(alt_df, R_NamesSymbol, alt_names);
    SET_VECTOR_ELT(out, 4, alt_df);
    // contig_fields
    SEXP contig_df = PROTECT(allocVector(VECSXP, 4));
    SEXP contig_names = PROTECT(allocVector(STRSXP, 4));
    SET_STRING_ELT(contig_names, 0, mkChar("id"));
    SET_STRING_ELT(contig_names, 1, mkChar("length"));
    SET_STRING_ELT(contig_names, 2, mkChar("md5"));
    SET_STRING_ELT(contig_names, 3, mkChar("url"));
    SEXP cid_col = PROTECT(allocVector(STRSXP, meta->n_contig));
    SEXP clen_col = PROTECT(allocVector(STRSXP, meta->n_contig));
    SEXP cmd5_col = PROTECT(allocVector(STRSXP, meta->n_contig));
    SEXP curl_col = PROTECT(allocVector(STRSXP, meta->n_contig));
    for (int i = 0; i < meta->n_contig; i++) {
        SET_STRING_ELT(cid_col, i, mkChar(meta->contig_fields[i].id ? meta->contig_fields[i].id : ""));
        SET_STRING_ELT(clen_col, i, mkChar(meta->contig_fields[i].length ? meta->contig_fields[i].length : ""));
        SET_STRING_ELT(cmd5_col, i, mkChar(meta->contig_fields[i].md5 ? meta->contig_fields[i].md5 : ""));
        SET_STRING_ELT(curl_col, i, mkChar(meta->contig_fields[i].url ? meta->contig_fields[i].url : ""));
    }
    SET_VECTOR_ELT(contig_df, 0, cid_col);
    SET_VECTOR_ELT(contig_df, 1, clen_col);
    SET_VECTOR_ELT(contig_df, 2, cmd5_col);
    SET_VECTOR_ELT(contig_df, 3, curl_col);
    setAttrib(contig_df, R_NamesSymbol, contig_names);
    SET_VECTOR_ELT(out, 5, contig_df);
    // sample_names
    SEXP samples = PROTECT(allocVector(STRSXP, meta->n_samples));
    for (int i = 0; i < meta->n_samples; i++)
        SET_STRING_ELT(samples, i, mkChar(meta->sample_names[i] ? meta->sample_names[i] : ""));
    SET_VECTOR_ELT(out, 6, samples);
    setAttrib(out, R_NamesSymbol, names);
    UNPROTECT(12);
    return out;
}

// Query by contiguous indices using existing context
SEXP RC_VBI_query_by_indices_ctx(SEXP vbi_vcf_ctx, SEXP start_idx, SEXP end_idx, SEXP include_info, SEXP include_format, SEXP include_genotypes) {
    VBIVcfContextPtr ctx = (VBIVcfContextPtr) R_ExternalPtrAddr(vbi_vcf_ctx);
    if (!ctx) Rf_error("[VBI] VCF context pointer is NULL");
    if (!ctx->vbi_idx) Rf_error("[VBI] No VBI index available in context");
    int inc_info = asLogical(include_info);
    int inc_format = asLogical(include_format);
    int inc_geno = asLogical(include_genotypes);
    int start = asInteger(start_idx) - 1;
    int end = asInteger(end_idx) - 1;
    if (start < 0) start = 0;
    if (end >= ctx->vbi_idx->num_marker) end = ctx->vbi_idx->num_marker - 1;
    if (end < start) {
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
    int *hits = (int*) malloc(sizeof(int)*nfound); if(!hits) Rf_error("[VBI] OOM");
    for(int i=0;i<nfound;i++) hits[i]= start + i;
    SEXP out = vbi_query_variants_basic(ctx, hits, nfound, inc_info, inc_format, inc_geno);
    free(hits);
    return out;
}

// CGRanges-optimized region query using existing VBI VCF context
SEXP RC_VBI_query_region_cgranges_ctx(SEXP vbi_vcf_ctx, SEXP region_str, SEXP include_info, SEXP include_format, SEXP include_genotypes) {
    VBIVcfContextPtr ctx = (VBIVcfContextPtr) R_ExternalPtrAddr(vbi_vcf_ctx);
    if (!ctx) Rf_error("[VBI] VCF context pointer is NULL");
    if (!ctx->vbi_idx) Rf_error("[VBI] No VBI index available in context");
    
    int inc_info = asLogical(include_info);
    int inc_format = asLogical(include_format);
    int inc_genotypes = asLogical(include_genotypes);
    const char *reg = CHAR(STRING_ELT(region_str, 0));
    
    // Use cgranges for fast interval tree query
    int nfound = 0;
    int *indices = vbi_index_query_region_cgranges(ctx->vbi_idx, reg, &nfound);
    if (!indices || nfound == 0) {
        if (indices) free(indices);
        SEXP df = PROTECT(allocVector(VECSXP, 0));
        SEXP names = PROTECT(allocVector(STRSXP, 0));
        SEXP rn = PROTECT(allocVector(INTSXP, 0));
        setAttrib(df, R_NamesSymbol, names);
        setAttrib(df, R_RowNamesSymbol, rn);
        setAttrib(df, R_ClassSymbol, mkString("data.frame"));
        UNPROTECT(3);
        return df;
    }
    
    // Use existing VCF context - no need to reopen files!
    SEXP out = vbi_query_variants_basic(ctx, indices, nfound, inc_info, inc_format, inc_genotypes);
    free(indices);
    return out;
}


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