#include <stdint.h>
#include <limits.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include "vbi_index_capi.h"
#include <R_ext/Print.h>

// Parse a single region string (chr, chr:pos, chr:start-end)
int parse_region(const char *str, region_t *reg) {
    reg->is_point = 0;
    reg->start = 0; reg->end = 0;
    const char *colon = strchr(str, ':');
    if (!colon) {
        strncpy(reg->chrom, str, sizeof(reg->chrom)-1);
        reg->chrom[sizeof(reg->chrom)-1] = 0;
        reg->start = 0; reg->end = INT64_MAX;
        return 0;
    }
    size_t clen = colon-str;
    strncpy(reg->chrom, str, clen);
    reg->chrom[clen] = 0;
    const char *dash = strchr(colon+1, '-');
    if (!dash) {
        reg->start = atoll(colon+1);
        reg->end = reg->start;
        reg->is_point = 1;
    } else {
        reg->start = atoll(colon+1);
        reg->end = atoll(dash+1);
    }
    return 0;
}

int parse_regions(const char *str, region_t **regions, int *nregions) {
    int count = 1;
    for (const char *p = str; *p; ++p) if (*p == ',') count++;
    *regions = malloc(count * sizeof(region_t));
    *nregions = 0;
    char *tmp = strdup(str);
    char *saveptr = NULL;
    char *tok = strtok_r(tmp, ",", &saveptr);
    while (tok) {
        if (parse_region(tok, &(*regions)[*nregions]) == 0) (*nregions)++;
        tok = strtok_r(NULL, ",", &saveptr);
    }
    free(tmp);
    return 0;
}


// Free VBI index
void vbi_index_free(vbi_index_t *idx) {
    if (!idx) return;
    if (idx->chrom_ids) free(idx->chrom_ids);
    if (idx->positions) free(idx->positions);
    if (idx->offsets) free(idx->offsets);
    if (idx->chrom_names) {
        for (int i = 0; i < idx->n_chroms; ++i) free(idx->chrom_names[i]);
        free(idx->chrom_names);
    }
    if (idx->cr) cr_destroy(idx->cr);
    free(idx);
}


// Query by region string (returns malloc'd array of indices)
int *vbi_index_query_region(vbi_index_t *idx, const char *region_str, int *nfound) {
    // parse region string (support comma-separated)
    region_t *regions = NULL;
    int nregions = 0;
    if (parse_regions(region_str, &regions, &nregions) != 0) {
        *nfound = 0;
        return NULL;
    }
    int *hits = malloc(idx->num_marker * sizeof(int));
    int nhits = 0;
    for (size_t i = 0; i < idx->num_marker; ++i) {
        const char *chr = idx->chrom_names[idx->chrom_ids[i]];
        for (int r = 0; r < nregions; ++r) {
            if (strcmp(chr, regions[r].chrom) == 0 &&
                idx->positions[i] >= regions[r].start &&
                idx->positions[i] <= regions[r].end) {
                hits[nhits++] = i;
                break;
            }
        }
    }
    free(regions);
    *nfound = nhits;
    return hits;
}

// Query by index range (nth to kth, inclusive)
int *vbi_index_query_index_range(vbi_index_t *idx, int start, int end, int *nfound) {
    if (start < 0) start = 0;
    if (end >= idx->num_marker) end = idx->num_marker - 1;
    if (end < start) { *nfound = 0; return NULL; }
    int n = end - start + 1;
    int *arr = malloc(n * sizeof(int));
    for (int i = 0; i < n; ++i) arr[i] = start + i;
    *nfound = n;
    return arr;
}

// Get BGZF offsets for given variant indices
int64_t *vbi_index_offsets_for_indices(vbi_index_t *idx, const int *indices, int n) {
    int64_t *arr = malloc(n * sizeof(int64_t));
    for (int i = 0; i < n; ++i) arr[i] = idx->offsets[indices[i]];
    return arr;
}

const char *vbi_index_chrom_name(vbi_index_t *idx, int idx_var) {
    return idx->chrom_names[idx->chrom_ids[idx_var]];
}

int64_t vbi_index_position(vbi_index_t *idx, int idx_var) {
    return idx->positions[idx_var];
}


int do_index(const char *infile, const char *outfile, int n_threads) {
    htsFile *fp = NULL;
    if (n_threads > 1) {
        char optstr[64];
        snprintf(optstr, sizeof(optstr), "r:threads=%d", n_threads);
        fp = hts_open(infile, optstr);
    } else {
        fp = hts_open(infile, "r");
    }
    if (!fp) {
        Rf_error("Error: cannot open %s\n", infile);
        return 1;
    }
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (!hdr) {
        Rf_error("Error: failed to read VCF/BCF header\n");
        hts_close(fp);
        return 1;
    }
    int64_t num_sample = bcf_hdr_nsamples(hdr);
    int64_t num_marker = 0;
    int max_chroms = 64;
    char **chrom_names = malloc(max_chroms * sizeof(char*));
    int chrom_count = 0;
    int chrom_id = -1;
    int32_t *chrom_ids = NULL;
    int64_t *positions = NULL;
    int64_t *offsets = NULL;
    size_t alloc = 1024, n = 0;
    chrom_ids = malloc(alloc * sizeof(int32_t));
    positions = malloc(alloc * sizeof(int64_t));
    offsets = malloc(alloc * sizeof(int64_t));
    bcf1_t *rec = bcf_init();
    while (1) {
        BGZF *bg = (BGZF *)fp->fp.bgzf;
        int64_t this_offset = bgzf_tell(bg);
        int ret = bcf_read(fp, hdr, rec);
        if (ret != 0) break;
        const char *chr = bcf_hdr_id2name(hdr, rec->rid);
        // Find or add chrom
        int found = 0;
        for (int i = 0; i < chrom_count; ++i) {
            if (strcmp(chr, chrom_names[i]) == 0) { chrom_id = i; found = 1; break; }
        }
        if (!found) {
            if (chrom_count == max_chroms) {
                max_chroms *= 2;
                chrom_names = realloc(chrom_names, max_chroms * sizeof(char*));
            }
            chrom_names[chrom_count] = strdup(chr);
            chrom_id = chrom_count++;
        }
        if (n == alloc) {
            alloc *= 2;
            chrom_ids = realloc(chrom_ids, alloc * sizeof(int32_t));
            positions = realloc(positions, alloc * sizeof(int64_t));
            offsets = realloc(offsets, alloc * sizeof(int64_t));
        }
        chrom_ids[n] = chrom_id;
        positions[n] = rec->pos + 1;
        offsets[n] = this_offset;
        n++;
    }
    num_marker = n;
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    hts_close(fp);

    // Write index
    FILE *fidx = fopen(outfile, "wb");
    if (!fidx) { for (int i = 0; i < chrom_count; ++i) free(chrom_names[i]); free(chrom_names); free(chrom_ids); free(positions); free(offsets); return 1; }
    fwrite(&num_sample, sizeof(int64_t), 1, fidx);
    fwrite(&num_marker, sizeof(int64_t), 1, fidx);
    int32_t n_chroms32 = chrom_count;
    fwrite(&n_chroms32, sizeof(int32_t), 1, fidx);
    for (int i = 0; i < chrom_count; ++i) {
        int32_t len = strlen(chrom_names[i]);
        fwrite(&len, sizeof(int32_t), 1, fidx);
        fwrite(chrom_names[i], 1, len, fidx);
    }
    for (size_t i = 0; i < n; ++i) {
        fwrite(&chrom_ids[i], sizeof(int32_t), 1, fidx);
        fwrite(&positions[i], sizeof(int64_t), 1, fidx);
        fwrite(&offsets[i], sizeof(int64_t), 1, fidx);
    }
    fclose(fidx);
    for (int i = 0; i < chrom_count; ++i) free(chrom_names[i]);
    free(chrom_names); free(chrom_ids); free(positions); free(offsets);
    Rprintf("Indexing  finished: %" PRId64 " samples, %" PRId64 " markers, %d chromosomes\n", num_sample, num_marker, chrom_count);
    return 0;
}


// Load VBI index into memory
vbi_index_t *vbi_index_load(const char *filename) {
    FILE *f = fopen(filename, "rb");
    if (!f) return NULL;
    vbi_index_t *idx = calloc(1, sizeof(vbi_index_t));
    if (!idx) { fclose(f); return NULL; }
    if (fread(&idx->num_sample, sizeof(int64_t), 1, f) != 1) goto fail;
    if (fread(&idx->num_marker, sizeof(int64_t), 1, f) != 1) goto fail;
    int32_t n_chroms = 0;
    if (fread(&n_chroms, sizeof(int32_t), 1, f) != 1) goto fail;
    idx->n_chroms = n_chroms;
    idx->chrom_names = malloc(n_chroms * sizeof(char*));
    for (int i = 0; i < n_chroms; ++i) {
        int32_t len = 0;
        if (fread(&len, sizeof(int32_t), 1, f) != 1) goto fail;
        idx->chrom_names[i] = malloc(len+1);
        if (fread(idx->chrom_names[i], 1, len, f) != (size_t)len) goto fail;
        idx->chrom_names[i][len] = 0;
    }
    size_t n = idx->num_marker;
    idx->chrom_ids = malloc(n * sizeof(int32_t));
    idx->positions = malloc(n * sizeof(int64_t));
    idx->offsets = malloc(n * sizeof(int64_t));
    if (!idx->chrom_ids || !idx->positions || !idx->offsets) goto fail;
    for (size_t i = 0; i < n; ++i) {
        if (fread(&idx->chrom_ids[i], sizeof(int32_t), 1, f) != 1) goto fail;
        if (fread(&idx->positions[i], sizeof(int64_t), 1, f) != 1) goto fail;
        if (fread(&idx->offsets[i], sizeof(int64_t), 1, f) != 1) goto fail;
    }
    fclose(f);
    // Build cgranges for fast overlap
    idx->cr = cr_init();
    for (size_t i = 0; i < n; ++i) {
        cr_add(idx->cr, idx->chrom_names[idx->chrom_ids[i]], (int32_t)idx->positions[i], (int32_t)idx->positions[i], i);
    }
    cr_index(idx->cr);
    return idx;
fail:
    if (idx->chrom_ids) free(idx->chrom_ids);
    if (idx->positions) free(idx->positions);
    if (idx->offsets) free(idx->offsets);
    if (idx->chrom_names) {
        for (int i = 0; i < idx->n_chroms; ++i) free(idx->chrom_names[i]);
        free(idx->chrom_names);
    }
    free(idx);
    fclose(f);
    return NULL;
}


// Print the first n lines of the VBI index for debugging
void vbi_index_print(const vbi_index_t *idx, int n) {
    if (!idx) {
        Rprintf("[VBI] Index is NULL\n");
        return;
    }
    int max = n > 0 ? n : idx->num_marker;
    if (max > idx->num_marker) max = idx->num_marker;
    Rprintf("[VBI] num_marker: %ld\n", (long)idx->num_marker);
    for (int i = 0; i < max; ++i) {
        const char *chr = idx->chrom_names[idx->chrom_ids[i]];
        Rprintf("[VBI] %d: %s\t%ld\toffset=%ld\n", i+1, chr, (long)idx->positions[i], (long)idx->offsets[i]);
    }
    return ;
}