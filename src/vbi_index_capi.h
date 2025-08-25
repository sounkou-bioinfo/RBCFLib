#ifndef VBI_INDEX_CAPI_H
#define VBI_INDEX_CAPI_H

#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>
#include "htslib/hts.h"
#include "htslib/bgzf.h"
#include "htslib/vcf.h"
#include "htslib/hfile.h"
#include "cgranges.h"
#include <stdint.h>
#include <stddef.h>
#include <Rinternals.h>

typedef struct {
    int64_t num_sample;
    int64_t num_marker;
    int32_t *chrom_ids;
    int64_t *positions;
    int64_t *offsets;
    char **chrom_names;
    int n_chroms;
    void *cr; // cgranges_t *cr; (opaque)
} vbi_index_t;

typedef struct {
    char chrom[128];
    int64_t start;
    int64_t end;
    int is_point;
} region_t;

vbi_index_t *vbi_index_load(const char *filename);
void vbi_index_free(vbi_index_t *idx);
int parse_regions(const char *str, region_t **regions, int *nregions);


// Load VBI index from file
vbi_index_t *vbi_index_load(const char *filename);
void vbi_index_free(vbi_index_t *idx);

// Query by region string (e.g., "chr1:1000-2000"). Returns array of variant indices (0-based), sets *nfound.
int *vbi_index_query_region(vbi_index_t *idx, const char *region_str, int *nfound);

// Query by index range (nth to kth, inclusive). Returns array of indices, sets *nfound.
int *vbi_index_query_index_range(vbi_index_t *idx, int start, int end, int *nfound);

// Get BGZF offsets for given variant indices (array of int), returns array of int64_t offsets.
int64_t *vbi_index_offsets_for_indices(vbi_index_t *idx, const int *indices, int n);

// Get chromosome name for a given variant index
const char *vbi_index_chrom_name(vbi_index_t *idx, int idx_var);
// Get position for a given variant index
int64_t vbi_index_position(vbi_index_t *idx, int idx_var);

// Print the first n lines of the VBI index for debugging
void vbi_index_print(const vbi_index_t *idx, int n);

int *vbi_index_query_region_cgranges(vbi_index_t *idx, const char *region_str, int *nfound);

#endif // VBI_INDEX_CAPI_H
