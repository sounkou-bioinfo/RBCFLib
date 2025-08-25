/* vbi_index.c
 *
 * Create a VBI index for BCF or VCF.GZ files using htslib.
 * The index is a 2-column matrix of int64_t:
 *   line1:  num_sample  num_marker
 *   line2:  0           bgzf_offset_for_#CHROM_line
 *   line3:  var_1_pos   bgzf_offset_for_var_1
 *   ...
 * Usage: ./vbi_index input.bcf output.vbi
 */

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


void print_usage(const char *prog) {
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "  %s index <input.bcf|vcf.gz> <output.vbi>\n", prog);
    fprintf(stderr, "  %s query --vcf <input.bcf|vcf.gz> --vbi <index.vbi> [--threads N] <region1>[,<region2>...]\n", prog);
    fprintf(stderr, "    <region> format: chr, chr:pos, chr:start-end (e.g. 1, 1:1000, 1:1000-2000, 1:1000-2000,2:500-800)\n");
}



typedef struct {
    int64_t num_sample;
    int64_t num_marker;
    int32_t *chrom_ids;   // [num_marker]
    int64_t *positions;   // [num_marker]
    int64_t *offsets;     // [num_marker]
    char **chrom_names;   // [n_chroms]
    int n_chroms;
    cgranges_t *cr;       // for fast overlap
} vbi_index_t;

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

// Parse a single region string (chr, chr:pos, chr:start-end)
typedef struct {
    char chrom[128];
    int64_t start;
    int64_t end;
    int is_point; // 1 if single position
} region_t;

int parse_region(const char *str, region_t *reg) {
    // Accept: chr, chr:pos, chr:start-end
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

// Parse comma-separated regions
int parse_regions(const char *str, region_t **regions, int *nregions) {
    // Count commas
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

// Map chrom name to numeric ID using VCF header
int chrom_to_id(bcf_hdr_t *hdr, const char *chrom) {
    int id = bcf_hdr_name2id(hdr, chrom);
    return id;
}

// Binary search index for region
void vbi_find_range(const vbi_index_t *idx, int chrom_id, int64_t start, int64_t end, int *first, int *last) {
    // For single-chromosome index, just binary search positions
    int l = 0, r = idx->num_marker-1, m;
    *first = -1; *last = -1;
    // Find first >= start
    while (l <= r) {
        m = (l+r)/2;
        if (idx->positions[m] < start) l = m+1;
        else r = m-1;
    }
    if (l < idx->num_marker && idx->positions[l] <= end) *first = l;
    // Find last <= end
    l = 0; r = idx->num_marker-1;
    while (l <= r) {
        m = (l+r)/2;
        if (idx->positions[m] > end) r = m-1;
        else l = m+1;
    }
    if (r >= 0 && idx->positions[r] >= start) *last = r;
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
        fprintf(stderr, "Error: cannot open %s\n", infile);
        return 1;
    }
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (!hdr) {
        fprintf(stderr, "Error: failed to read VCF/BCF header\n");
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
    while (bcf_read(fp, hdr, rec) == 0) {
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
        BGZF *bg = (BGZF *)fp->fp.bgzf;
        offsets[n] = bgzf_tell(bg);
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
    printf("Indexing finished: %" PRId64 " samples, %" PRId64 " markers, %d chromosomes\n", num_sample, num_marker, chrom_count);
    return 0;
}

// Query VBI index for a range (chr:start-end)

// Query VBI index and VCF/BCF for regions
int do_query(const char *vcf_file, const char *vbi_file, const char *regions_str, int n_threads) {
    vbi_index_t *idx = vbi_index_load(vbi_file);
    if (!idx) { fprintf(stderr, "Error loading VBI index\n"); return 1; }

    htsFile *fp = NULL;
    if (n_threads > 1) {
        char optstr[64];
        snprintf(optstr, sizeof(optstr), "r:threads=%d", n_threads);
        fp = hts_open(vcf_file, optstr);
    } else {
        fp = hts_open(vcf_file, "r");
    }
    if (!fp) { fprintf(stderr, "Error: cannot open %s\n", vcf_file); vbi_index_free(idx); return 1; }
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (!hdr) { fprintf(stderr, "Error: failed to read VCF/BCF header\n"); hts_close(fp); vbi_index_free(idx); return 1; }

    region_t *regions = NULL;
    int nregions = 0;
    if (parse_regions(regions_str, &regions, &nregions) != 0) {
        fprintf(stderr, "Error parsing regions\n");
        bcf_hdr_destroy(hdr); hts_close(fp); vbi_index_free(idx); return 1;
    }

    // Build cgranges for query regions
    cgranges_t *crq = cr_init();
    for (int i = 0; i < nregions; ++i) {
        cr_add(crq, regions[i].chrom, (int32_t)regions[i].start, (int32_t)regions[i].end, i);
    }
    cr_index(crq);

    bcf1_t *rec = bcf_init();
    if (!rec) { fprintf(stderr, "Error: bcf_init failed\n"); cr_destroy(crq); bcf_hdr_destroy(hdr); hts_close(fp); vbi_index_free(idx); return 1; }

    // For each marker in index, check overlap
    int64_t *hits = NULL; int64_t m_hits = 0;
    htsFile *out = hts_open("-", "w"); // write to stdout
    if (!out) {
        fprintf(stderr, "Error: failed to open stdout for writing VCF/BCF\n");
        bcf_destroy(rec); cr_destroy(crq); free(regions); bcf_hdr_destroy(hdr); hts_close(fp); vbi_index_free(idx); return 1;
    }
    // Write header once
    if (bcf_hdr_write(out, hdr) < 0) {
        fprintf(stderr, "Error: failed to write VCF/BCF header\n");
        hts_close(out); bcf_destroy(rec); cr_destroy(crq); free(regions); bcf_hdr_destroy(hdr); hts_close(fp); vbi_index_free(idx); return 1;
    }
    for (size_t i = 0; i < idx->num_marker; ++i) {
        const char *chr = idx->chrom_names[idx->chrom_ids[i]];
        int64_t n = cr_overlap(crq, chr, (int32_t)idx->positions[i], (int32_t)idx->positions[i], &hits, &m_hits);
        if (n > 0) {
            // Seek and output
            int seek_ok = 0;
            if (fp->format.compression == bgzf) {
                BGZF *bg = (BGZF *)fp->fp.bgzf;
                seek_ok = (bgzf_seek(bg, idx->offsets[i], SEEK_SET) == 0);
            } else {
                hFILE *hf = (hFILE *)fp->fp.hfile;
                seek_ok = (hseek(hf, (off_t)idx->offsets[i], SEEK_SET) == 0);
            }
            if (!seek_ok) {
                fprintf(stderr, "Warning: seek failed for marker %zu\n", i);
                continue;
            }
            int r = bcf_read(fp, hdr, rec);
            if (r < 0) continue;
            bcf_unpack(rec, BCF_UN_STR);
            if (vcf_write1(out, hdr, rec) < 0) {
                fprintf(stderr, "Warning: vcf_write1 failed for marker %zu\n", i);
            }
        }
    }
    hts_close(out);
    free(hits);
    bcf_destroy(rec);
    cr_destroy(crq);
    free(regions);
    bcf_hdr_destroy(hdr);
    hts_close(fp);
    vbi_index_free(idx);
    return 0;
}


int main(int argc, char **argv) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }
    if (strcmp(argv[1], "index") == 0) {
        // Usage: index <input.bcf|vcf.gz> <output.vbi> [--threads N]
        const char *infile = NULL, *outfile = NULL;
        int n_threads = 1;
        int argi = 2;
        if (argi < argc) infile = argv[argi++];
        if (argi < argc) outfile = argv[argi++];
        while (argi < argc) {
            if (strcmp(argv[argi], "--threads") == 0 && argi+1 < argc) {
                n_threads = atoi(argv[++argi]);
            }
            argi++;
        }
        if (!infile || !outfile) {
            print_usage(argv[0]);
            return 1;
        }
        return do_index(infile, outfile, n_threads);
    } else if (strcmp(argv[1], "query") == 0) {
        // Usage: query --vcf <file> --vbi <index> [--threads N] <regions>
        const char *vcf = NULL, *vbi = NULL, *regions = NULL;
        int n_threads = 1;
        for (int i = 2; i < argc; ++i) {
            if (strcmp(argv[i], "--vcf") == 0 && i+1 < argc) { vcf = argv[++i]; }
            else if (strcmp(argv[i], "--vbi") == 0 && i+1 < argc) { vbi = argv[++i]; }
            else if (strcmp(argv[i], "--threads") == 0 && i+1 < argc) { n_threads = atoi(argv[++i]); }
            else if (!regions) { regions = argv[i]; }
        }
        if (!vcf || !vbi || !regions) {
            print_usage(argv[0]);
            return 1;
        }
        return do_query(vcf, vbi, regions, n_threads);
    } else {
        print_usage(argv[0]);
        return 1;
    }
}
