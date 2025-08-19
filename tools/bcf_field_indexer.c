/* bcf_field_indexer.c
 *
 * Build field-specific indexes for BCF files, storing offsets for each record and for specific fields (e.g., INFO, FORMAT, genotype data) in a directory.
 *
 * Usage:
 *   ./bcf_field_indexer yourfile.bcf index_dir/
 *
 * This is a sketch for demonstration. Only INFO and FORMAT offsets are indexed for each record.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <inttypes.h>
#include "htslib/hts.h"
#include "htslib/bgzf.h"
#include "htslib/vcf.h"

#define FIELD_INFO 1
#define FIELD_FORMAT 2

static void ensure_dir(const char *dir) {
    struct stat st = {0};
    if (stat(dir, &st) == -1) {
        mkdir(dir, 0777);
    }
}

int main(int argc, char **argv) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s yourfile.bcf index_dir/\n", argv[0]);
        return 1;
    }
    const char *bcf_path = argv[1];
    const char *index_dir = argv[2];
    ensure_dir(index_dir);

    char info_index_path[4096], format_index_path[4096], record_index_path[4096];
    snprintf(info_index_path, sizeof(info_index_path), "%s/info_offsets.bin", index_dir);
    snprintf(format_index_path, sizeof(format_index_path), "%s/format_offsets.bin", index_dir);
    snprintf(record_index_path, sizeof(record_index_path), "%s/record_offsets.bin", index_dir);

    FILE *info_idx = fopen(info_index_path, "wb");
    FILE *format_idx = fopen(format_index_path, "wb");
    FILE *record_idx = fopen(record_index_path, "wb");
    if (!info_idx || !format_idx || !record_idx) {
        fprintf(stderr, "Error: cannot open index files for writing\n");
        return 1;
    }

    htsFile *fp = hts_open(bcf_path, "r");
    if (!fp) { fprintf(stderr, "Error: cannot open %s\n", bcf_path); return 1; }
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (!hdr) { fprintf(stderr, "Error: failed to read BCF header\n"); hts_close(fp); return 1; }
    bcf1_t *rec = bcf_init();

    size_t nrec = 0;
    int64_t file_offset = 0;
    while (1) {
        file_offset = bgzf_tell((BGZF *)fp->fp.bgzf);
        int r = bcf_read(fp, hdr, rec);
        if (r < 0) break;
        // Write record offset
        fwrite(&file_offset, sizeof(file_offset), 1, record_idx);

        // BCF record layout: [CHROM][POS][rlen][QUAL][n_info][n_allele][n_fmt][n_sample] ...
        // We want to find the offset of INFO and FORMAT fields within rec->shared and rec->indiv
        // For simplicity, we use the offsets within the decompressed record buffer
        // INFO is at the start of shared (after fixed fields)
        // FORMAT is at the start of indiv
        // For a more precise index, parse the shared/indiv buffer
        int32_t info_offset = 0; // relative to rec->shared.s
        int32_t format_offset = 0; // relative to rec->indiv.s
        // Write offsets (for demo, just 0)
        fwrite(&info_offset, sizeof(info_offset), 1, info_idx);
        fwrite(&format_offset, sizeof(format_offset), 1, format_idx);
        nrec++;
    }
    printf("Indexed %zu records.\n", nrec);

    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    hts_close(fp);
    fclose(info_idx);
    fclose(format_idx);
    fclose(record_idx);
    return 0;
}
