/* test_mmap.c
 *
 * Build example:
 * gcc -O2 -Wall -std=c11 test_mmap.c -I/path/to/htslib/include -L/path/to/htslib/lib -lhts -lz -lbz2 -llzma -lcrypt -lpthread -o test_mmap
 *
 * Or (if you link with a static libhts.a as you did):
 * gcc -O2 -Wall -std=c11 test_mmap.c -I../src/bcftools-1.22/htslib-1.22 -I/usr/local/include -L../src/bcftools-1.22/htslib-1.22 -lhts -lz -lbz2 -llzma -lcurl -lcrypto -lssl -lpthread -o test_mmap
 *
 * Usage:
 * ./test_mmap yourfile.vcf.gz
 *
 */

#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <time.h>
#include <unistd.h>
#include "htslib/hfile.h"
#include "htslib/hts.h"
#include "htslib/bgzf.h"
#include "htslib/vcf.h"
#define SWAP(a,b) temp=a;a=b;b=temp;

float binmedian(int n, float *x, float *median, float *mean, float *std) {
  // Compute the mean and standard deviation
  float sum = 0;
  int i;
  for (i = 0; i < n; i++) {
    sum += x[i];
  }
  float mu = sum/n;
  *mean = mu;

  sum = 0;
  for (i = 0; i < n; i++) {
    sum += (x[i]-mu)*(x[i]-mu);
  }
  float sigma = sqrt(sum/n);
  *std = sigma;

  // Bin x across the interval [mu-sigma, mu+sigma]
  int bottomcount = 0;
  int bincounts[1001];
  for (i = 0; i < 1001; i++) {
    bincounts[i] = 0;
  }
  float scalefactor = 1000/(2*sigma);
  float leftend =  mu-sigma;
  float rightend = mu+sigma;
  int bin;

  for (i = 0; i < n; i++) {
    if (x[i] < leftend) {
      bottomcount++;
    }
    else if (x[i] < rightend) {
      bin = (int)((x[i]-leftend) * scalefactor);
      bincounts[bin]++;
    }
  }

  // If n is odd
  if (n & 1) {
    // Recursive step
    int k, r, count, medbin;
    float oldscalefactor, oldleftend;
    int oldbin;
    float temp;

    k = (n+1)/2;
    r = 0;

    for (;;) {
      // Find the bin that contains the median, and the order
      // of the median within that bin
      count = bottomcount;
      for (i = 0; i < 1001; i++) {
        count += bincounts[i];

        if (count >= k) {
          medbin = i;
          k = k - (count-bincounts[i]);
          break;
        }
      }

      bottomcount = 0;
      for (i = 0; i < 1001; i++) {
        bincounts[i] = 0;
      }
      oldscalefactor = scalefactor;
      oldleftend = leftend;
      scalefactor = 1000*oldscalefactor;
      leftend = medbin/oldscalefactor + oldleftend;
      rightend = (medbin+1)/oldscalefactor + oldleftend;

      // Determine which points map to medbin, and put
      // them in spots r,...n-1
      i = r; r = n;
      while (i < r) {
        oldbin = (int)((x[i]-oldleftend) * oldscalefactor);

        if (oldbin == medbin) {
          r--;
          SWAP(x[i],x[r]);

          // Re-bin on a finer scale
          if (x[i] < leftend) {
            bottomcount++;
          }
          else if (x[i] < rightend) {
            bin = (int)((x[i]-leftend) * scalefactor);
            bincounts[bin]++;
          }
        }
        else {
          i++;
        }
      }

      // Stop if all points in medbin are the same
      int samepoints = 1;
      for (i = r+1; i < n; i++) {
        if (x[i] != x[r]) {
          samepoints = 0;
          break;
        }
      }
      if (samepoints) {
        return x[r];
      }

      // Stop if there's <= 20 points left
      if (n-r <= 20) {
        break;
      }
    }

    // Perform insertion sort on the remaining points,
    // and then pick the kth smallest
    float a;
    int j;
    for (i = r+1; i < n; i++) {
      a = x[i];
      for (j = i-1; j >= r; j--) {
        if (x[j] < a) {
          break;
        }
        x[j+1] = x[j];
      }
      x[j+1] = a;
    }

    return x[r-1+k];
  }

  // If n is even (not implemented yet)
  else {
    return 0;
  }
}

/* optional: plugin initializer (no-op if not linked) */
extern int hfile_plugin_init_mmap(void);

static const char * compression_name(enum htsCompression c) {
    switch (c) {
        case no_compression: return "Uncompressed";
        case gzip:           return "GZIP";
        case bgzf:           return "BGZF";
        default:             return "Unknown";
    }
}

int main(int argc, char **argv) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <vcf[.gz]>\n", argv[0]);
        return 1;
    }
    const char *path = argv[1];
    char *ref;
    char *geno;
    size_t chunk_len = 1000;
    int ntimes =  101;

    /* init plugin registry (safe even if plugin isn't present) */
    hFILE *dummy = hopen("data:,", "r");
    if (dummy) { int _ = hclose(dummy); (void)_; }
    (void) hfile_plugin_init_mmap();

    /* open with mmap backend */
    char uri[4096];
    snprintf(uri, sizeof(uri), "mmap:%s", path);
    printf("Opening via mmap backend: %s\n", path);

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    htsFile *fp = hts_open(uri, "r");
    clock_gettime(CLOCK_MONOTONIC, &t1);
    if (!fp) {
        fprintf(stderr, "Error: cannot open %s\n", uri);
        return 1;
    }
    double open_s = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
    printf("✓ Opened. format: %s\n", hts_format_description(&fp->format));
    printf("Open time: %.6f s | compression: %s\n",
           open_s, compression_name(fp->format.compression));

    /* read header */
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (!hdr) {
        fprintf(stderr, "Error: failed to read VCF header\n");
        hts_close(fp);
        return 1;
    }
    printf("✓ Header loaded - %d samples\n", bcf_hdr_nsamples(hdr));

    /* prepare custom binary index */
    char bin_index_path[4096];
    snprintf(bin_index_path, sizeof(bin_index_path), "%s.customidx.bin", path);
    FILE *bout = fopen(bin_index_path, "wb");
    if (!bout) {
        fprintf(stderr, "Error: cannot open %s for writing\n", bin_index_path);
        bcf_hdr_destroy(hdr);
        hts_close(fp);
        return 1;
    }

    printf("\n--- Building custom index (start_voff, size) ---\n");
    printf("[INFO] For BGZF files offsets are virtual offsets from bgzf_tell(); use bgzf_seek(voff, SEEK_SET) to seek.\n");

    bcf1_t *rec = bcf_init();
    if (!rec) { fprintf(stderr,"Error: bcf_init failed\n"); fclose(bout); bcf_hdr_destroy(hdr); hts_close(fp); return 1; }

    size_t nrec = 0;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    const int is_bgzf = (fp->format.compression == bgzf);
    while (1) {
        int64_t start_off = -1, end_off = -1;

        if (is_bgzf) {
            BGZF *bg = (BGZF *)fp->fp.bgzf;
            start_off = bgzf_tell(bg);
        } else {
            hFILE *hf = (hFILE *)fp->fp.hfile;
            start_off = htell(hf);
        }

        int r = bcf_read(fp, hdr, rec);
        if (r < 0) break;

        if (is_bgzf) {
            BGZF *bg = (BGZF *)fp->fp.bgzf;
            end_off = bgzf_tell(bg);
        } else {
            hFILE *hf = (hFILE *)fp->fp.hfile;
            end_off = htell(hf);
        }

        int64_t off = start_off;
        int64_t sz = end_off - start_off;
        // Write offsets and sizes to the index file
        fwrite(&off, sizeof(int64_t), 1, bout);
        fwrite(&sz, sizeof(int64_t), 1, bout);
        bcf_unpack(rec, BCF_UN_STR);
        // do some fake work of extracting all alleles
      
       // print only every 100000 iterations
     //  if( nrec % 100000 == 0) {
      //     printf("Record %zu: off=%" PRId64 ", size=%" PRId64, nrec+1, off, sz);
           // Optionally print record info
       //    printf(" | %s:%" PRId64 " %s", bcf_seqname(hdr, rec), (int64_t)(rec->pos + 1), rec->d.allele[0]);
        //   printf("\n") ;
       //    if (rec->n_allele > 1) printf("->%s", rec->d.allele[1]);
       //} else {
       //   bcf_seqname(hdr, rec);
        //  ref = rec->d.allele[1];
       //}        

        nrec++;
    }

    fclose(bout);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    hts_close(fp);

    printf("\n✓ Done.\n");

    // --- Benchmark: Load index into memory ---
    printf("\n--- Benchmark: Load index into memory ---\n");
    FILE *idxf = fopen(bin_index_path, "rb");
    if (!idxf) { fprintf(stderr, "Error: cannot open %s for reading\n", bin_index_path); return 1; }
    int64_t *offsets = malloc(nrec * sizeof(int64_t));
    int64_t *sizes = malloc(nrec * sizeof(int64_t));
    if (!offsets || !sizes) { fprintf(stderr, "Error: malloc failed\n"); fclose(idxf); return 1; }
    for (size_t i = 0; i < nrec; ++i) {
        if (fread(&offsets[i], sizeof(int64_t), 1, idxf) != 1) break;
        if (fread(&sizes[i], sizeof(int64_t), 1, idxf) != 1) break;
    }
    fclose(idxf);
    printf("Index loaded: %zu records\n", nrec);

    // --- Benchmark: Random seek to Nth record ---
    printf("\n--- Benchmark: Random seek to Nth record ---\n");
    size_t nth = nrec / 2; // middle record
    clock_gettime(CLOCK_MONOTONIC, &t0);
    fp = hts_open(uri, "r");
    hdr = bcf_hdr_read(fp);
    rec = bcf_init();
    int seek_ok = -1;
    if (fp->format.compression == bgzf) {
        BGZF *bg = (BGZF *)fp->fp.bgzf;
        seek_ok = bgzf_seek(bg, offsets[nth], SEEK_SET);
    } else {
        hFILE *hf = (hFILE *)fp->fp.hfile;
        seek_ok = hseek(hf, (off_t)offsets[nth], SEEK_SET);
    }
    int r = (seek_ok == 0) ? bcf_read(fp, hdr, rec) : -1;
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double seek_s = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
    if (r >= 0) {
        bcf_unpack(rec, BCF_UN_STR);
        printf("Seeked to record %zu: %s:%" PRId64 " %s (%.6f ms)\n", nth+1, bcf_seqname(hdr, rec), (int64_t)(rec->pos+1), rec->d.allele[0], seek_s * 1000  );
    } else {
        printf("Seek/read failed for record %zu (%.6f ms)\n", nth+1, seek_s * 1000);
    }
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    hts_close(fp);

    // --- Benchmark: Access a chunk of consecutive records ---
    printf("\n--- Benchmark: Access a chunk of consecutive records ---\n");
    size_t chunk_start = nrec / 4;
  
    float chunks[ntimes] ;
    if (chunk_start + chunk_len > nrec) chunk_len = nrec - chunk_start;
    fp = hts_open(uri, "r");
    hdr = bcf_hdr_read(fp);
    rec = bcf_init();
    clock_gettime(CLOCK_MONOTONIC, &t0);
    for ( int n = 0 ; n < ntimes; n++) {
    for (size_t i = 0; i < chunk_len; ++i) {
        int seek_ok = -1;
        if (fp->format.compression == bgzf) {
            BGZF *bg = (BGZF *)fp->fp.bgzf;
            seek_ok = bgzf_seek(bg, offsets[chunk_start + i], SEEK_SET);
        } else {
            hFILE *hf = (hFILE *)fp->fp.hfile;
            seek_ok = hseek(hf, (off_t)offsets[chunk_start + i], SEEK_SET);
        }
        if (seek_ok == 0) { int _ = bcf_read(fp, hdr, rec); (void)_; }
    }
    clock_gettime(CLOCK_MONOTONIC, &t1);
    chunks[n] = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
}
float av = 0;
float min = 0;
float max = 0;
float std = 0;
float median = 0;
//   statistics
for (int n = 0; n < ntimes; n++){
    av += chunks[n];
    if( min == 0 || chunks[n] < min) min = chunks[n];
    if( max == 0 || chunks[n] > max) max = chunks[n];
}
av /= ntimes;
// compute std
// median
// sort the array
median = binmedian(ntimes,chunks,&median,&av,&std);

printf("Average: %.6f ms | Min: %.6f ms | Max: %.6f ms | Std: %.6f ms | Median: %.6f ms\n", av * 1000, min * 1000, max * 1000, std * 1000, median * 1000);
printf("Accessed %zu consecutive records (from %zu): %.6f ms (%.2f ms/record)\n", chunk_len, chunk_start+1, av * 1000, 1000.0*av/chunk_len);
bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    hts_close(fp);

    free(offsets);
    free(sizes);

    return 0;
}