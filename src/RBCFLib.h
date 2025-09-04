#ifndef RBCFLIB_H
#define RBCFLIB_H

#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "htslib/hts.h"
#include "htslib/faidx.h"


/* 
    * version
*/
// C level
extern char *bcftools_version(void);
//  R Linkage
extern SEXP RC_HTSLibVersion(void);
extern SEXP RC_BCFToolsVersion(void);

/* 
     * faidx
*/
extern SEXP RC_FaidxIndexFasta(SEXP fasta_path);
extern SEXP RC_FaidxFetchRegion(SEXP fasta_path, SEXP seqname, SEXP start, SEXP end);

/*

  * BCFTools Wrapper Functions

*/

/* Binary path storage */

/* Binary path storage and env vars*/
extern char *cached_bcftools_path;
extern char *cached_bcftools_plugins_path;


/* Binary path function */

const char* BCFToolsBinaryPath(void);
/* PLUGINS PATH */

const char* BCFToolsPluginsPath(void);

/* Unified bcftools pipeline function */

extern SEXP RC_bcftools_pipeline(SEXP commands, SEXP args, SEXP n_commands,
                    SEXP capture_stdout, SEXP capture_stderr, SEXP stdout_file, SEXP stderr_file);

/* 

 * VBI index and query functions

*/
extern SEXP RC_VBI_index(SEXP vcf_path, SEXP vbi_path, SEXP threads);
extern SEXP RC_VBI_query_range(SEXP vcf_path, SEXP vbi_path, SEXP region, SEXP threads);
extern SEXP RC_VBI_query_by_indices(SEXP vcf_path, SEXP vbi_path, SEXP start_idx, SEXP end_idx, SEXP threads);
extern SEXP RC_VBI_print_index(SEXP vbi_path, SEXP n);
extern SEXP RC_VBI_load_index(SEXP vbi_path);
extern SEXP RC_VBI_query_region_cgranges(SEXP idx_ptr, SEXP region_str);
// VBI extract ranges
extern SEXP RC_VBI_extract_ranges(SEXP idx_ptr, SEXP n);
// cleanup
extern SEXP RC_VBI_index_memory_usage(SEXP extPtr);


/*

 * CGRanges functions

*/
extern SEXP RC_cgranges_create();
extern SEXP RC_cgranges_add(SEXP cr_ptr, SEXP chrom, SEXP start, SEXP end, SEXP label);
extern SEXP RC_cgranges_index(SEXP cr_ptr);
extern SEXP RC_cgranges_overlap(SEXP cr_ptr, SEXP chrom, SEXP start, SEXP end);
extern SEXP RC_cgranges_destroy(SEXP cr_ptr);
extern SEXP RC_cgranges_overlap(SEXP cr_ptr, SEXP chrom, SEXP start, SEXP end);
extern SEXP RC_cgranges_extract_by_index(SEXP cr_ptr, SEXP indices);
#endif /* RBCFLIB_H */