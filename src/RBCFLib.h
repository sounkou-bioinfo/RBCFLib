#ifndef RBCFLIB_H
#define RBCFLIB_H
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* Binary path function */
const char* BCFToolsBinaryPath(void);

/* Version function declarations */
SEXP RC_HTSLibVersion(void);
SEXP RC_BCFToolsVersion(void);

/* Unified bcftools pipeline function */
#ifndef _WIN32

SEXP RC_bcftools_pipeline(SEXP commands, SEXP args, SEXP n_commands,
                    SEXP capture_stdout, SEXP capture_stderr, SEXP stdout_file, SEXP stderr_file);
#endif
/* FASTA index and retrieval functions */
SEXP RC_FaidxIndexFasta(SEXP fasta_path);
SEXP RC_FaidxFetchRegion(SEXP fasta_path, SEXP seqname, SEXP start, SEXP end);

/* VBI index and query functions */
extern SEXP RC_VBI_index(SEXP vcf_path, SEXP vbi_path, SEXP threads);
extern SEXP RC_VBI_query_range(SEXP vcf_path, SEXP vbi_path, SEXP region, SEXP threads);
extern SEXP RC_VBI_query_index(SEXP vcf_path, SEXP vbi_path, SEXP start_idx, SEXP end_idx, SEXP threads);
extern SEXP RC_VBI_print_index(SEXP vbi_path, SEXP n);
SEXP RC_VBI_load_index(SEXP vbi_path);
SEXP RC_VBI_query_region_cgranges(SEXP idx_ptr, SEXP region_str);

#endif /* RBCFLIB_H */