#ifndef RBCFLIB_H
#define RBCFLIB_H
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
/* Version function declarations */
SEXP RC_HTSLibVersion(void);
SEXP RC_BCFToolsVersion(void);
SEXP RC_bcftools_run(SEXP command, SEXP args, SEXP capture_stdout, SEXP capture_stderr, 
                    SEXP stdout_file, SEXP stderr_file);
SEXP RC_bcftools_munge(SEXP args, SEXP capture_stdout, SEXP capture_stderr, SEXP stdout_file,
                  SEXP stderr_file);
SEXP RC_bcftools_score(SEXP args, SEXP capture_stdout, SEXP capture_stderr, SEXP stdout_file,
                  SEXP stderr_file);
SEXP RC_bcftools_liftover(SEXP args, SEXP capture_stdout, SEXP capture_stderr, SEXP stdout_file,
                  SEXP stderr_file);
SEXP RC_bcftools_metal(SEXP args, SEXP capture_stdout, SEXP capture_stderr, SEXP stdout_file,
                  SEXP stderr_file);
SEXP RC_bcftools_pgs(SEXP args, SEXP capture_stdout, SEXP capture_stderr, SEXP stdout_file,
                  SEXP stderr_file);
SEXP RC_bcftools_blup(SEXP args, SEXP capture_stdout, SEXP capture_stderr, SEXP stdout_file,
                  SEXP stderr_file);

/* FASTA index and retrieval functions */
SEXP RC_FaidxIndexFasta(SEXP fasta_path);
SEXP RC_FaidxFetchRegion(SEXP fasta_path, SEXP seqname, SEXP start, SEXP end);

/* ALTBCF functions */
SEXP RC_InitALTBCF(SEXP max_memory_mb);
SEXP RC_CreateALTBCFIntColumn(SEXP filename, SEXP field_name, SEXP region);
SEXP RC_CreateALTBCFStrColumn(SEXP filename, SEXP field_name, SEXP region);
SEXP RC_GetALTBCFStats(void);
SEXP RC_GetDetailedCacheStats(void);
SEXP RC_ResetCacheStats(void);
SEXP RC_PrefetchBlocks(SEXP column, SEXP block_ids);
SEXP RC_SetCacheParams(SEXP max_memory_mb, SEXP prefetch_size);

#endif /* RBCFLIB_H */