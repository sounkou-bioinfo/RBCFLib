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

#endif /* RBCFLIB_H */