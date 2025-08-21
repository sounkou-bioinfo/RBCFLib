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

#endif /* RBCFLIB_H */