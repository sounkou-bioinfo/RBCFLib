#ifndef RBCFLIB_H
#define RBCFLIB_H
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "bcftools.RBCFLIB.h"
/* Version function declarations */
SEXP RC_HTSLibVersion(void);
SEXP RC_BCFToolsVersion(void);
SEXP RC_bcftools_run(SEXP command, SEXP args, SEXP capture_stdout, SEXP capture_stderr, 
                    SEXP stdout_file, SEXP stderr_file, SEXP is_usage);

SEXP RC_bcftools_munge(SEXP args, SEXP capture_stdout, SEXP capture_stderr, 
                       SEXP stdout_file, SEXP stderr_file, SEXP is_usage);
#endif /* RBCFLIB_H */