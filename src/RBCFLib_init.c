#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>
#include "RBCFLib.h"


extern SEXP RC_bcftools_run(SEXP command, SEXP args, SEXP capture_stdout, SEXP capture_stderr, 
                    SEXP stdout_file, SEXP stderr_file, SEXP is_usage);
/* R native routine registration */
static const R_CallMethodDef CallEntries[] = {
    {"RC_HTSLibVersion", (DL_FUNC) &RC_HTSLibVersion, 0},
    {"RC_BCFToolsVersion", (DL_FUNC) &RC_BCFToolsVersion, 0},
    {"RC_bcftools_run", (DL_FUNC) &RC_bcftools_run, 7},
    {NULL, NULL, 0}
};

void R_init_RBCFLib(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}