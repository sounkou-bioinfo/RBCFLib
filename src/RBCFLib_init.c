#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>
#include "RBCFLib.h"

/* R native routine registration */
static const R_CallMethodDef CallEntries[] = {
    {"RC_HTSLibVersion", (DL_FUNC) &RC_HTSLibVersion, 0},
    {"RC_BCFToolsVersion", (DL_FUNC) &RC_BCFToolsVersion, 0},
    {"RC_bcftools_run", (DL_FUNC) &RC_bcftools_run, 7},
    {"RC_bcftools_munge", (DL_FUNC) &RC_bcftools_munge, 6},
    {"RC_bcftools_score", (DL_FUNC) &RC_bcftools_score, 6},
    {"RC_bcftools_liftover", (DL_FUNC) &RC_bcftools_liftover, 6},
    {"RC_bcftools_metal", (DL_FUNC) &RC_bcftools_metal, 6},
    {"RC_bcftools_pgs", (DL_FUNC) &RC_bcftools_pgs, 6},
    {"RC_FaidxIndexFasta", (DL_FUNC) &RC_FaidxIndexFasta, 1},
    {"RC_FaidxFetchRegion", (DL_FUNC) &RC_FaidxFetchRegion, 4},
    {NULL, NULL, 0}
};

void R_init_RBCFLib(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}