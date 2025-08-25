#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>
#include "RBCFLib.h"

/* R native routine registration */
static const R_CallMethodDef CallEntries[] = {
    {"RC_HTSLibVersion", (DL_FUNC) &RC_HTSLibVersion, 0},
    {"RC_BCFToolsVersion", (DL_FUNC) &RC_BCFToolsVersion, 0},
    #ifndef _WIN32
    {"RC_bcftools_pipeline", (DL_FUNC) &RC_bcftools_pipeline, 7},
    #endif
    {"RC_FaidxIndexFasta", (DL_FUNC) &RC_FaidxIndexFasta, 1},
    {"RC_FaidxFetchRegion", (DL_FUNC) &RC_FaidxFetchRegion, 4},
    /* vbi*/
    {"RC_VBI_index", (DL_FUNC) &RC_VBI_index, 3},
    {"RC_VBI_query_range", (DL_FUNC) &RC_VBI_query_range, 4},
    {"RC_VBI_query_index", (DL_FUNC) &RC_VBI_query_index, 5},
    {NULL, NULL, 0}
};

void R_init_RBCFLib(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}