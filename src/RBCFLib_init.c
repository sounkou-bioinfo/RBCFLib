#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>
#include "RBCFLib.h"

/* R native routine registration */
static const R_CallMethodDef CallEntries[] = {
    {"RC_VBI_extract_ranges", (DL_FUNC) &RC_VBI_extract_ranges, 2},
    {"RC_cgranges_create", (DL_FUNC) &RC_cgranges_create, 0},
    {"RC_cgranges_add", (DL_FUNC) &RC_cgranges_add, 5},
    {"RC_cgranges_index", (DL_FUNC) &RC_cgranges_index, 1},
    {"RC_cgranges_overlap", (DL_FUNC) &RC_cgranges_overlap, 4},
    {"RC_cgranges_destroy", (DL_FUNC) &RC_cgranges_destroy, 1},
    {"RC_VBI_load_index", (DL_FUNC) &RC_VBI_load_index, 1},
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
    {"RC_VBI_query_by_indices", (DL_FUNC) &RC_VBI_query_by_indices, 5},
    {"RC_VBI_print_index", (DL_FUNC) &RC_VBI_print_index, 2},
    {"RC_VBI_query_region_cgranges", (DL_FUNC) &RC_VBI_query_region_cgranges, 3},
    {NULL, NULL, 0}
};

void R_init_RBCFLib(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}