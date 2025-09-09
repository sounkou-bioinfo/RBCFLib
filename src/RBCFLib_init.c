#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>
#include "RBCFLib.h"

/* R native routine registration */


static const R_CallMethodDef CallEntries[] = {
    /* CGRanges */
    {"RC_cgranges_create", (DL_FUNC) &RC_cgranges_create, 0},
    {"RC_cgranges_add", (DL_FUNC) &RC_cgranges_add, 5},
    {"RC_cgranges_index", (DL_FUNC) &RC_cgranges_index, 1},
    {"RC_cgranges_overlap", (DL_FUNC) &RC_cgranges_overlap, 4},
    {"RC_cgranges_destroy", (DL_FUNC) &RC_cgranges_destroy, 1},
    /* BCFTOOLS */
    {"RC_HTSLibVersion", (DL_FUNC) &RC_HTSLibVersion, 0},
    {"RC_BCFToolsVersion", (DL_FUNC) &RC_BCFToolsVersion, 0},
    #ifndef _WIN32
    {"RC_bcftools_pipeline", (DL_FUNC) &RC_bcftools_pipeline, 7},
    #endif
    /* FASTA */ 
    {"RC_FaidxIndexFasta", (DL_FUNC) &RC_FaidxIndexFasta, 1},
    {"RC_FaidxFetchRegion", (DL_FUNC) &RC_FaidxFetchRegion, 4},
    /* vbi*/
    {"RC_VBI_index", (DL_FUNC) &RC_VBI_index, 3},
    {"RC_VBI_query_range", (DL_FUNC) &RC_VBI_query_range, 5},
    {"RC_VBI_query_by_indices", (DL_FUNC) &RC_VBI_query_by_indices_ctx, 6},
    {"RC_VBI_print_index", (DL_FUNC) &RC_VBI_print_index, 2},
    {"RC_VBI_query_region_cgranges", (DL_FUNC) &RC_VBI_query_region_cgranges_ctx, 5},
    {"RC_cgranges_extract_by_index", (DL_FUNC) &RC_cgranges_extract_by_index, 2},
    {"RC_VBI_index_memory_usage", (DL_FUNC) &RC_VBI_index_memory_usage, 1},
    {"RC_VBI_extract_ranges", (DL_FUNC) &RC_VBI_extract_ranges, 2},
    {"RC_VBI_load_index", (DL_FUNC) &RC_VBI_load_index, 1},
    {"RC_VBI_vcf_load", (DL_FUNC) &RC_VBI_vcf_load, 2},
    {"RC_VBI_query_region", (DL_FUNC) &RC_VBI_query_region, 5},
    {"RC_VBI_samples", (DL_FUNC) &RC_VBI_samples, 1},
    {"RC_VBI_nsamples", (DL_FUNC) &RC_VBI_nsamples, 1},
    {"RC_VBI_sample_at", (DL_FUNC) &RC_VBI_sample_at, 2},
    {"RC_VBI_sample2index", (DL_FUNC) &RC_VBI_sample2index, 2},
    {"RC_VBI_infos", (DL_FUNC) &RC_VBI_infos, 1},
    {"RC_VBI_formats", (DL_FUNC) &RC_VBI_formats, 1},
    {"RC_VBI_filters", (DL_FUNC) &RC_VBI_filters, 1},
    /*rbcf functions*/

    {"RC_RBcfFileClose", (DL_FUNC) &RBcfFileClose, 1},
    {"RC_RBcfFileOpen", (DL_FUNC) &RBcfFileOpen, 2},
    {"RC_RBcfNewWriter", (DL_FUNC) &RBcfNewWriter, 2},
    {"RC_RBcfFileWriteCtx", (DL_FUNC) &RBcfFileWriteCtx, 2},
    {"RC_RBcfSeqNames", (DL_FUNC) &RBcfSeqNames, 1},
    {"RC_RBcfNSamples", (DL_FUNC) &RBcfNSamples, 1},
    {"RC_RBcfSamples", (DL_FUNC) &RBcfSamples, 1},
    {"RC_RBcfSampleAtIndex0", (DL_FUNC) &RBcfSampleAtIndex0, 2},
    {"RC_BcfFilterTable", (DL_FUNC) &BcfFilterTable, 1},
    {"RC_BcfInfoTable", (DL_FUNC) &BcfInfoTable, 1},
    {"RC_BcfFormatTable", (DL_FUNC) &BcfFormatTable, 1},
    {"RC_RBcfHeaderDict", (DL_FUNC) &RBcfHeaderDict, 1},
    {"RC_BcfConvertSampleToIndex0", (DL_FUNC) &BcfConvertSampleToIndex0, 2},
    {"RC_RBcfQueryRegion", (DL_FUNC) &RBcfQueryRegion, 2},
    {"RC_RBcfNextLine", (DL_FUNC) &RBcfNextLine, 1},
    {"RC_RBcfCtxRid", (DL_FUNC) &RBcfCtxRid, 1},
    {"RC_RBcfCtxSeqName", (DL_FUNC) &RBcfCtxSeqName, 1},
    {"RC_RBcfCtxPos", (DL_FUNC) &RBcfCtxPos, 1},
    {"RC_RBcfCtxHasId", (DL_FUNC) &RBcfCtxHasId, 1},
    {"RC_RBcfCtxId", (DL_FUNC) &RBcfCtxId, 1},
    {"RC_RBcfCtxEnd", (DL_FUNC) &RBcfCtxEnd, 1},
    {"RC_RBcfCtxNAlleles", (DL_FUNC) &RBcfCtxNAlleles, 1},
    {"RC_RBcfCtxAlleles", (DL_FUNC) &RBcfCtxAlleles, 1},
    {"RC_RBcfCtxReference", (DL_FUNC) &RBcfCtxReference, 1},
    {"RC_RBcfCtxAlternateAlleles", (DL_FUNC) &RBcfCtxAlternateAlleles, 1},
    {"RC_RBcfCtxHasQual", (DL_FUNC) &RBcfCtxHasQual, 1},
    {"RC_RBcfCtxQual", (DL_FUNC) &RBcfCtxQual, 1},
    {"RC_RBcfCtxFiltered", (DL_FUNC) &RBcfCtxFiltered, 1},
    {"RC_RBcfCtxFilters", (DL_FUNC) &RBcfCtxFilters, 1},
    {"RC_VariantHasFilter", (DL_FUNC) &VariantHasFilter, 2},
    {"RC_RBcfCtxVariantTypes", (DL_FUNC) &RBcfCtxVariantTypes, 1},
    {"RC_RBcfCtxVariantIsSnp", (DL_FUNC) &RBcfCtxVariantIsSnp, 1},
    {"RC_VariantNSamples", (DL_FUNC) &VariantNSamples, 1},
    {"RC_RBcfCtxVariantMaxPloidy", (DL_FUNC) &RBcfCtxVariantMaxPloidy, 1},
    {"RC_VariantGetGenotype", (DL_FUNC) &VariantGetGenotype, 2},
    {"RC_RBcfCtxVariantGtAllelesIndexes0", (DL_FUNC) &RBcfCtxVariantGtAllelesIndexes0, 1},
    {"RC_RBcfCtxVariantAllGtAllelesIndexes0", (DL_FUNC) &RBcfCtxVariantAllGtAllelesIndexes0, 1},
    {"RC_RBcfCtxVariantAllGtAllelesAlleleCounts", (DL_FUNC) &RBcfCtxVariantAllGtAllelesAlleleCounts, 2},
    {"RC_RBcfCtxVariantAllGtStrings", (DL_FUNC) &RBcfCtxVariantAllGtStrings, 1},
    {"RC_GenotypeSample", (DL_FUNC) &GenotypeSample, 1},
    {"RC_RBcfCtxVariantGtPhased", (DL_FUNC) &RBcfCtxVariantGtPhased, 1},
    {"RC_VariantHasAttribute", (DL_FUNC) &VariantHasAttribute, 2},
    {"RC_VariantGetInfoKeySet", (DL_FUNC) &VariantGetInfoKeySet, 1},
    {"RC_VariantGetFormatKeySet", (DL_FUNC) &VariantGetFormatKeySet, 1},
    {"RC_VariantVepTable", (DL_FUNC) &VariantVepTable, 1},
    {"RC_VariantSnpEffTable", (DL_FUNC) &VariantSnpEffTable, 1},
    {"RC_VariantStringAttribute", (DL_FUNC) &VariantStringAttribute, 2},
    {"RC_VariantIntAttribute", (DL_FUNC) &VariantIntAttribute, 2},
    {"RC_VariantFloatAttribute", (DL_FUNC) &VariantFloatAttribute, 2},
    {"RC_VariantFlagAttribute", (DL_FUNC) &VariantFlagAttribute, 2},
    {"RC_GenotypeStringAttribute", (DL_FUNC) &GenotypeStringAttribute, 2},
    {"RC_GenotypeInt32Attribute", (DL_FUNC) &GenotypeInt32Attribute, 2},
    {"RC_GenotypeFloatAttribute", (DL_FUNC) &GenotypeFloatAttribute, 2},
    {"RC_VariantGenotypesFlagAttribute", (DL_FUNC) &VariantGenotypesFlagAttribute, 2},
    {"RC_VariantGenotypesInt32Attribute", (DL_FUNC) &VariantGenotypesInt32Attribute, 2},
    {"RC_VariantGenotypesFloatAttribute", (DL_FUNC) &VariantGenotypesFloatAttribute, 2},
    /* end */
    {NULL, NULL, 0}
};

void R_init_RBCFLib(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}