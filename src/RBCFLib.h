#ifndef RBCFLIB_H
#define RBCFLIB_H

#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "htslib/hts.h"
#include "htslib/faidx.h"


/* 
    * version
*/
// C level
extern char *bcftools_version(void);
extern char *bcftools_score_version(void);
//  R Linkage
extern SEXP RC_HTSLibVersion(void);
extern SEXP RC_BCFToolsVersion(void);
extern SEXP RC_BCFToolsScoreVersion(void);

/* 
     * faidx
*/
extern SEXP RC_FaidxIndexFasta(SEXP fasta_path);
extern SEXP RC_FaidxFetchRegion(SEXP fasta_path, SEXP seqname, SEXP start, SEXP end);

/*

  * BCFTools Wrapper Functions

*/

/* Binary path storage */

/* Binary path storage and env vars*/
extern char *cached_bcftools_path;
extern char *cached_bcftools_plugins_path;


/* Binary path function */

const char* BCFToolsBinaryPath(void);
/* PLUGINS PATH */

const char* BCFToolsPluginsPath(void);

/* Unified bcftools pipeline function */

extern SEXP RC_bcftools_pipeline(SEXP commands, SEXP args, SEXP n_commands,
                    SEXP capture_stdout, SEXP capture_stderr, SEXP stdout_file, SEXP stderr_file);

/* 

 * VBI index and query functions

*/
extern SEXP RC_VBI_index(SEXP vcf_path, SEXP vbi_path, SEXP threads);
extern SEXP RC_VBI_query_range(SEXP vbi_vcf_ctx, SEXP region_str, SEXP include_info, SEXP include_format, SEXP include_genotypes);
extern SEXP RC_VBI_query_by_indices(SEXP vbi_vcf_ctx, SEXP start_idx, SEXP end_idx, SEXP include_info, SEXP include_format, SEXP include_genotypes);
extern SEXP RC_VBI_query_by_indices_ctx(SEXP vbi_vcf_ctx, SEXP start_idx, SEXP end_idx, SEXP include_info, SEXP include_format, SEXP include_genotypes);
extern SEXP RC_VBI_print_index(SEXP vbi_ptr, SEXP n);
extern SEXP RC_VBI_load_index(SEXP vbi_path);
extern SEXP RC_VBI_vcf_load(SEXP vcf_path, SEXP vbi_path);
extern SEXP RC_VBI_query_region(SEXP vbi_vcf_ctx, SEXP region_str, SEXP include_info, SEXP include_format, SEXP include_genotypes);
extern SEXP RC_VBI_query_region_cgranges_ctx(SEXP vbi_vcf_ctx, SEXP region_str, SEXP include_info, SEXP include_format, SEXP include_genotypes);
extern SEXP RC_VBI_samples(SEXP vbi_vcf_ctx);
extern SEXP RC_VBI_nsamples(SEXP vbi_vcf_ctx);
extern SEXP RC_VBI_sample_at(SEXP vbi_vcf_ctx, SEXP index);
extern SEXP RC_VBI_sample2index(SEXP vbi_vcf_ctx, SEXP sample_name);
extern SEXP RC_VBI_infos(SEXP vbi_vcf_ctx);
extern SEXP RC_VBI_formats(SEXP vbi_vcf_ctx);
extern SEXP RC_VBI_filters(SEXP vbi_vcf_ctx);
// VBI extract ranges
extern SEXP RC_VBI_extract_ranges(SEXP idx_ptr, SEXP n);
// cleanup
extern SEXP RC_VBI_index_memory_usage(SEXP extPtr);


/*

 * CGRanges functions

*/
extern SEXP RC_cgranges_create();
extern SEXP RC_cgranges_add(SEXP cr_ptr, SEXP chrom, SEXP start, SEXP end, SEXP label);
extern SEXP RC_cgranges_index(SEXP cr_ptr);
extern SEXP RC_cgranges_overlap(SEXP cr_ptr, SEXP chrom, SEXP start, SEXP end);
extern SEXP RC_cgranges_destroy(SEXP cr_ptr);
extern SEXP RC_cgranges_overlap(SEXP cr_ptr, SEXP chrom, SEXP start, SEXP end);
extern SEXP RC_cgranges_extract_by_index(SEXP cr_ptr, SEXP indices);

/*
 * rbcf R API functions
*/
extern SEXP RBcfFileClose(SEXP sexpFile);
extern SEXP RBcfFileOpen(SEXP Rfilename,SEXP sexpRequireIdx);
extern SEXP RBcfNewWriter(SEXP sexpIn,SEXP Rfilename);
extern SEXP RBcfFileWriteCtx(SEXP sexpOut,SEXP sexpCtx);
extern SEXP RBcfSeqNames(SEXP sexpFile);
extern SEXP RBcfNSamples(SEXP sexpFile);
extern SEXP RBcfSamples(SEXP sexpFile);
extern SEXP RBcfSampleAtIndex0(SEXP sexpFile,SEXP sexpIndex);
extern SEXP BcfFilterTable(SEXP sexpFile);
extern SEXP BcfInfoTable(SEXP sexpFile);
extern SEXP BcfFormatTable(SEXP sexpFile);
extern SEXP RBcfHeaderDict(SEXP sexpFile);
extern SEXP BcfConvertSampleToIndex0(SEXP sexpFile,SEXP sexpSample);
extern SEXP RBcfQueryRegion(SEXP sexpFile,SEXP sexpInterval);
extern SEXP RBcfNextLine(SEXP sexpFile);
extern SEXP RBcfCtxRid(SEXP sexpCtx);
extern SEXP RBcfCtxSeqName(SEXP sexpCtx);
extern SEXP RBcfCtxPos(SEXP sexpCtx);
extern SEXP RBcfCtxHasId(SEXP sexpCtx);
extern SEXP RBcfCtxId(SEXP sexpCtx);
extern SEXP RBcfCtxEnd(SEXP sexpCtx);
extern SEXP RBcfCtxNAlleles(SEXP sexpCtx);
extern SEXP RBcfCtxAlleles(SEXP sexpCtx);
extern SEXP RBcfCtxReference(SEXP sexpCtx);
extern SEXP RBcfCtxAlternateAlleles(SEXP sexpCtx);
extern SEXP RBcfCtxHasQual(SEXP sexpCtx);
extern SEXP RBcfCtxQual(SEXP sexpCtx);
extern SEXP RBcfCtxFiltered(SEXP sexpCtx);
extern SEXP RBcfCtxFilters(SEXP sexpCtx);
extern SEXP VariantHasFilter(SEXP sexpCtx,SEXP sexpFilterName);
extern SEXP RBcfCtxVariantTypes(SEXP sexpCtx);
extern SEXP RBcfCtxVariantIsSnp(SEXP sexpCtx);
extern SEXP VariantNSamples(SEXP sexpCtx);
extern SEXP RBcfCtxVariantMaxPloidy(SEXP sexpCtx);
extern SEXP VariantGetGenotype(SEXP sexpCtx,SEXP sexpgtidx);
extern SEXP RBcfCtxVariantGtAllelesIndexes0(SEXP sexpGt);
extern SEXP RBcfCtxVariantAllGtAllelesIndexes0(SEXP sexpCtx);
extern SEXP RBcfCtxVariantAllGtAllelesAlleleCounts(SEXP sexpCtx, SEXP sexpAlleleIndex);
extern SEXP RBcfCtxVariantAllGtStrings(SEXP sexpCtx);
extern SEXP GenotypeSample(SEXP sexpGt);
extern SEXP RBcfCtxVariantGtPhased(SEXP sexpGt);
extern SEXP VariantHasAttribute(SEXP sexpCtx,SEXP sexpatt);
extern SEXP VariantGetInfoKeySet(SEXP sexpCtx);
extern SEXP VariantGetFormatKeySet(SEXP sexpCtx);
extern SEXP VariantVepTable(SEXP sexpCtx);
extern SEXP VariantSnpEffTable(SEXP sexpCtx);
extern SEXP VariantStringAttribute(SEXP sexpCtx,SEXP sexpatt);
extern SEXP VariantIntAttribute(SEXP sexpCtx,SEXP sexpatt);
extern SEXP VariantFloatAttribute(SEXP sexpCtx,SEXP sexpatt);
extern SEXP VariantFlagAttribute(SEXP sexpCtx,SEXP sexpatt);
extern SEXP GenotypeStringAttribute(SEXP sexpGt,SEXP sexpatt);
extern SEXP GenotypeInt32Attribute(SEXP sexpGt,SEXP sexpatt);
extern SEXP GenotypeFloatAttribute(SEXP sexpGt,SEXP sexpatt);
extern SEXP VariantGenotypesFlagAttribute(SEXP sexpCtx,SEXP sexpatt);
extern SEXP VariantGenotypesInt32Attribute(SEXP sexpCtx,SEXP sexpatt);
extern SEXP VariantGenotypesFloatAttribute(SEXP sexpCtx,SEXP sexpatt); 
#endif /* RBCFLIB_H */