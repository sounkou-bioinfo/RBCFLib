#' Check if an object looks like a VCF context
#'
#' @name looks_like_vcf_context
#' @title Check if object is a VCF context
#' @param vcf the object to check
#' @return TRUE if object is a VCF context
#' @export
looks_like_vcf_context <- function(vcf) {
  is.list(vcf) &&
    length(vcf) == 3 &&
    class(vcf[[1]]) == "externalptr" &&
    class(vcf[[2]]) == "externalptr" &&
    is.character(vcf[[3]]) &&
    length(vcf[[3]]) == 1
}

#' Check if an object looks like a variant context
#'
#' @name looks_like_variant_context
#' @title Check if object is a variant context
#' @param vc the object to check
#' @return TRUE if object is a variant context
#' @export
looks_like_variant_context <- function(vc) {
  is.list(vc) &&
    length(vc) == 2 &&
    class(vc[[1]]) == "externalptr" &&
    class(vc[[2]]) == "externalptr"
}

#' Check if an object looks like a genotype context
#'
#' @name looks_like_gt_context
#' @title Check if object is a genotype context
#' @param gt the object to check
#' @return TRUE if object is a genotype context
#' @export
looks_like_gt_context <- function(gt) {
  is.list(gt) &&
    length(gt) == 3 &&
    class(gt[[1]]) == "externalptr" &&
    class(gt[[2]]) == "externalptr" &&
    class(gt[[3]]) == "integer"
}


#' Open a VCF or a BCF file
#'
#' fp<-BCFOpen("my.vcf.gz",TRUE)
#' BCFClose(fp)
#' fp<-BCFOpen("my.bcf")
#' BCFClose(fp)
#'
#' @param filename the path to the vcf file
#' @param requireIndex load associated vcf index
#' @return the new VCF reader
#' @title Open a VCF or BCF file
#' @export
#' @examples
#' \dontrun{
#' fp <- BCFOpen("my.vcf.gz", TRUE)
#' BCFClose(fp)
#' fp <- BCFOpen("my.bcf")
#' BCFClose(fp)
#' }
BCFOpen <- function(filename, requireIndex = TRUE) {
  stopifnot(is.character(filename))
  stopifnot(is.logical(requireIndex))
  .Call(RC_RBcfFileOpen, filename, requireIndex, PACKAGE = "RBCFLib")
}

#' Close a VCF reader
#'
#' fp<-BCFOpen("my.vcf.gz",TRUE)
#' BCFClose(fp)
#'
#' @param fp the vcf reader
#' @return true on success
#' @title Close a VCF reader
#' @export
#' @examples
#' \dontrun{
#' fp <- BCFOpen("my.vcf.gz", TRUE)
#' BCFClose(fp)
#' }
BCFClose <- function(fp) {
  stopifnot(looks_like_vcf_context(fp))
  if (!is.null(fp)) {
    .Call(RC_RBcfFileClose, fp, PACKAGE = "RBCFLib")
  }
  invisible(NULL)
}


#' Open a new VCF writer.
#' Must be closed with BCFClose
#'
#' @param fp the vcf reader
#' @param fname the name of the output vcf file
#' @return the writer or NULL on failure
#' @title Open a new VCF writer
#' @export
#' @examples
#' \dontrun{
#' fp <- BCFOpen("my.vcf.gz", TRUE)
#' out <- BCFNewWriter(fp, "output.vcf")
#' BCFClose(out)
#' }
BCFNewWriter <- function(fp, fname = "-") {
  stopifnot(looks_like_vcf_context(fp))
  .Call(RC_RBcfNewWriter, fp, fname, PACKAGE = "RBCFLib")
}

#' Save a Variant in a VCF writer
#'
#' @param fp the vcf reader
#' @param vc the variant
#' @return true on success
#' @title Save a Variant in a VCF writer
#' @export
#' @examples
#' \dontrun{
#' fp <- BCFOpen("my.vcf.gz", TRUE)
#' out <- BCFNewWriter(fp, "output.vcf")
#' while (!is.null(vc <- BCFNext(fp))) {
#'   BCFWriteVariant(out, vc)
#' }
#' BCFClose(fp)
#' BCFClose(out)
#' }
BCFWriteVariant <- function(fp, vc) {
  stopifnot(looks_like_vcf_context(fp))
  .Call(RC_RBcfFileWriteCtx, fp, vc, PACKAGE = "RBCFLib")
}

#' Get number of samples in VCF/BCF file
#'
#' @name BCFNSamples
#' @title Get number of samples
#' @param fp the vcf reader
#' @return the number of samples
#' @export
BCFNSamples <- function(fp) {
  stopifnot(looks_like_vcf_context(fp))
  .Call(RC_RBcfNSamples, fp, PACKAGE = "RBCFLib")
}

#' Get sample names from VCF/BCF file
#'
#' @name BCFSamples
#' @title Get sample names
#' @param fp the vcf reader
#' @return samples
#' @export
BCFSamples <- function(fp) {
  stopifnot(looks_like_vcf_context(fp))
  .Call(RC_RBcfSamples, fp, PACKAGE = "RBCFLib")
}

#' Convert sample names to indices
#'
#' @name BCFSample2Index
#' @title Convert sample names to indices
#' @param fp the vcf reader
#' @param sn the sample name
#' @return 1-based sample indices
#' @export
BCFSample2Index <- function(fp, sn) {
  stopifnot(looks_like_vcf_context(fp))
  sapply(sn, function(S) {
    .Call(RC_BcfConvertSampleToIndex0, fp, S, PACKAGE = "RBCFLib") + 1
  })
}
#' list the indexed chromosomes
#' @param fp the vcf reader
#' @return a list of chromosome
#' @title Get indexed chromosomes
#' @export
BCFChromosomes <- function(fp) {
  stopifnot(looks_like_vcf_context(fp))
  .Call(RC_RBcfSeqNames, fp, PACKAGE = "RBCFLib")
}

#' alias of BCFChromosomes
#' @param fp the vcf reader
#' @return a list of chromosome
#' @title Get indexed contigs
#' @export
BCFContigs <- function(fp) {
  stopifnot(looks_like_vcf_context(fp))
  BCFChromosomes(fp)
}

#' Get sample name by index
#'
#' @name BCFSampleAt
#' @title Get sample name by index
#' @param fp the vcf reader
#' @param idx the 1-based index
#' @return the idx-th sample (1-based)
#' @export
BCFSampleAt <- function(fp, idx) {
  stopifnot(looks_like_vcf_context(fp))
  stopifnot(idx > 0)
  .Call(RC_RBcfSampleAtIndex0, fp, idx - 1, PACKAGE = "RBCFLib")
}

#' Close a VCF reader
#'
#' fp<-BCFOpen("my.vcf.gz",TRUE)
#' BCFClose(fp)
#'
#' @param fp the vcf reader
#' @return true on success
BCFClose <- function(fp) {
  stopifnot(looks_like_vcf_context(fp))
  .Call(RC_RBcfFileClose, fp, PACKAGE = "RBCFLib")
}

#' get dictionary from a VCF reader.
#'
#'
#' fp<-BCFOpen("my.vcf.gz",TRUE)
#' dict<-BCFDictionary(fp)
#' print(dict)
#'      chrom size
#' RF01  RF01 3302
#' RF02  RF02 2687
#' RF03  RF03 2592
#' RF04  RF04 2362
#' RF05  RF05 1579
#' RF06  RF06 1356
#' RF07  RF07 1074
#' RF08  RF08 1059
#' RF09  RF09 1062
#' RF10  RF10  751
#' RF11  RF11  666
#'
#' @param fp the vcf reader
#' @return the vcf dictionary as a table
#' @title Get dictionary from VCF reader
#' @export
BCFDictionary <- function(fp) {
  stopifnot(looks_like_vcf_context(fp))
  .Call(RC_RBcfHeaderDict, fp, PACKAGE = "RBCFLib")
}


#' prepare the VCF reader for a new vcf iteration over a given interval.
#' VCF reader must be associated and opened with a valid index.
#'
#' @param fp the vcf reader
#' @param interval the genomic interval to be scanned
#' @param collect if TRUE return a list of variants in the region or NULL on failure
#' @return TRUE on success or a list of variant if collect=TRUE
#' @title Query VCF by interval
#' @export
BCFQuery <- function(fp, interval, collect = FALSE) {
  stopifnot(looks_like_vcf_context(fp))
  stopifnot(is.character(interval))
  stopifnot(is.logical(collect))
  ret <- .Call(RC_RBcfQueryRegion, fp, interval, PACKAGE = "RBCFLib")

  if (all(collect)) {
    if (!ret) {
      return(NULL)
    }
    # initialize
    variants <- list()
    n <- 1
    # loop while we can read a variant
    while (!is.null(vc <- BCFNext(fp))) {
      #append
      variants[[n]] = vc
      n <- n + 1
    }
    return(variants)
  }
  ret
}

#' read the next Variant Context in the VCF reader
#'
#' @param fp the vcf reader
#' @return an opaque vcf context or NULL
#' @title Read next variant context
#' @examples
#' \dontrun{
#' fp <- BCFOpen("in.vcf.gz")
#' BCFQuery(fp,"RF02:1-1000")
#' while(!is.null(vc<-BCFNext(fp))) {
#'      ## do something with vc
#'      }
#' BCFClose(fp)
#' }
#' @export
BCFNext <- function(fp) {
  stopifnot(looks_like_vcf_context(fp))
  .Call(RC_RBcfNextLine, fp, PACKAGE = "RBCFLib")
}

#' get the numeric index (tid) of the chromosome for this variant
#'
#' @param vc the variant
#' @return the numeric index for this variant
#' @title Get chromosome numeric index (tid)
#' @export
VariantTid <- function(vc) {
  stopifnot(looks_like_variant_context(vc))
  .Call(RC_RBcfCtxRid, vc, PACKAGE = "RBCFLib")
}
#' get chromosome name for this variant
#'
#' @param vc the variant
#' @return the chromosome name
#' @title Get chromosome name for variant
#' @export
VariantContig <- function(vc) {
  stopifnot(looks_like_variant_context(vc))
  .Call(RC_RBcfCtxSeqName, vc, PACKAGE = "RBCFLib")
}
#' alias of VariantContig
#'
#' @param vc the variant
#' @return the chromosome name
#' @title Get chromosome name (alias)
#' @export
VariantChrom <- function(vc) {
  stopifnot(looks_like_variant_context(vc))
  VariantContig(vc)
}
#' get chromosome POS for this variant
#'
#' @param vc the variant
#' @return the starting position for this variant
#' @title Get variant position
#' @export
VariantPos <- function(vc) {
  stopifnot(looks_like_variant_context(vc))
  .Call(RC_RBcfCtxPos, vc, PACKAGE = "RBCFLib")
}
#' alias of VariantPos
#'
#' @param vc the variant
#' @return the starting position for this variant
#' @title Get variant start position (alias)
#' @export
VariantStart <- function(vc) {
  stopifnot(looks_like_variant_context(vc))
  VariantPos(vc)
}

#' return whether variant have got an ID
#'
#' @param vc the variant
#' @return the TRUE if variant have got an ID
#' @title Check if variant has ID
#' @export
VariantHasId <- function(vc) {
  stopifnot(looks_like_variant_context(vc))
  .Call(RC_RBcfCtxHasId, vc, PACKAGE = "RBCFLib")
}

#' return the variant ID
#'
#' @param vc the variant
#' @return the variant ID or NULL
#' @title Get variant ID
#' @export
VariantId <- function(vc) {
  stopifnot(looks_like_variant_context(vc))
  .Call(RC_RBcfCtxId, vc, PACKAGE = "RBCFLib")
}

#' get chromosome END for this variant
#'
#' @param vc the variant
#' @return the END position for this variant
#' @title Get variant end position
#' @export
VariantEnd <- function(vc) {
  stopifnot(looks_like_variant_context(vc))
  .Call(RC_RBcfCtxEnd, vc, PACKAGE = "RBCFLib")
}

#' alias of VariantEnd
#'
#' @param vc the variant
#' @return the END position for this variant
#' @title Get variant stop position (alias)
#' @export
VariantStop <- function(vc) {
  stopifnot(looks_like_variant_context(vc))
  VariantEnd(vc)
}


#' get the number of alleles for this variant
#'
#' @param vc the variant
#' @return the number of alleles for this variant
#' @title Get number of alleles for variant
#' @export
VariantNAlleles <- function(vc) {
  stopifnot(looks_like_variant_context(vc))
  .Call(RC_RBcfCtxNAlleles, vc, PACKAGE = "RBCFLib")
}

#' get the alleles for this variant
#'
#' @param vc the variant
#' @return the alleles for this variant
#' @title Get alleles for variant
#' @export
VariantAlleles <- function(vc) {
  stopifnot(looks_like_variant_context(vc))
  .Call(RC_RBcfCtxAlleles, vc, PACKAGE = "RBCFLib")
}

#' get the reference allele for this variant
#'
#' @param vc the variant
#' @return the reference allele for this variant
#' @title Get reference allele for variant
#' @export
VariantReference <- function(vc) {
  stopifnot(looks_like_variant_context(vc))
  .Call(RC_RBcfCtxReference, vc, PACKAGE = "RBCFLib")
}

#' get the alternate alleles for this variant
#'
#' @param vc the variant
#' @return the alternate alleles for this variant
#' @title Get alternate alleles for variant
#' @export
VariantAltAlleles <- function(vc) {
  stopifnot(looks_like_variant_context(vc))
  .Call(RC_RBcfCtxAlternateAlleles, vc, PACKAGE = "RBCFLib")
}


#' test if QUAL is available for this variant
#'
#' @param vc the variant
#' @return QUAL is available for this variant
#' @title Check if QUAL is available for variant
#' @export
VariantHasQual <- function(vc) {
  stopifnot(looks_like_variant_context(vc))
  .Call(RC_RBcfCtxHasQual, vc, PACKAGE = "RBCFLib")
}

#' test the QUAL for this variant
#'
#' @param vc the variant
#' @return QUAL or NULL
#' @title Get QUAL for variant
#' @export
VariantQual <- function(vc) {
  stopifnot(looks_like_variant_context(vc))
  .Call(RC_RBcfCtxQual, vc, PACKAGE = "RBCFLib")
}


#' test is variant is filtered
#'
#' @param vc the variant
#' @return  variant is filtered
#' @title Check if variant is filtered
#' @export
VariantIsFiltered <- function(vc) {
  stopifnot(looks_like_variant_context(vc))
  .Call(RC_RBcfCtxFiltered, vc, PACKAGE = "RBCFLib")
}

#' return the FILTERS for this variant
#'
#' @param vc the variant
#' @return the filters
#' @title Get filters for variant
#' @export
VariantFilters <- function(vc) {
  stopifnot(looks_like_variant_context(vc))
  .Call(RC_RBcfCtxFilters, vc, PACKAGE = "RBCFLib")
}


#' Check if variant has specific filter
#'
#' @name VariantHasFilter
#' @title Check if variant has specific filter
#' @param vc the variant
#' @param fn filter name
#' @return true if variant is filtered with 'fn'
#' @export
VariantHasFilter <- function(vc, fn) {
  stopifnot(looks_like_variant_context(vc))
  .Call(RC_VariantHasFilter, vc, fn, PACKAGE = "RBCFLib")
}


#' return the types for this variant (as defined in htslib)
#'
#' @param vc the variant
#' @return the types for this variant
#' @title Get types for variant
#' @export
VariantTypes <- function(vc) {
  stopifnot(looks_like_variant_context(vc))
  .Call(RC_RBcfCtxVariantTypes, vc, PACKAGE = "RBCFLib")
}


#' return true if variant is a SNP
#'
#' @param vc the variant
#' @return true if the variant is a SNP
#' @title Check if variant is SNP
#' @export
VariantIsSnp <- function(vc) {
  stopifnot(looks_like_variant_context(vc))
  .Call(RC_RBcfCtxVariantIsSnp, vc, PACKAGE = "RBCFLib")
}

#' max ploidy for this variant
#'
#' @param vc the variant
#' @return max ploidy
#' @title Get max ploidy for variant
#' @export
VariantMaxPloidy <- function(vc) {
  stopifnot(looks_like_variant_context(vc))
  .Call(RC_RBcfCtxVariantMaxPloidy, vc, PACKAGE = "RBCFLib")
}

#' Get number of samples for variant
#'
#' @name VariantNSamples
#' @title Get number of samples for variant
#' @param vc the variant
#' @return the number of samples/genotypes for this variant
#' @export
VariantNSamples <- function(vc) {
  stopifnot(looks_like_variant_context(vc))
  .Call(RC_VariantNSamples, vc, PACKAGE = "RBCFLib")
}

#' Get all genotypes for variant
#'
#' @name VariantGenotypes
#' @title Get all genotypes for variant
#' @param vc the variant
#' @return the number of samples/genotypes for this variant
#' @export
VariantGenotypes <- function(vc) {
  stopifnot(looks_like_variant_context(vc))
  sapply(1:VariantNSamples(vc), function(idx) {
    VariantGenotype(vc, idx)
  })
}


#' return the genotype for a variant
#'
#' @param vc the variant
#' @param nameOrIdx the sample name (slower) or the 1-based sample index
#' @return the given genotype
#' @title Get genotype for variant
#' @export
VariantGenotype <- function(vc, nameOrIdx) {
  stopifnot(looks_like_variant_context(vc))
  if (is.numeric(nameOrIdx)) {
    stopifnot(nameOrIdx > 0)
    nameOrIdx <- as.integer(nameOrIdx) - 1
  } else {
    stopifnot(is.character(nameOrIdx))
  }
  .Call(RC_VariantGetGenotype, vc, nameOrIdx, PACKAGE = "RBCFLib")
}

#' return the 0-based alleles indexes for the given genotype
#'
#' @param gt the genotype
#' @return max ploidy
#' @title Get 0-based alleles indexes for genotype
#' @export
GenotypeAllelesIdx0 <- function(gt) {
  stopifnot(looks_like_gt_context(gt))
  .Call(RC_RBcfCtxVariantGtAllelesIndexes0, gt, PACKAGE = "RBCFLib")
}

#' return the 0-based allele indizes for all genotypes of the variant
#'
#' The returned vector contains the allele-indizes for all genotypes in sample order
#'
#' For example, for a 3 sample VCF on a diploid variant, the returned vector could be
#' ```
#'   GT1  GT2  GT3
#' c(0,1, 0,0, 1,1)
#' ```
#'
#' @param vc the variant
#' @return an integer vector of length `ploidy * n_samples`
#' @title Get 0-based allele indexes for all genotypes
#' @export
VariantGenotypesAlleleIdx0 <- function(vc) {
  stopifnot(looks_like_variant_context(vc))
  .Call(RC_RBcfCtxVariantAllGtAllelesIndexes0, vc, PACKAGE = "RBCFLib")
}


#' Counts the allele-counts per genotype for a given allele
#'
#' The returned vector contains the number of occurance of a specific allele
#' in all genotypes.
#'
#' For example, for a 3 sample VCF on a diploid variant with genotypes
#' ```
#'   GT1  GT2  GT3
#'   0/1  0/0  1/1
#' ```
#' the result will be
#' ```
#' VariantGenotypesAlleleCounts(vc, 0) => c(1,2,0)
#' VariantGenotypesAlleleCounts(vc, 1) => c(1,0,2)
#' ```
#'
#' @param vc the variant
#' @param allele_index The index of the allele to count (0-based)
#' @return an integer vector of length `n_samples`
#' @title Get allele counts per genotype
#' @export
VariantGenotypesAlleleCounts <- function(vc, allele_index = 1) {
  stopifnot(looks_like_variant_context(vc))
  stopifnot(is.numeric(allele_index))
  stopifnot(allele_index >= 0)
  allele_index <- as.integer(allele_index)
  .Call(
    RC_RBcfCtxVariantAllGtAllelesAlleleCounts,
    vc,
    allele_index,
    PACKAGE = "RBCFLib"
  )
}

#' return the genotype strings for all genotypes of the variant
#'
#' The returned vector contains the genotype strings in sample order.
#'
#' For example, for a 3 sample VCF on a diploid variant, the returned vector could be
#' ```
#'    GT1    GT2    GT3
#' c("0/1", "0/0", "1/1")
#' ```
#'
#' @param vc the variant
#' @return a character vector of length `n_samples`
#' @title Get genotype strings for all genotypes
#' @export
VariantGenotypesAlleleStrings <- function(vc) {
  stopifnot(looks_like_variant_context(vc))
  .Call(RC_RBcfCtxVariantAllGtStrings, vc, PACKAGE = "RBCFLib")
}

#' Get ploidy of genotype
#'
#' @name GenotypePloidy
#' @title Get ploidy of genotype
#' @param gt the genotype
#' @return the number of alleles for the genotypes
#' @export
GenotypePloidy <- function(gt) {
  stopifnot(looks_like_gt_context(gt))
  alleles <- GenotypeAllelesIdx0(gt)
  v = 0
  if (!is.null(alleles)) {
    v = length(alleles)
  }
  v
}

#' Check if genotype is homozygous reference
#'
#' @name GenotypeHomRef
#' @title Check if genotype is homozygous reference
#' @param gt the genotype
#' @return TRUE if genotype is diploid and all alleles are reference
#' @export
GenotypeHomRef <- function(gt) {
  stopifnot(looks_like_gt_context(gt))
  alleles <- GenotypeAllelesIdx0(gt)
  !is.null(alleles) &&
    length(alleles) == 2 &&
    alleles[1] == 0 &&
    alleles[2] == 0
}

#' Check if genotype is heterozygous
#'
#' @name GenotypeHet
#' @title Check if genotype is heterozygous
#' @param gt the genotype
#' @return TRUE if genotype is diploid and heterozygous
#' @export
GenotypeHet <- function(gt) {
  stopifnot(looks_like_gt_context(gt))
  alleles <- GenotypeAllelesIdx0(gt)
  !is.null(alleles) &&
    length(alleles) == 2 &&
    alleles[1] != alleles[2] &&
    alleles[1] >= 0 &&
    alleles[2] >= 0
}

#' Check if genotype is homozygous variant
#'
#' @name GenotypeHomVar
#' @title Check if genotype is homozygous variant
#' @param gt the genotype
#' @return TRUE if genotype is diploid and homozygous on alt allele
#' @export
GenotypeHomVar <- function(gt) {
  stopifnot(looks_like_gt_context(gt))
  alleles <- GenotypeAllelesIdx0(gt)
  !is.null(alleles) &&
    length(alleles) == 2 &&
    alleles[1] > 0 &&
    alleles[2] > 0 &&
    alleles[1] == alleles[2]
}

#' Check if genotype is heterozygous non-reference
#'
#' @name GenotypeHetNonRef
#' @title Check if genotype is heterozygous non-reference
#' @param gt the genotype
#' @return TRUE if genotype is diploid and heterozygous and doesn't contains the alt allele
#' @export
GenotypeHetNonRef <- function(gt) {
  stopifnot(looks_like_gt_context(gt))
  alleles <- GenotypeAllelesIdx0(gt)
  !is.null(alleles) &&
    length(alleles) == 2 &&
    alleles[1] != alleles[2] &&
    alleles[1] > 0 &&
    alleles[2] > 0
}

#' Check if genotype is no-call
#'
#' @name GenotypeNoCall
#' @title Check if genotype is no-call
#' @param gt the genotype
#' @return TRUE if genotypes contains no allele '' or any is no call '.'
#' @export
GenotypeNoCall <- function(gt) {
  stopifnot(looks_like_gt_context(gt))
  alleles <- GenotypeAllelesIdx0(gt)
  is.null(alleles) || length(alleles) == 0 || -1 %in% alleles
}

#' Check if genotype is phased
#'
#' @name GenotypePhased
#' @title Check if genotype is phased
#' @param gt the genotype
#' @return TRUE if genotype is phased
#' @export
GenotypePhased <- function(gt) {
  stopifnot(looks_like_gt_context(gt))
  .Call(RC_RBcfCtxVariantGtPhased, gt, PACKAGE = "RBCFLib")
}

#' Get sample name for genotype
#'
#' @name GenotypeSample
#' @title Get sample name for genotype
#' @param gt the genotype
#' @return the name of the sample associated to the genotype
#' @export
GenotypeSample <- function(gt) {
  stopifnot(looks_like_gt_context(gt))
  .Call(RC_GenotypeSample, gt, PACKAGE = "RBCFLib")
}


#' Check if variant has INFO attribute
#'
#' @name VariantHasAttribute
#' @title Check if variant has INFO attribute
#' @param vc the variant
#' @param att the INFO/Attribute
#' @return true if variant has INFO attribute
#' @export
VariantHasAttribute <- function(vc, att) {
  stopifnot(looks_like_variant_context(vc))
  stopifnot(is.character(att))
  .Call(RC_VariantHasAttribute, vc, att, PACKAGE = "RBCFLib")
}


#' Get string attribute from variant INFO
#'
#' @name VariantStringAttribute
#' @title Get string attribute from variant INFO
#' @param vc the variant
#' @param att the INFO/Attribute
#' @param split split strings using commas
#' @return the the INFO attribute for the given key.
#' @export
VariantStringAttribute <- function(vc, att, split = TRUE) {
  stopifnot(looks_like_variant_context(vc))
  stopifnot(is.character(att))
  stopifnot(length(att) == 1)
  s <- .Call(RC_VariantStringAttribute, vc, att, PACKAGE = "RBCFLib")
  if (!is.null(s) && length(s) > 0 && split) {
    s <- unlist(strsplit(s, split = ","))
  }
  s
}

#' Get integer attribute from variant INFO
#'
#' @name VariantIntAttribute
#' @title Get integer attribute from variant INFO
#' @param vc the variant
#' @param att the INFO/Attribute
#' @return the the INFO attribute for the given key.
#' @export
VariantIntAttribute <- function(vc, att) {
  stopifnot(looks_like_variant_context(vc))
  stopifnot(is.character(att))
  stopifnot(length(att) == 1)
  .Call(RC_VariantIntAttribute, vc, att, PACKAGE = "RBCFLib")
}

#' Get float attribute from variant INFO
#'
#' @name VariantFloatAttribute
#' @title Get float attribute from variant INFO
#' @param vc the variant
#' @param att the INFO/Attribute
#' @return the the INFO attribute for the given key.
#' @export
VariantFloatAttribute <- function(vc, att) {
  stopifnot(looks_like_variant_context(vc))
  stopifnot(is.character(att))
  stopifnot(length(att) == 1)
  .Call(RC_VariantFloatAttribute, vc, att, PACKAGE = "RBCFLib")
}

#' Get flag attribute from variant INFO
#'
#' @name VariantFlagAttribute
#' @title Get flag attribute from variant INFO
#' @param vc the variant
#' @param att the INFO/Attribute
#' @return the the INFO attribute for the given key.
#' @export
VariantFlagAttribute <- function(vc, att) {
  stopifnot(looks_like_variant_context(vc))
  stopifnot(is.character(att))
  stopifnot(length(att) == 1)
  .Call(RC_VariantFlagAttribute, vc, att, PACKAGE = "RBCFLib")
}


#' Get filter table from VCF/BCF header
#'
#' @name BCFFilters
#' @title Get filter table from VCF/BCF header
#' @param fp the vcf reader
#' @return a table of filters
#' @export
BCFFilters <- function(fp) {
  .Call(RC_BcfFilterTable, fp, PACKAGE = "RBCFLib")
}

#' Get INFO table from VCF/BCF header
#'
#' @name BCFInfos
#' @title Get INFO table from VCF/BCF header
#' @param fp the vcf reader
#' @return a table of filters
#' @export
BCFInfos <- function(fp) {
  .Call(RC_BcfInfoTable, fp, PACKAGE = "RBCFLib")
}

#' Get FORMAT table from VCF/BCF header
#'
#' @name BCFFormats
#' @title Get FORMAT table from VCF/BCF header
#' @param fp the vcf reader
#' @return a table of filters
#' @export
BCFFormats <- function(fp) {
  .Call(RC_BcfFormatTable, fp, PACKAGE = "RBCFLib")
}

#' Get INFO keys for variant
#'
#' @name VariantInfoIds
#' @title Get INFO keys for variant
#' @param vc the variant
#' @return the list INFOs for this variant
#' @export
VariantInfoIds <- function(vc) {
  stopifnot(looks_like_variant_context(vc))
  .Call(RC_VariantGetInfoKeySet, vc, PACKAGE = "RBCFLib")
}

#' Get FORMAT keys for variant
#'
#' @name VariantFormatIds
#' @title Get FORMAT keys for variant
#' @param vc the variant
#' @return the list FORMATs for this variant
#' @export
VariantFormatIds <- function(vc) {
  stopifnot(looks_like_variant_context(vc))
  .Call(RC_VariantGetFormatKeySet, vc, PACKAGE = "RBCFLib")
}
#' Get VEP annotation table for variant
#'
#' @name VariantVep
#' @title Get VEP annotation table for variant
#' @param vc the variant
#' @return VEP table for this variant or NULL
#' @export
VariantVep <- function(vc) {
  stopifnot(looks_like_variant_context(vc))
  .Call(RC_VariantVepTable, vc, PACKAGE = "RBCFLib")
}
#' Get SnpEff annotation table for variant
#'
#' @name VariantSnpEff
#' @title Get SnpEff annotation table for variant
#' @param vc the variant
#' @return SNPEFF table for this variant or NULL
#' @export
VariantSnpEff <- function(vc) {
  stopifnot(looks_like_variant_context(vc))
  .Call(RC_VariantSnpEffTable, vc, PACKAGE = "RBCFLib")
}


#' Return a specific FORMAT flag value on all genotypes
#'
#' The returned vector contains the flags by sample and attribute number (if the attribute comprises multiple flags).
#'
#' For example, for a 3 sample VCF with a flag having a single logical, the returned vector could be
#' ```
#'    GT1    GT2   GT3
#' c(TRUE, FALSE, TRUE)
#' ```
#' @param vc the variant
#' @param att the attribute name to retrieve
#'
#' @return vector of logicals containing attribute values for all genotypes
#'
#' @seealso
#'    \link{VariantGenotypesIntAttribute},
#'    \link{VariantGenotypesFloatAttribute}
#' @title Get FORMAT flag value on all genotypes
#' @export
VariantGenotypesFlagAttribute <- function(vc, att) {
  stopifnot(looks_like_variant_context(vc))
  stopifnot(is.character(att))
  stopifnot(length(att) == 1)
  as.logical(.Call(
    RC_VariantGenotypesFlagAttribute,
    vc,
    att,
    PACKAGE = "RBCFLib"
  ))
}

#' Return a specific FORMAT integer value on all genotypes
#'
#' The returned vector contains the values by sample and attribute number (if the attribute comprises multiple integer values).
#'
#' For example, for a 3 sample VCF extracting the allelic read depth (AD) on a singleton variant, the results could look like:
#' ```
#'   GT1-REF, GT1-ALT, GT2-REF, GT2-ALT, GT3-REF, GT3-ALT
#' c(     10,     100,      25,      94,      45,       7)
#' ```
#'
#' @param vc the variant
#' @param att the attribute name to retrieve
#'
#' @return vector of numerics containing attribute values for all genotypes
#'
#' @seealso
#'    \link{VariantGenotypesFlagAttribute},
#'    \link{VariantGenotypesFloatAttribute}
#' @title Get FORMAT integer value on all genotypes
#' @export
VariantGenotypesIntAttribute <- function(vc, att) {
  stopifnot(looks_like_variant_context(vc))
  stopifnot(is.character(att))
  stopifnot(length(att) == 1)
  .Call(RC_VariantGenotypesInt32Attribute, vc, att, PACKAGE = "RBCFLib")
}

#' Return a specific FORMAT numeric value on all genotypes
#'
#' The returned vector contains the values by sample and attribute number (if the attribute comprises multiple integer values).
#'
#' For example, for a 3 sample VCF extracting the allelic read fraction (AF) on a multi-allelic variant, the results could look like:
#' ```
#'   GT1-ALT1, GT1-ALT2, GT2-ALT1, GT2-ALT2, GT3-ALT1, GT3-ALT2
#' c(     0.9,      0.1,      0.6,      0.1,      0.1,      0.4)
#' ```
#'
#' @param vc the variant
#' @param att the attribute name to retrieve
#'
#' @return vector of numerics containing attribute values for all genotypes
#'
#' @seealso
#'    \link{VariantGenotypesFlagAttribute},
#'    \link{VariantGenotypesIntAttribute}
#' @title Get FORMAT numeric value on all genotypes
#' @export
VariantGenotypesFloatAttribute <- function(vc, att) {
  stopifnot(looks_like_variant_context(vc))
  stopifnot(is.character(att))
  stopifnot(length(att) == 1)
  .Call(RC_VariantGenotypesFloatAttribute, vc, att, PACKAGE = "RBCFLib")
}

#' Get integer attribute from genotype
#'
#' @name GenotypeIntAttribute
#' @title Get integer attribute from genotype
#' @param gt the genotype
#' @param att the key
#' @return the values for this key+genotype
#' @export
GenotypeIntAttribute <- function(gt, att) {
  .Call(RC_GenotypeInt32Attribute, gt, att, PACKAGE = "RBCFLib")
}


#' Get string attribute from genotype
#'
#' @name GenotypeStringAttribute
#' @title Get string attribute from genotype
#' @param gt the genotype
#' @param att the key
#' @return a String for this flag+genotype or NULL
#' @export
GenotypeStringAttribute <- function(gt, att) {
  .Call(RC_GenotypeStringAttribute, gt, att, PACKAGE = "RBCFLib")
}

#' Check if genotype is filtered
#'
#' @name GenotypeFiltered
#' @title Check if genotype is filtered
#' @param gt the genotype
#' @return TRUE if genotype is filtered (FORMAT/FT is set)
#' @export
GenotypeFiltered <- function(gt) {
  flt <- GenotypeStringAttribute(gt, "FT")
  !(is.null(flt) || length(flt) == 0 || flt == ".")
}

#' Get read depth for genotype
#'
#' @name GenotypeDp
#' @title Get read depth for genotype
#' @param gt the genotype
#' @return the DP or -1
#' @export
GenotypeDp <- function(gt) {
  v <- GenotypeIntAttribute(gt, "DP")
  if (length(v) != 1) {
    return(-1)
  }
  v[1]
}
#' Check if genotype has read depth
#'
#' @name GenotypeHasDp
#' @title Check if genotype has read depth
#' @param gt the genotype
#' @return TRUE if DP is available
#' @export
GenotypeHasDp <- function(gt) {
  GenotypeDp(gt) >= 0
}
#' Get genotype quality for genotype
#'
#' @name GenotypeGq
#' @title Get genotype quality for genotype
#' @param gt the genotype
#' @return the GQ
#' @export
GenotypeGq <- function(gt) {
  v <- GenotypeIntAttribute(gt, "GQ")
  if (length(v) != 1) {
    return(-1)
  }
  v[1]
}
#' Check if genotype has quality score
#'
#' @name GenotypeHasGq
#' @title Check if genotype has quality score
#' @param gt the genotype
#' @return TRUE if GQ is available
#' @export
GenotypeHasGq <- function(gt) {
  GenotypeGq(gt) >= 0
}
#' Get PL (Phred-scaled likelihoods) for genotype
#'
#' @name GenotypePl
#' @title Get PL for genotype
#' @param gt the genotype
#' @return PL or empty vector
GenotypePl <- function(gt) {
  GenotypeIntAttribute(gt, "PL")
}
#' Check if PL is available for genotype
#'
#' @name GenotypeHasPl
#' @title Check if PL is available for genotype
#' @param gt the genotype
#' @return TRUE if PL is available
GenotypeHasPl <- function(gt) {
  length(GenotypePl(gt)) > 0
}

#' Get AD (allelic depths) for genotype
#'
#' @name GenotypeAd
#' @title Get AD for genotype
#' @param gt the genotype
#' @return AD or empty vector
GenotypeAd <- function(gt) {
  GenotypeIntAttribute(gt, "AD")
}

#' Check if AD is available for genotype
#'
#' @name GenotypeHasAd
#' @title Check if AD is available for genotype
#' @param gt the genotype
#' @return TRUE if AD is available
GenotypeHasAd <- function(gt) {
  length(GenotypeAd(gt)) > 0
}
