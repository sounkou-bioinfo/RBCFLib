#' Extract variant ranges from a VBI index pointer
#' @param VbiPtr External pointer to VBI index
#' @param n Number of ranges to extract (default: all)
#' @return data.frame with chrom, start, end, label
#' @export
VBIExtractRanges <- function(VbiPtr, n = NA) {
    .Call(RC_VBI_extract_ranges, VbiPtr, n)
}

#' cgranges R binding: create, add, index, overlap, destroy
#' @export
CGRangesCreate <- function() .Call(RC_cgranges_create)
#' @export
CGRangesAdd <- function(cr, chrom, start, end, label) {
    .Call(
        RC_cgranges_add,
        cr,
        as.character(chrom),
        as.integer(start),
        as.integer(end),
        as.integer(label)
    )
}
#' @export
CGRangesIndex <- function(cr) .Call(RC_cgranges_index, cr)
#' @export
CGRangesOverlap <- function(cr, chrom, start, end) {
    .Call(
        RC_cgranges_overlap,
        cr,
        as.character(chrom),
        as.integer(start),
        as.integer(end)
    )
}
#' @export
CGRangesDestroy <- function(cr) .Call(RC_cgranges_destroy, cr)
#' Load a VBI index file and return an external pointer
#'
#' @param vbi_path Path to the VBI index file
#' @return External pointer to the loaded VBI index
#' @export
VBIIndexLoad <- function(vbi_path) {
    .Call(RC_VBI_load_index, as.character(vbi_path), PACKAGE = "RBCFLib")
}
#' Query VBI index by region using cgranges (fast interval tree)
#'
#' @param VcfPath Path to the VCF/BCF file
#' @param vbi_ptr External pointer to loaded VBI index
#' @param region Region string (e.g., "chr1:1000-2000")
#' @return Character vector of VCF records for the region
#' @export
VBIQueryRegionCGRanges <- function(VcfPath, vbi_ptr, region) {
    .Call(
        RC_VBI_query_region_cgranges,
        as.character(VcfPath),
        vbi_ptr,
        as.character(region),
        PACKAGE = "RBCFLib"
    )
}
#' Print the first n lines of a VBI index (from an external pointer)
#'
#' @param vbi_ptr External pointer to loaded VBI index
#' @param n Number of lines to print (default: 10)
#' @export
VBIPrintIndex <- function(vbi_ptr, n = 10) {
    invisible(.Call(
        RC_VBI_print_index,
        vbi_ptr,
        as.integer(n),
        PACKAGE = "RBCFLib"
    ))
}
#' VBI Index
#'
#' Create a VBI index for a VCF/BCF file.

#' @param VcfPath Path to the VCF/BCF file
#' @param VbiPath Path to the VBI index file
#' @param Threads Number of threads to use (default: 1)
#' @return Path to index file
#' @export
VBIIndex <- function(VcfPath, VbiPath, Threads = 1) {
    .Call(
        RC_VBI_index,
        as.character(VcfPath),
        as.character(VbiPath),
        as.integer(Threads)
    )
}

#' Query VBI index by region

#' @param VcfPath Path to the VCF/BCF file
#' @param VbiPath Path to the VBI index file
#' @param Region Region string (e.g., "chr1:1000-2000") for query
#' @param Threads Number of threads to use (default: 1)
#' @return List of results
#' @export
VBIQueryRange <- function(VcfPath, VbiPath, Region, Threads = 1) {
    .Call(
        RC_VBI_query_range,
        as.character(VcfPath),
        VbiPath,
        as.character(Region),
        as.integer(Threads),
        PACKAGE = "RBCFLib"
    )
}

#' Query VBI index by marker index range

#' @param VcfPath Path to the VCF/BCF file
#' @param VbiPath Path to the VBI index file
#' @param StartIdx Start index for index-based query

#' @param EndingIndice End index for index-based query
#' @param Threads Number of threads to use (default: 1)
#' @return List of results
#' @export
VBIQueryByIndices <- function(
    VcfPath,
    VbiPath,
    StartingIndice,
    EndingIndice,
    Threads = 1
) {
    .Call(
        RC_VBI_query_by_indices,
        as.character(VcfPath),
        VbiPath,
        as.integer(StartingIndice),
        as.integer(EndingIndice),
        as.integer(Threads),
        PACKAGE = "RBCFLib"
    )
}
