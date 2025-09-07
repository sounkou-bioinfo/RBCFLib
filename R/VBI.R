#' Get C-level memory usage of a VBI index pointer (in bytes)
#' @param vbi_ptr External pointer to VBI index
#' @return Numeric, memory usage in bytes
#' @export
VBIIndexMemoryUsage <- function(vbi_ptr) {
    .Call(RC_VBI_index_memory_usage, vbi_ptr, PACKAGE = "RBCFLib")
}
#' Vectorized overlap query for CGRanges
#'
#' Performs overlap queries against a CGRanges index for multiple intervals simultaneously.
#' This is more efficient than calling CGRangesOverlap multiple times.
#'
#' @param cr External pointer to CGRanges object
#' @param chrom Character vector of chromosome names
#' @param start Integer vector of start positions (1-based)
#' @param end Integer vector of end positions (1-based, inclusive)
#' @return List of integer vectors, each containing the indices of overlapping intervals
#' @export
CGRangesOverlapVec <- function(cr, chrom, start, end) {
    .Call(
        RC_cgranges_overlap,
        cr,
        as.character(chrom),
        as.integer(start),
        as.integer(end),
        PACKAGE = "RBCFLib"
    )
}

#' Extract intervals by index from CGRanges
#'
#' Retrieves specific intervals from a CGRanges object using their 1-based indices.
#' Useful for extracting intervals found through overlap queries.
#'
#' @param cr External pointer to CGRanges object
#' @param indices Integer vector of 1-based indices of intervals to extract
#' @return Data frame with columns: chrom, start, end, label
#' @export
CGRangesExtractByIndex <- function(cr, indices) {
    .Call(RC_cgranges_extract_by_index, cr, as.integer(indices))
}
#' Extract variant ranges from a VBI index pointer
#' @param VbiPtr External pointer to VBI index
#' @param n Number of ranges to extract (default: all)
#' @return data.frame with chrom, start, end, label
#' @export
VBIExtractRanges <- function(VbiPtr, n = NA) {
    storage.mode(n) <- "integer"
    .Call(RC_VBI_extract_ranges, VbiPtr, n, PACKAGE = "RBCFLib")
}

#' Create a new CGRanges object
#'
#' Creates a new empty CGRanges object for storing genomic intervals.
#' CGRanges provides fast interval overlap queries using an interval tree data structure.
#'
#' @return External pointer to a new CGRanges object
#' @seealso \code{\link{CGRangesAdd}}, \code{\link{CGRangesIndex}}, \code{\link{CGRangesOverlap}}
#' @export
#' @examples
#' \dontrun{
#' # Create a new CGRanges object
#' cr <- CGRangesCreate()
#'
#' # Add some intervals
#' CGRangesAdd(cr, "chr1", 100, 200, 1)
#' CGRangesAdd(cr, "chr1", 150, 250, 2)
#'
#' # Index for fast queries
#' CGRangesIndex(cr)
#'
#' # Query overlaps
#' overlaps <- CGRangesOverlap(cr, "chr1", 180, 220)
#' }
CGRangesCreate <- function() .Call(RC_cgranges_create)
#' Add an interval to CGRanges
#'
#' Adds a genomic interval to an existing CGRanges object. After adding all intervals,
#' call \code{\link{CGRangesIndex}} to build the interval tree for fast queries.
#'
#' @param cr External pointer to CGRanges object
#' @param chrom Character, chromosome name
#' @param start Integer, start position (1-based)
#' @param end Integer, end position (1-based, inclusive)
#' @param label Integer, user-defined label for this interval
#' @return Invisible NULL
#' @seealso \code{\link{CGRangesCreate}}, \code{\link{CGRangesIndex}}
#' @export
CGRangesAdd <- function(cr, chrom, start, end, label) {
    .Call(
        RC_cgranges_add,
        cr,
        as.character(chrom),
        as.integer(start),
        as.integer(end),
        as.integer(label),
        PACKAGE = "RBCFLib"
    )
}
#' Index CGRanges for fast queries
#'
#' Builds an interval tree index for the CGRanges object to enable fast overlap queries.
#' This must be called after adding all intervals and before performing any queries.
#'
#' @param cr External pointer to CGRanges object
#' @return Invisible NULL
#' @seealso \code{\link{CGRangesAdd}}, \code{\link{CGRangesOverlap}}
#' @export
CGRangesIndex <- function(cr) .Call(RC_cgranges_index, cr)
#' Query CGRanges for overlapping intervals
#'
#' Finds all intervals in the CGRanges object that overlap with the specified query region.
#' The CGRanges object must be indexed first using \code{\link{CGRangesIndex}}.
#'
#' @param cr External pointer to indexed CGRanges object
#' @param chrom Character, chromosome name to query
#' @param start Integer, start position of query region (1-based)
#' @param end Integer, end position of query region (1-based, inclusive)
#' @return Integer vector of indices (1-based) of overlapping intervals
#' @seealso \code{\link{CGRangesIndex}}, \code{\link{CGRangesOverlapVec}}
#' @export
CGRangesOverlap <- function(cr, chrom, start, end) {
    .Call(
        RC_cgranges_overlap,
        cr,
        as.character(chrom),
        as.integer(start),
        as.integer(end),
        PACKAGE = "RBCFLib"
    )
}
#' Destroy a CGRanges object
#'
#' Frees the memory associated with a CGRanges object. After calling this function,
#' the CGRanges object should not be used anymore.
#'
#' @param cr External pointer to CGRanges object
#' @return Invisible NULL
#' @note This function is automatically called when the R object is garbage collected,
#'   but can be called explicitly to free memory immediately.
#' @export
CGRangesDestroy <- function(cr) {
    .Call(RC_cgranges_destroy, cr, PACKAGE = "RBCFLib")
}
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
        as.integer(Threads),
        PACKAGE = "RBCFLib"
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
