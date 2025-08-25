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
#' @param vbi_ptr External pointer to loaded VBI index
#' @param region Region string (e.g., "chr1:1000-2000")
#' @return Integer vector of 1-based indices
#' @export
VBIQueryRegionCGRanges <- function(vbi_ptr, region) {
    .Call(
        RC_VBI_query_region_cgranges,
        vbi_ptr,
        as.character(region),
        PACKAGE = "RBCFLib"
    )
}
#' Print the first n lines of a VBI index file
#'
#' @param vbi_path Path to the VBI index file
#' @param n Number of lines to print (default: 10)
#' @export
VBIPrintIndex <- function(vbi_path, n = 10) {
    invisible(.Call(
        RC_VBI_print_index,
        as.character(vbi_path),
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
        as.character(VbiPath),
        as.character(Region),
        as.integer(Threads),
        PACKAGE = "RBCFLib"
    )
}

#' Query VBI index by marker index range

#' @param VcfPath Path to the VCF/BCF file
#' @param VbiPath Path to the VBI index file
#' @param StartIdx Start index for index-based query
#' @param EndIdx End index for index-based query
#' @param Threads Number of threads to use (default: 1)
#' @return List of results
#' @export
VBIQueryIndex <- function(VcfPath, VbiPath, StartIdx, EndIdx, Threads = 1) {
    .Call(
        RC_VBI_query_index,
        as.character(VcfPath),
        as.character(VbiPath),
        as.integer(StartIdx),
        as.integer(EndIdx),
        as.integer(Threads),
        PACKAGE = "RBCFLib"
    )
}
