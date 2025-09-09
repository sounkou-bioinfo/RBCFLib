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
#' @param vbi_vcf_ctx VBI VCF context object from VCFLoad()
#' @param region Region string (e.g., "chr1:1000-2000")
#' @return data.frame with variant information
#' @export
VBIQueryRegionCGRanges <- function(vbi_vcf_ctx, region) {
    .Call(
        RC_VBI_query_region_cgranges,
        vbi_vcf_ctx,
        as.character(region),
        as.logical(FALSE), # include_info
        as.logical(FALSE), # include_format
        as.logical(FALSE), # include_genotypes
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

#' Query VBI index by region (legacy function)
#'
#' @param vbi_vcf_ctx VBI VCF context object from VCFLoad()
#' @param region Region string (e.g., "chr1:1000-2000") for query
#' @return data.frame with variant information
#' @export
VBIQueryRange <- function(vbi_vcf_ctx, region) {
    # Use the standardized query function instead of the legacy one
    VBIQueryRegion(vbi_vcf_ctx, region, include_info = FALSE)
}

#' Query VBI index by contiguous variant index range
#'
#' @param vbi_vcf_ctx VBI VCF context object from VCFLoad()
#' @param start_index Start index for index-based query (1-based)
#' @param end_index End index for index-based query (1-based)
#' @return data.frame with variant information
#' @export
VBIQueryByIndices <- function(vbi_vcf_ctx, start_index, end_index) {
    .Call(
        RC_VBI_query_by_indices,
        vbi_vcf_ctx,
        as.integer(start_index),
        as.integer(end_index),
        as.logical(FALSE), # include_info
        as.logical(FALSE), # include_format
        as.logical(FALSE), # include_genotypes
        PACKAGE = "RBCFLib"
    )
}
#' Load a VCF file with VBI index integration
#'
#' Creates a unified VCF object that includes the VCF file, header information,
#' and VBI index for fast queries. This follows the RBCF pattern but with VBI integration.
#'
#' @param vcf_path Path to the VCF/BCF file
#' @param vbi_path Path to the VBI index file (optional, will auto-detect if NULL)
#' @return External pointer to VBI VCF context object
#' @export
#' @examples
#' \dontrun{
#' vcf_file <- system.file("exdata", "imputed.gt.vcf.gz", package = "RBCFLib")
#' vcf_obj <- VCFLoad(vcf_file)
#'
#' # Query with basic fields
#' hits <- VBIQueryRegion(vcf_obj, "chr21:5030082-5030356")
#'
#' # Query with INFO fields
#' hits_with_info <- VBIQueryRegion(vcf_obj, "chr21:5030082-5030356",
#'                                  include_info = TRUE)
#' }
VCFLoad <- function(vcf_path, vbi_path = NULL) {
    .Call(
        RC_VBI_vcf_load,
        as.character(vcf_path),
        if (is.null(vbi_path)) NULL else as.character(vbi_path),
        PACKAGE = "RBCFLib"
    )
}

#' VBI query by region with full INFO parsing
#'
#' Query VCF variants by genomic region with comprehensive field extraction.
#' This function provides access to all VCF fields including INFO data,
#' unlike the basic VBI query functions.
#'
#' @param vbi_vcf_ctx VBI VCF context object from VCFLoad()
#' @param region Region string (e.g., "chr1:1000-2000")
#' @param include_info Logical, whether to include INFO fields (default: FALSE)
#' @return data.frame with comprehensive variant information including all INFO fields when requested
#' @export
#' @examples
#' \dontrun{
#' vcf_file <- system.file("exdata", "imputed.gt.vcf.gz", package = "RBCFLib")
#' vcf_obj <- VCFLoad(vcf_file)
#'
#' # Basic query
#' hits <- VBIQueryRegion(vcf_obj, "chr21:5030082-5030356")
#'
#' # Query with INFO fields
#' hits_info <- VBIQueryRegion(vcf_obj, "chr21:5030082-5030356",
#'                            include_info = TRUE)
#' }

#' VBI query with genotype data access
#'
#' Query VCF variants by genomic region with genotype data access.
#' Returns variant context that can be used with genotype accessor functions.
#'
#' @param vbi_vcf_ctx VBI VCF context object from VCFLoad()
#' @param region Region string (e.g., "chr1:1000-2000")
#' @param include_info Logical, whether to include INFO fields (default: FALSE)
#' @param include_format Logical, whether to include FORMAT fields (default: FALSE)
#' @param include_genotypes Logical, whether to include genotype data (default: FALSE)
#' @return data.frame with comprehensive variant information
#' @export
#' @examples
#' \dontrun{
#' vcf_file <- system.file("exdata", "imputed.gt.vcf.gz", package = "RBCFLib")
#' vcf_obj <- VCFLoad(vcf_file)
#'
#' # Basic query
#' hits <- VBIQueryRegion(vcf_obj, "chr21:5030082-5030356")
#'
#' # Query with INFO and genotype data
#' hits_full <- VBIQueryRegion(vcf_obj, "chr21:5030082-5030356",
#'                            include_info = TRUE, include_genotypes = TRUE)
#' }
VBIQueryRegion <- function(
    vbi_vcf_ctx,
    region,
    include_info = FALSE,
    include_format = FALSE,
    include_genotypes = FALSE
) {
    .Call(
        RC_VBI_query_region,
        vbi_vcf_ctx,
        as.character(region),
        as.logical(include_info),
        as.logical(include_format),
        as.logical(include_genotypes),
        PACKAGE = "RBCFLib"
    )
}

#' Get VBI VCF context samples
#'
#' Retrieve sample names from a VBI VCF context object.
#'
#' @param vbi_vcf_ctx VBI VCF context object from VCFLoad()
#' @return Character vector of sample names
#' @export
VBISamples <- function(vbi_vcf_ctx) {
    .Call(RC_VBI_samples, vbi_vcf_ctx, PACKAGE = "RBCFLib")
}

#' Get number of samples in VBI VCF context
#'
#' @param vbi_vcf_ctx VBI VCF context object from VCFLoad()
#' @return Integer, number of samples
#' @export
VBINSamples <- function(vbi_vcf_ctx) {
    .Call(RC_VBI_nsamples, vbi_vcf_ctx, PACKAGE = "RBCFLib")
}

#' Get sample at specific index from VBI VCF context
#'
#' @param vbi_vcf_ctx VBI VCF context object from VCFLoad()
#' @param index 1-based sample index
#' @return Character, sample name
#' @export
VBISampleAt <- function(vbi_vcf_ctx, index) {
    .Call(RC_VBI_sample_at, vbi_vcf_ctx, as.integer(index), PACKAGE = "RBCFLib")
}

#' Convert sample name to index in VBI VCF context
#'
#' @param vbi_vcf_ctx VBI VCF context object from VCFLoad()
#' @param sample_name Character, sample name
#' @return Integer, 1-based sample index (0 if not found)
#' @export
VBISample2Index <- function(vbi_vcf_ctx, sample_name) {
    .Call(
        RC_VBI_sample2index,
        vbi_vcf_ctx,
        as.character(sample_name),
        PACKAGE = "RBCFLib"
    )
}

#' Get INFO fields table from VBI VCF context
#'
#' @param vbi_vcf_ctx VBI VCF context object from VCFLoad()
#' @return data.frame with INFO field definitions
#' @export
VBIInfos <- function(vbi_vcf_ctx) {
    .Call(RC_VBI_infos, vbi_vcf_ctx, PACKAGE = "RBCFLib")
}

#' Get FORMAT fields table from VBI VCF context
#'
#' @param vbi_vcf_ctx VBI VCF context object from VCFLoad()
#' @return data.frame with FORMAT field definitions
#' @export
VBIFormats <- function(vbi_vcf_ctx) {
    .Call(RC_VBI_formats, vbi_vcf_ctx, PACKAGE = "RBCFLib")
}

#' Get FILTER fields table from VBI VCF context
#'
#' @param vbi_vcf_ctx VBI VCF context object from VCFLoad()
#' @return data.frame with FILTER field definitions
#' @export
VBIFilters <- function(vbi_vcf_ctx) {
    .Call(RC_VBI_filters, vbi_vcf_ctx, PACKAGE = "RBCFLib")
}
