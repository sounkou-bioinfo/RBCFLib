#' Index a FASTA file using htslib faidx
#'
#' @param fasta_path Path to the FASTA file
#' @return Path to the generated .fai index file
#' @export
FaidxIndexFasta <- function(fasta_path) {
  if (!file.exists(fasta_path)) {
    stop("FASTA file does not exist: ", fasta_path)
  }
  .Call(RC_FaidxIndexFasta, fasta_path)
}

#' Fetch a sequence region from a FASTA file using htslib faidx
#'
#' @param fasta_path Path to the FASTA file
#' @param seqname Chromosome/contig name (e.g., "chr1")
#' @param start 1-based start position (inclusive)
#' @param end 1-based end position (inclusive)
#' @return Character vector of the fetched sequence
#' @export
FaidxFetchRegion <- function(fasta_path, seqname, start, end) {
  if (!file.exists(fasta_path)) {
    stop("FASTA file does not exist: ", fasta_path)
  }
  if (!file.exists(paste0(fasta_path, ".fai"))) {
    FaidxIndexFasta(fasta_path)
  }
  .Call(
    RC_FaidxFetchRegion,
    fasta_path,
    seqname,
    as.integer(start),
    as.integer(end)
  )
}

# should probably add zstd

#' Download human reference genomes (GRCh37 and GRCh38)
#'
#' Downloads the GRCh37 and GRCh38 human reference FASTA files and their indexes to the specified directories.
#' @importFrom utils download.file
#' @param grch37_dir Directory to store GRCh37 reference (default: ~/GRCh37)
#' @param grch38_dir Directory to store GRCh38 reference (default: ~/GRCh38)
#' @param urls Named list of URLs for genome resources (with defaults)
#' @param cytoband Logical; also download cytoband files (default: FALSE)
#' @param chain Logical; also download chain files for liftover (default: FALSE)
#' @param method Method to be used for downloading files. Possible values are:
#'   \itemize{
#'     \item \code{"wget"} (default) - Uses the external wget utility which provides better functionality
#'           for large files including automatic retries and the ability to resume interrupted downloads
#'     \item \code{"curl"} - Uses the external curl utility
#'     \item \code{"auto"} - Let R choose the best method available
#'     \item \code{"internal"} - Uses R's internal download method
#'     \item \code{"libcurl"} - Uses the libcurl library
#'     \item \code{"lynx"} - Uses the lynx browser in text mode
#'   }
#' @param extra Character vector of additional command-line arguments for the 'wget' and 'curl' methods.
#'   For example, to use wget's continue feature, you can provide \code{extra = c("-C", "on")} to resume
#'   partially downloaded files. See \code{?download.file} for more details.
#'
#'   The \code{wget} method is recommended for large genome files as it has excellent support for
#'   resuming interrupted downloads, which is important when downloading multi-gigabyte files.
#'   See \code{?download.file} for more details on these methods.
#' @return Named list with paths to downloaded files
#' @export
DownloadHumanReferenceGenomes <- function(
  grch37_dir = file.path(Sys.getenv("HOME"), "GRCh37"),
  grch38_dir = file.path(Sys.getenv("HOME"), "GRCh38"),
  urls = list(
    grch37_fasta = "http://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz",
    grch38_fasta = paste0(
      "http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/",
      "GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/",
      "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
    ),
    grch37_cytoband = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz",
    grch38_cytoband = "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz",
    grch37_chain = "http://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/hg18ToHg19.over.chain.gz",
    grch38_chain_18 = "http://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/hg18ToHg38.over.chain.gz",
    grch38_chain_19 = "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz"
  ),
  cytoband = FALSE,
  chain = FALSE,
  method = "wget",
  extra = getOption("download.file.extra")
) {
  dir.create(grch37_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(grch38_dir, showWarnings = FALSE, recursive = TRUE)

  # GRCh37
  grch37_fasta_gz <- file.path(grch37_dir, basename(urls$grch37_fasta))

  grch37_fasta <- gsub(".gz$", "", grch37_fasta_gz)

  grch37_fai <- paste0(grch37_fasta, ".fai")

  if (!file.exists(grch37_fasta)) {
    tryCatch(
      {
        download.file(
          urls$grch37_fasta,
          grch37_fasta_gz,
          mode = "wb",
          method = method,
          extra = extra
        )
        DecompressFile(grch37_fasta_gz, grch37_fasta)
      },
      error = function(e) {
        warning(paste(
          "Error downloading/decompressing GRCh37 reference:",
          e$message,
          "\nContinuing with execution..."
        ))
      }
    )
  }
  if (file.exists(grch37_fasta) && !file.exists(grch37_fai)) {
    tryCatch(
      {
        FaidxIndexFasta(grch37_fasta)
      },
      error = function(e) {
        warning(paste(
          "Error indexing GRCh37 reference:",
          e$message,
          "\nContinuing with execution..."
        ))
      }
    )
  }

  # GRCh38
  grch38_fasta_gz <- file.path(
    grch38_dir,
    basename(urls$grch38_fasta)
  )
  grch38_fasta <- gsub(".gz$", "", grch38_fasta_gz)
  grch38_fai <- paste0(grch38_fasta, ".fai")

  if (!file.exists(grch38_fasta)) {
    tryCatch(
      {
        download.file(
          urls$grch38_fasta,
          grch38_fasta_gz,
          mode = "wb",
          method = method,
          extra = extra
        )
        DecompressFile(grch38_fasta_gz, grch38_fasta)
      },
      error = function(e) {
        warning(paste(
          "Error downloading/decompressing GRCh38 reference:",
          e$message,
          "\nContinuing with execution..."
        ))
      }
    )
  }
  if (file.exists(grch38_fasta) && !file.exists(grch38_fai)) {
    tryCatch(
      {
        FaidxIndexFasta(grch38_fasta)
      },
      error = function(e) {
        warning(paste(
          "Error indexing GRCh38 reference:",
          e$message,
          "\nContinuing with execution..."
        ))
      }
    )
  }

  # Optionally download cytoband and chain files
  cytoband_files <- NULL
  chain_files <- NULL
  if (cytoband) {
    cytoband_files <- list(
      grch37 = file.path(grch37_dir, basename(urls$grch37_cytoband)),
      grch38 = file.path(grch38_dir, basename(urls$grch38_cytoband))
    )
    if (!file.exists(cytoband_files$grch37)) {
      tryCatch(
        {
          download.file(
            urls$grch37_cytoband,
            cytoband_files$grch37,
            mode = "wb",
            method = method,
            extra = extra
          )
        },
        error = function(e) {
          warning(paste(
            "Error downloading GRCh37 cytoband file:",
            e$message,
            "\nContinuing with execution..."
          ))
        }
      )
    }
    if (!file.exists(cytoband_files$grch38)) {
      tryCatch(
        {
          download.file(
            urls$grch38_cytoband,
            cytoband_files$grch38,
            mode = "wb",
            method = method,
            extra = extra
          )
        },
        error = function(e) {
          warning(paste(
            "Error downloading GRCh38 cytoband file:",
            e$message,
            "\nContinuing with execution..."
          ))
        }
      )
    }
  }
  if (chain) {
    chain_files <- list(
      grch37 = file.path(grch37_dir, basename(urls$grch37_chain)),
      grch38_18 = file.path(grch38_dir, basename(urls$grch38_chain_18)),
      grch38_19 = file.path(grch38_dir, basename(urls$grch38_chain_19))
    )
    if (!file.exists(chain_files$grch37)) {
      tryCatch(
        {
          download.file(
            urls$grch37_chain,
            chain_files$grch37,
            mode = "wb",
            method = method,
            extra = extra
          )
        },
        error = function(e) {
          warning(paste(
            "Error downloading GRCh37 chain file:",
            e$message,
            "\nContinuing with execution..."
          ))
        }
      )
    }
    if (!file.exists(chain_files$grch38_18)) {
      tryCatch(
        {
          download.file(
            urls$grch38_chain_18,
            chain_files$grch38_18,
            mode = "wb",
            method = method,
            extra = extra
          )
        },
        error = function(e) {
          warning(paste(
            "Error downloading GRCh38 chain file (hg18):",
            e$message,
            "\nContinuing with execution..."
          ))
        }
      )
    }
    if (!file.exists(chain_files$grch38_19)) {
      tryCatch(
        {
          download.file(
            urls$grch38_chain_19,
            chain_files$grch38_19,
            mode = "wb",
            method = method,
            extra = extra
          )
        },
        error = function(e) {
          warning(paste(
            "Error downloading GRCh38 chain file (hg19):",
            e$message,
            "\nContinuing with execution..."
          ))
        }
      )
    }
  }

  list(
    grch37_fasta = grch37_fasta,
    grch37_fai = grch37_fai,
    grch38_fasta = grch38_fasta,
    grch38_fai = grch38_fai,
    cytoband = cytoband_files,
    chain = chain_files
  )
}

#' Decompress a compressed file
#'
#' Decompresses various compressed file formats (.gz, .bz2, .xz) using base R functionality, borrowing
#' from the R.utils package for the underlying logic. This function is designed to work with
#' without relying on external dependencies.
#' ref : https://github.com/HenrikBengtsson/R.utils/blob/74def095eaa244e355d05fdf790ee6393dad1d99/R/compressFile.R#L15
#' @param input_file Path to the compressed file to decompress
#' @param output_file Path for the decompressed output file. If NULL, the input filename without extension is used
#' @param remove_input Logical; whether to remove the input file after successful decompression (default: FALSE)
#' @param block_size Size of data chunks to process at once (in bytes, default: 1e8)
#' @return Path to the decompressed file
#' @export
DecompressFile <- function(
  input_file,
  output_file = NULL,
  remove_input = FALSE,
  block_size = 1e8
) {
  if (!file.exists(input_file)) {
    stop("Input file does not exist: ", input_file)
  }

  # Determine the file type based on extension
  ext <- tolower(tools::file_ext(input_file))

  # Create output filename if not provided
  if (is.null(output_file)) {
    output_file <- sub(paste0("\\.", ext, "$"), "", input_file)
  }

  # Choose the appropriate connection function based on file extension
  if (ext == "gz") {
    conn_fun <- gzfile
  } else if (ext == "bz2") {
    conn_fun <- bzfile
  } else if (ext == "xz") {
    conn_fun <- xzfile
  } else {
    stop(
      "Unsupported file extension: ",
      ext,
      ". Supported extensions are .gz, .bz2, and .xz"
    )
  }

  # Use tryCatch to handle potential connection errors
  tryCatch(
    {
      # Open connections and decompress
      input_conn <- conn_fun(input_file, "rb")
      on.exit(
        if (!is.null(input_conn) && isOpen(input_conn)) close(input_conn),
        add = TRUE
      )

      output_conn <- file(output_file, "wb")
      on.exit(
        if (!is.null(output_conn) && isOpen(output_conn)) close(output_conn),
        add = TRUE
      )

      # Read and write in chunks
      data <- readBin(input_conn, "raw", n = block_size)
      while (length(data) > 0) {
        writeBin(data, output_conn)
        data <- readBin(input_conn, "raw", n = block_size)
      }

      # Close connections explicitly (on.exit will handle this as a fallback)
      close(input_conn)
      close(output_conn)
    },
    error = function(e) {
      # Log the error but continue without failing
      warning(paste(
        "Connection error during decompression:",
        e$message,
        "\nContinuing with execution..."
      ))
    }
  )

  # Remove input file if requested
  if (remove_input && file.exists(input_file)) {
    file.remove(input_file)
  }

  # Return the output file path even if there was an error
  return(output_file)
}


# Function to collect output from file (handles both text and binary)
collect_output <- function(file_path, ...) {
  if (!file.exists(file_path) || file.info(file_path)$size == 0) {
    return(character(0))
  }

  result <- tryCatch(
    {
      con <- file(file_path, "r")
      on.exit(close(con), add = TRUE)
      readLines(con)
    },
    error = function(e) {
      # Handle binary output by reading as raw if text reading fails
      con <- file(file_path, "rb")
      on.exit(close(con), add = TRUE)
      readBin(con, what = "raw", n = file.info(file_path)$size)
    }
  )

  return(result)
}

#' Get the path to the bcftools binary
#' @export
BCFTOOLS_PLUGINS <- function() {
  return(system.file("bin", "plugins", package = "RBCFLib"))
}
