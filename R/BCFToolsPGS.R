#' Compute Polygenic Score Loadings Using BCFTools PGS
#'
#' Runs the BCFTools pgs plugin to compute polygenic score loadings.
#' This function converts GWAS summary statistics into polygenic score loadings
#' by accounting for linkage disequilibrium (LD) between variants.
#' Requires CHOLMOD library for sparse matrix operations.
#'
#' @param InputFileName Character; Path to GWAS-VCF file with summary statistics.
#' @param LDMatrix Character; Path to LDGM file with sparse LD matrix.
#' @param Regions Character; Restrict to comma-separated list of regions.
#' @param RegionsFile Character; Restrict to regions listed in file.
#' @param Targets Character; Similar to Regions but streams rather than index-jumps.
#' @param TargetsFile Character; Similar to RegionsFile but streams rather than index-jumps.
#' @param Prior Character; Prior specification for effect size variance.
#' @param LDScoreTag Character; INFO tag containing LD score values.
#' @param QScoreThreshold Numeric; Minimum quality score threshold for variants.
#' @param MaxFileSize Integer; Maximum file size in MB for memory mapping.
#' @param MaxChunkSize Integer; Maximum chunk size in MB for memory mapping.
#' @param AverageEffects Logical; Average effects across proxies for redundant sites.
#' @param AverageLDScore Numeric; Average genome-wide LD score for variants.
#' @param ExpectedRatio Numeric; Expected ratio of associations to null.
#' @param MAFThreshold Numeric; Remove variants with MAF below this threshold.
#' @param NoNormalize Logical; Do not normalize by allele frequency.
#' @param SampleNames Character; Comma-separated list of sample names to process.
#' @param SamplesFile Character; File with list of samples to include.
#' @param IncludeFilter Character; Include sites for which the expression is true.
#' @param ExcludeFilter Character; Exclude sites for which the expression is true.
#' @param OutputFile Character; Path to output file.
#' @param OutputType Character; b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF.
#' @param NumThreads Integer; Number of extra output compression threads.
#' @param WriteIndex Logical or Character; Automatically index the output file (optionally specify index format).
#' @param CatchStdout Logical; Capture standard output.
#' @param CatchStderr Logical; Capture standard error.
#' @param SaveStdout Character; Path to save standard output to.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{status}{Integer exit status of the command (0 for success, non-zero for errors)}
#'   \item{stdout}{Character vector of captured standard output lines, or NULL if not captured}
#'   \item{stderr}{Character vector of captured standard error lines, or NULL if not captured}
#'   \item{command}{Character vector representing the full command invoked}
#' }
#'
#' @note This function requires the CHOLMOD library for sparse matrix operations.
#' If the package was not compiled with CHOLMOD support, an error will be thrown.
#'
#' @export
BCFToolsPGS <- function(
  InputFileName,
  LDMatrix,
  Regions = NULL,
  RegionsFile = NULL,
  Targets = NULL,
  TargetsFile = NULL,
  Prior = NULL,
  LDScoreTag = NULL,
  QScoreThreshold = NULL,
  MaxFileSize = NULL,
  MaxChunkSize = NULL,
  AverageEffects = FALSE,
  AverageLDScore = NULL,
  ExpectedRatio = NULL,
  MAFThreshold = NULL,
  NoNormalize = FALSE,
  SampleNames = NULL,
  SamplesFile = NULL,
  IncludeFilter = NULL,
  ExcludeFilter = NULL,
  OutputFile = NULL,
  OutputType = NULL,
  NumThreads = NULL,
  WriteIndex = FALSE,
  CatchStdout = TRUE,
  CatchStderr = TRUE,
  SaveStdout = NULL
) {
  # Validate required parameters
  if (missing(InputFileName)) {
    stop("InputFileName is required")
  }

  if (missing(LDMatrix)) {
    stop("LDMatrix is required")
  }

  # Initialize arguments vector
  args <- character()

  # Build the command arguments
  args <- c(args, "--ldgm", LDMatrix)

  if (!is.null(Regions)) {
    args <- c(args, "--regions", Regions)
  }

  if (!is.null(RegionsFile)) {
    args <- c(args, "--regions-file", RegionsFile)
  }

  if (!is.null(Targets)) {
    args <- c(args, "--targets", Targets)
  }

  if (!is.null(TargetsFile)) {
    args <- c(args, "--targets-file", TargetsFile)
  }

  if (!is.null(Prior)) {
    args <- c(args, "--prior", Prior)
  }

  if (!is.null(LDScoreTag)) {
    args <- c(args, "--lds-tag", LDScoreTag)
  }

  if (!is.null(QScoreThreshold)) {
    args <- c(args, "--q-score-thr", as.character(QScoreThreshold))
  }

  if (!is.null(MaxFileSize)) {
    args <- c(args, "--max-file-size", as.character(MaxFileSize))
  }

  if (!is.null(MaxChunkSize)) {
    args <- c(args, "--max-chunk-size", as.character(MaxChunkSize))
  }

  if (AverageEffects) {
    args <- c(args, "--avg-effects")
  }

  if (!is.null(AverageLDScore)) {
    args <- c(args, "--avg-lds", as.character(AverageLDScore))
  }

  if (!is.null(ExpectedRatio)) {
    args <- c(args, "--er", as.character(ExpectedRatio))
  }

  if (!is.null(MAFThreshold)) {
    args <- c(args, "--min-maf", as.character(MAFThreshold))
  }

  if (NoNormalize) {
    args <- c(args, "--no-normalize")
  }

  if (!is.null(SampleNames)) {
    args <- c(args, "--samples", SampleNames)
  }

  if (!is.null(SamplesFile)) {
    args <- c(args, "--samples-file", SamplesFile)
  }

  if (!is.null(IncludeFilter)) {
    args <- c(args, "--include", IncludeFilter)
  }

  if (!is.null(ExcludeFilter)) {
    args <- c(args, "--exclude", ExcludeFilter)
  }

  if (!is.null(OutputFile)) {
    args <- c(args, "--output-file", OutputFile)
  }

  if (!is.null(OutputType)) {
    args <- c(args, "--output-type", OutputType)
  }

  if (!is.null(NumThreads)) {
    args <- c(args, "--threads", as.character(NumThreads))
  }

  if (is.logical(WriteIndex) && WriteIndex) {
    args <- c(args, "--write-index")
  } else if (is.character(WriteIndex)) {
    args <- c(args, paste0("--write-index=", WriteIndex))
  }

  # Add input file as last argument
  args <- c(args, InputFileName)

  # Create temporary files for stderr (and stdout if needed)
  stderrFile <- tempfile("bcftools_stderr_")
  stdoutFile <- if (is.null(SaveStdout)) {
    tempfile("bcftools_stdout_")
  } else {
    SaveStdout
  }

  # Call the unified pipeline C function
  status <- tryCatch(
    {
      pipeline_result <- .Call(
        RC_bcftools_pipeline,
        list("+pgs"), # Plugin command wrapped in list
        list(args), # Args wrapped in list
        1L, # Number of commands = 1
        CatchStdout,
        CatchStderr,
        stdoutFile,
        stderrFile
      )
      # Extract the single exit code
      pipeline_result[1]
    },
    error = function(e) {
      if (grepl("CHOLMOD", e$message)) {
        stop(
          "BCFTools was compiled without CHOLMOD support. Please install CHOLMOD and recompile.",
          call. = FALSE
        )
      }
      stop(e)
    }
  )

  # Process the results
  command <- attr(status, "command")

  # Read captured output if needed
  stdout <- NULL
  if (CatchStdout && file.exists(stdoutFile)) {
    stdout <- readLines(stdoutFile, warn = FALSE)
    if (!is.null(SaveStdout)) {
      # If we saved to a user-specified file, don't delete it
      if (length(stdout) == 0) stdout <- NULL
    } else {
      # Clean up temp file
      unlink(stdoutFile)
    }
  }

  stderr <- NULL
  if (CatchStderr && file.exists(stderrFile)) {
    stderr <- readLines(stderrFile, warn = FALSE)
    # Clean up temp file
    unlink(stderrFile)
    if (length(stderr) == 0) stderr <- NULL
  }

  # Return results
  list(
    status = status,
    stdout = stdout,
    stderr = stderr,
    command = command
  )
}
