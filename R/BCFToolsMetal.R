#' Run Metal-like Meta-analysis on GWAS VCF Files
#'
#' Runs BCFTools metal plugin to perform meta-analysis of GWAS summary statistics
#' in GWAS-VCF format. This function is modeled after the original METAL tool
#' for meta-analyzing genome-wide association studies.
#'
#' @param InputFileNames Character; List of input GWAS-VCF files to meta-analyze.
#' @param WeightsFile Character; Path to file with sample weights.
#' @param Scheme Character; Meta-analysis scheme (fixed-effects or random-effects).
#' @param Heterogeneity Character; Heterogeneity test to use.
#' @param OutlierThreshold Numeric; Maximum number of standard deviations for outlier detection.
#' @param FreqImputation Character; Method for frequency imputation (none, meta, hapmap).
#' @param FreqImputationMinMAF Numeric; Minimum MAF for frequency imputation.
#' @param SampleName Character; Sample name for the output file.
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
#' @export
BCFToolsMetal <- function(
  InputFileNames,
  WeightsFile = NULL,
  Scheme = NULL,
  Heterogeneity = NULL,
  OutlierThreshold = NULL,
  FreqImputation = NULL,
  FreqImputationMinMAF = NULL,
  SampleName = NULL,
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
  if (missing(InputFileNames)) {
    stop("InputFileNames is required")
  }

  # Initialize arguments vector
  args <- character()

  # Build the command arguments
  if (!is.null(WeightsFile)) {
    args <- c(args, "--weights", WeightsFile)
  }

  if (!is.null(Scheme)) {
    args <- c(args, "--scheme", Scheme)
  }

  if (!is.null(Heterogeneity)) {
    args <- c(args, "--heterogeneity", Heterogeneity)
  }

  if (!is.null(OutlierThreshold)) {
    args <- c(args, "--outlier", as.character(OutlierThreshold))
  }

  if (!is.null(FreqImputation)) {
    args <- c(args, "--freq-imputation", FreqImputation)
  }

  if (!is.null(FreqImputationMinMAF)) {
    args <- c(
      args,
      "--freq-imputation-min-maf",
      as.character(FreqImputationMinMAF)
    )
  }

  if (!is.null(SampleName)) {
    args <- c(args, "--sample-name", SampleName)
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

  # Add input files as last arguments
  if (is.character(InputFileNames)) {
    if (length(InputFileNames) == 1) {
      args <- c(args, InputFileNames)
    } else {
      args <- c(args, InputFileNames)
    }
  }

  # Create temporary files for stderr (and stdout if needed)
  stderrFile <- tempfile("bcftools_stderr_")
  stdoutFile <- if (is.null(SaveStdout)) tempfile("bcftools_stdout_") else
    SaveStdout

  # Call the unified pipeline C function
  pipeline_result <- .Call(
    RC_bcftools_pipeline,
    list("+metal"),         # Plugin command wrapped in list
    list(args),             # Args wrapped in list
    1L,                     # Number of commands = 1
    CatchStdout,
    CatchStderr,
    stdoutFile,
    stderrFile
  )
  
  # Extract the single exit code and command attribute
  status <- pipeline_result[1]
  command <- attr(pipeline_result, "command")

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
