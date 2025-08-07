#' Compute Best Linear Unbiased Prediction Using BCFTools BLUP
#'
#' Runs the BCFTools blup plugin to compute best linear unbiased prediction.
#' This function allows applying BLUP to GWAS summary statistics to
#' generate improved polygenic score weights.
#'
#' @param InputFileName Character; Path to the VCF/BCF file containing GWAS summary statistics.
#' @param LDMatrix Character; Path to the LD matrix file.
#' @param Regions Character; Restrict to comma-separated list of regions.
#' @param RegionsFile Character; Restrict to regions listed in file.
#' @param Targets Character; Similar to Regions but streams rather than index-jumps.
#' @param TargetsFile Character; Similar to RegionsFile but streams rather than index-jumps.
#' @param Prior Character; Prior specification for effect size variance.
#' @param LDScoreTag Character; INFO tag containing LD score values.
#' @param QScoreThreshold Numeric; Apply weights only if quality score exceeds threshold.
#' @param MaxFileSize Integer; Maximum file size in MB for memory mapping.
#' @param MaxChunkSize Integer; Maximum chunk size in MB for memory mapping.
#' @param AverageEffects Logical; Average effects across proxies for redundant sites.
#' @param AverageLDScore Numeric; Average genome-wide LD score for variants.
#' @param ExpectedRatio Numeric; Expected ratio of associations to null.
#' @param MAFThreshold Numeric; Remove variants with MAF below this threshold.
#' @param NoNormalize Logical; Do not normalize by allele frequency.
#' @param IncludeFilter Character; Include sites for which the expression is true.
#' @param ExcludeFilter Character; Exclude sites for which the expression is true.
#' @param OutputFile Character; Path to output file.
#' @param OutputType Character; b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF.
#' @param NumThreads Integer; Number of extra output compression threads.
#' @param WriteIndex Logical; Automatically index the output file.
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
BCFToolsBLUP <- function(
  InputFileName,
  LDMatrix = NULL,
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

  # Initialize arguments vector
  args <- character()

  # Build the command arguments
  if (!is.null(LDMatrix)) {
    args <- c(args, "--ldgm", LDMatrix)
  }

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

  if (is.logical(AverageEffects) && AverageEffects) {
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

  if (is.logical(NoNormalize) && NoNormalize) {
    args <- c(args, "--no-normalize")
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
  }

  # Add input file as last argument
  args <- c(args, InputFileName)

  # Create temporary files for stderr (and stdout if needed)
  stderrFile <- tempfile("bcftools_stderr_")
  stdoutFile <- if (is.null(SaveStdout)) tempfile("bcftools_stdout_") else
    SaveStdout

  # Call the C function
  status <- .Call(
    RC_bcftools_blup,
    args,
    CatchStdout,
    CatchStderr,
    stdoutFile,
    stderrFile
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
