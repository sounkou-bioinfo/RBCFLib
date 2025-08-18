#' Liftover Variants Between Genome Assemblies
#'
#' Runs BCFTools liftover to convert variants coordinates between genome assemblies.
#' This function provides a way to remap genomic coordinates from one genome assembly to another.
#'
#' @param InputFileName Character; Path to input VCF/BCF file with variants to lift over.
#' @param ChainFile Character; Path to chain file that maps old assembly to new assembly.
#' @param FastaRef Character; Path to reference sequence in FASTA format.
#' @param Regions Character; Restrict to comma-separated list of regions.
#' @param RegionsFile Character; Restrict to regions listed in file.
#' @param Targets Character; Similar to Regions but streams rather than index-jumps.
#' @param TargetsFile Character; Similar to RegionsFile but streams rather than index-jumps.
#' @param FlipTag Character; INFO tag to mark reverse-complemented sites.
#' @param SwapTag Character; INFO tag to mark swapped REF/ALT sites.
#' @param DropTags Character; Comma-separated list of tags to drop (incompatible with lifting over).
#' @param ACTags Character; Comma-separated list of allele count tags to adjust.
#' @param AFTags Character; Comma-separated list of allele frequency tags to adjust.
#' @param DSTags Character; Comma-separated list of dosage tags to adjust.
#' @param GTTags Character; Comma-separated list of genotype tags to adjust.
#' @param ESTags Character; Comma-separated list of effect size tags to adjust.
#' @param Tags Character; List of tags to convert (can be specified multiple times).
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
BCFToolsLiftover <- function(
  InputFileName,
  ChainFile,
  FastaRef = NULL,
  Regions = NULL,
  RegionsFile = NULL,
  Targets = NULL,
  TargetsFile = NULL,
  FlipTag = NULL,
  SwapTag = NULL,
  DropTags = NULL,
  ACTags = NULL,
  AFTags = NULL,
  DSTags = NULL,
  GTTags = NULL,
  ESTags = NULL,
  Tags = NULL,
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

  if (missing(ChainFile)) {
    stop("ChainFile is required")
  }

  # Initialize arguments vector
  args <- character()

  # Build the command arguments
  args <- c(args, "--chain", ChainFile)

  if (!is.null(FastaRef)) {
    args <- c(args, "--fasta-ref", FastaRef)
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

  if (!is.null(FlipTag)) {
    args <- c(args, "--flip-tag", FlipTag)
  }

  if (!is.null(SwapTag)) {
    args <- c(args, "--swap-tag", SwapTag)
  }

  if (!is.null(DropTags)) {
    args <- c(args, "--drop-tags", DropTags)
  }

  if (!is.null(ACTags)) {
    args <- c(args, "--ac-tags", ACTags)
  }

  if (!is.null(AFTags)) {
    args <- c(args, "--af-tags", AFTags)
  }

  if (!is.null(DSTags)) {
    args <- c(args, "--ds-tags", DSTags)
  }

  if (!is.null(GTTags)) {
    args <- c(args, "--gt-tags", GTTags)
  }

  if (!is.null(ESTags)) {
    args <- c(args, "--es-tags", ESTags)
  }

  if (!is.null(Tags)) {
    if (is.character(Tags)) {
      for (tag in Tags) {
        args <- c(args, "--tags", tag)
      }
    }
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
  stdoutFile <- if (is.null(SaveStdout)) tempfile("bcftools_stdout_") else
    SaveStdout

  # Call the unified pipeline C function
  pipeline_result <- .Call(
    RC_bcftools_pipeline,
    list("+liftover"),      # Plugin command wrapped in list
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
