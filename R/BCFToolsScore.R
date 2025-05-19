#' Compute Polygenic Scores Using BCFTools Score
#'
#' Runs the BCFTools score plugin to compute polygenic scores.
#' This function allows applying weights for polygenic score calculation
#' to genotype data in VCF/BCF format.
#'
#' @param InputFileName Character; Path to the VCF/BCF file containing genotypes.
#' @param ScoresFile Character; Path to the file containing effect sizes/weights for variants.
#' @param SamplesFile Character; Path to file containing list of samples to include.
#' @param Regions Character; Restrict to comma-separated list of regions.
#' @param RegionsFile Character; Restrict to regions listed in file.
#' @param Targets Character; Similar to Regions but streams rather than index-jumps.
#' @param TargetsFile Character; Similar to RegionsFile but streams rather than index-jumps.
#' @param Samples Character; List of samples to include.
#' @param Format Character; Format field to use for scoring (GT/DS).
#' @param ScoresColumn Character; Column name or number (1-based) in ScoresFile containing weights.
#' @param OutputFile Character; Path to output file.
#' @param OutputType Character; b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF.
#' @param OutputColumns Character; Comma-separated list of columns to output.
#' @param NumThreads Integer; Number of extra output compression threads.
#' @param WriteIndex Logical; Automatically index the output file.
#' @param TSV Logical; Force output in TSV format.
#' @param IncludeFilter Character; Include sites for which the expression is true.
#' @param ExcludeFilter Character; Exclude sites for which the expression is true.
#' @param VariantID Logical; Use variant IDs instead of coordinates for alignment.
#' @param QScoreThreshold Numeric; Apply weights only if quality score exceeds threshold.
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
BCFToolsScore <- function(InputFileName,
                          ScoresFile = NULL,
                          SamplesFile = NULL,
                          Regions = NULL,
                          RegionsFile = NULL,
                          Targets = NULL,
                          TargetsFile = NULL,
                          Samples = NULL,
                          Format = NULL,
                          ScoresColumn = NULL,
                          OutputFile = NULL,
                          OutputType = NULL,
                          OutputColumns = NULL,
                          NumThreads = NULL,
                          WriteIndex = FALSE,
                          TSV = FALSE,
                          IncludeFilter = NULL,
                          ExcludeFilter = NULL,
                          VariantID = FALSE,
                          QScoreThreshold = NULL,
                          CatchStdout = TRUE,
                          CatchStderr = TRUE,
                          SaveStdout = NULL) {
    # Validate required parameters
    if (missing(InputFileName)) {
        stop("InputFileName is required")
    }

    # Initialize arguments vector
    args <- character()

    # Build the command arguments
    if (!is.null(ScoresFile)) {
        args <- c(args, "--scores", ScoresFile)
    }

    if (!is.null(SamplesFile)) {
        args <- c(args, "--samples-file", SamplesFile)
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

    if (!is.null(Samples)) {
        args <- c(args, "--samples", Samples)
    }

    if (!is.null(Format)) {
        args <- c(args, "--format", Format)
    }

    if (!is.null(ScoresColumn)) {
        args <- c(args, "--score-column", as.character(ScoresColumn))
    }

    if (!is.null(OutputFile)) {
        args <- c(args, "--output-file", OutputFile)
    }

    if (!is.null(OutputType)) {
        args <- c(args, "--output-type", OutputType)
    }

    if (!is.null(OutputColumns)) {
        args <- c(args, "--columns", OutputColumns)
    }

    if (!is.null(NumThreads)) {
        args <- c(args, "--threads", as.character(NumThreads))
    }

    if (is.logical(WriteIndex) && WriteIndex) {
        args <- c(args, "--write-index")
    }

    if (TSV) {
        args <- c(args, "--tsv")
    }

    if (!is.null(IncludeFilter)) {
        args <- c(args, "--include", IncludeFilter)
    }

    if (!is.null(ExcludeFilter)) {
        args <- c(args, "--exclude", ExcludeFilter)
    }

    if (VariantID) {
        args <- c(args, "--use-id")
    }

    if (!is.null(QScoreThreshold)) {
        args <- c(args, "--q-score-thr", as.character(QScoreThreshold))
    }

    # Add input file as last argument
    args <- c(args, InputFileName)

    # Create temporary files for stderr (and stdout if needed)
    stderrFile <- tempfile("bcftools_stderr_")
    stdoutFile <- if (is.null(SaveStdout)) tempfile("bcftools_stdout_") else SaveStdout

    # Call the C function
    status <- .Call(
        RC_bcftools_score,
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
