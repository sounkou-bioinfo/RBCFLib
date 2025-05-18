#' Convert summary statistics to GWAS-VCF
#'
#' Converts GWAS summary statistics to the standardized GWAS-VCF format.
#' This function wraps the BCFTools munge plugin which transforms various summary
#' statistics formats into a standardized VCF format following the GWAS-VCF specification.
#'
#' @param InputFileName Character; path to input summary statistics file
#' @param Columns Character; column headers preset (PLINK/PLINK2/REGENIE/SAIGE/BOLT/METAL/PGS/SSF)
#' @param ColumnsFile Character; path to file with column headers definitions
#' @param FastaRef Character; path to reference sequence in FASTA format
#' @param FaiFile Character; path to the reference sequence index (.fai) file
#' @param CacheSize Integer; size of FASTA cache in bytes
#' @param IffyTag Character; FILTER annotation tag for undetermined reference alleles
#' @param MismatchTag Character; FILTER annotation tag for reference-allele mismatches
#' @param SampleName Character; sample name for the phenotype
#' @param NumSamples Numeric; number of samples
#' @param NumCases Numeric; number of cases (for case-control studies)
#' @param EffSampleSize Numeric; effective sample size
#' @param NoVersion Logical; whether to omit version and command line in the header
#' @param OutputFile Character; path to output file
#' @param OutputType Character; output format (u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level)
#' @param NumThreads Integer; number of worker threads for multithreading
#' @param WriteIndex Character or Logical; whether to automatically index output files and format to use
#' @param CatchStdout Logical; whether to capture standard output (default: TRUE)
#' @param CatchStderr Logical; whether to capture standard error (default: TRUE)
#' @param SaveStdout Character; file path where to save standard output, or NULL (default: NULL)
#'
#' @return A named list with elements:
#' \describe{
#'   \item{status}{Integer exit status of the command (0 for success, non-zero for errors)}
#'   \item{stdout}{Character vector of captured standard output lines, or NULL if not captured}
#'   \item{stderr}{Character vector of captured standard error lines, or NULL if not captured}
#'   \item{command}{Character vector representing the full command invoked}
#' }
#'
#' @examples
#' \dontrun{
#' # Convert PLINK format summary statistics to BCF
#' fastaRef <- system.file("exdata", "Test.fa", package = "RBCFLib")
#' inputFile <- system.file("exdata", "test_plink.tsv", package = "RBCFLib")
#' outputFile <- tempfile(fileext = ".bcf")
#'
#' BCFToolsMunge(
#'     InputFileName = inputFile,
#'     Columns = "PLINK",
#'     FastaRef = fastaRef,
#'     OutputFile = outputFile,
#'     OutputType = "b"
#' )
#'
#' # Convert using custom column headers
#' colHeaders <- system.file("exdata", "colheaders.tsv", package = "RBCFLib")
#' BCFToolsMunge(
#'     InputFileName = inputFile,
#'     ColumnsFile = colHeaders,
#'     FastaRef = fastaRef,
#'     SampleName = "STUDY_2023",
#'     OutputFile = outputFile,
#'     OutputType = "z"
#' )
#' }
#'
#' @export
BCFToolsMunge <- function(InputFileName,
                          Columns = NULL,
                          ColumnsFile = NULL,
                          FastaRef = NULL,
                          FaiFile = NULL,
                          CacheSize = NULL,
                          IffyTag = NULL,
                          MismatchTag = NULL,
                          SampleName = NULL,
                          NumSamples = NULL,
                          NumCases = NULL,
                          EffSampleSize = NULL,
                          NoVersion = FALSE,
                          OutputFile = NULL,
                          OutputType = NULL,
                          NumThreads = NULL,
                          WriteIndex = FALSE,
                          CatchStdout = TRUE,
                          CatchStderr = TRUE,
                          SaveStdout = NULL) {
    # Validate essential parameters
    if (missing(InputFileName)) {
        stop("InputFileName is required")
    }

    if (is.null(Columns) && is.null(ColumnsFile)) {
        stop("Either Columns or ColumnsFile must be provided")
    }

    if (!is.null(Columns) && !is.null(ColumnsFile)) {
        stop("Only one of Columns or ColumnsFile should be provided, not both")
    }

    if (is.null(FastaRef) && is.null(FaiFile)) {
        stop("Either FastaRef or FaiFile must be provided")
    }

    # Initialize arguments vector
    args <- character()

    # Build the command arguments
    if (!is.null(Columns)) {
        args <- c(args, "-c", Columns)
    }

    if (!is.null(ColumnsFile)) {
        args <- c(args, "-C", ColumnsFile)
    }

    if (!is.null(FastaRef)) {
        args <- c(args, "-f", FastaRef)
    }

    if (!is.null(FaiFile)) {
        args <- c(args, "--fai", FaiFile)
    }

    if (!is.null(CacheSize)) {
        args <- c(args, "--set-cache-size", as.character(CacheSize))
    }

    if (!is.null(IffyTag)) {
        args <- c(args, "--iffy-tag", IffyTag)
    }

    if (!is.null(MismatchTag)) {
        args <- c(args, "--mismatch-tag", MismatchTag)
    }

    if (!is.null(SampleName)) {
        args <- c(args, "-s", SampleName)
    }

    if (!is.null(NumSamples)) {
        args <- c(args, "--ns", as.character(NumSamples))
    }

    if (!is.null(NumCases)) {
        args <- c(args, "--nc", as.character(NumCases))
    }

    if (!is.null(EffSampleSize)) {
        args <- c(args, "--ne", as.character(EffSampleSize))
    }

    if (NoVersion) {
        args <- c(args, "--no-version")
    }

    if (!is.null(OutputFile)) {
        args <- c(args, "-o", OutputFile)
    }

    if (!is.null(OutputType)) {
        args <- c(args, "-O", OutputType)
    }

    if (!is.null(NumThreads)) {
        args <- c(args, "--threads", as.character(NumThreads))
    }

    if (is.logical(WriteIndex) && WriteIndex) {
        args <- c(args, "-W")
    } else if (is.character(WriteIndex)) {
        args <- c(args, paste0("-W=", WriteIndex))
    }

    # Add input file as last argument
    args <- c(args, InputFileName)

    # Create temporary files for stderr (and stdout if needed)
    stderrFile <- tempfile("bcftools_stderr_")
    stdoutFile <- if (is.null(SaveStdout)) tempfile("bcftools_stdout_") else SaveStdout

    # Call the C function
    status <- .Call(
        RC_bcftools_munge,
        args,
        CatchStdout,
        CatchStderr,
        stdoutFile,
        stderrFile,
        FALSE
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
