#' Run bcftools commands
#'
#' Executes bcftools commands from within R using the embedded bcftools library.
#' This function allows running any bcftools command and its options directly from R.
#'
#' @param args A character vector of arguments to pass to bcftools, where the first
#'   element should be the bcftools command (e.g., "view", "query", "call", etc.)
#' @param captureStdout Logical; whether to capture standard output to a file (default: FALSE)
#' @param captureStderr Logical; whether to capture standard error to a file (default: FALSE)
#' @param stdoutFile Character; file path where to save standard output (if captureStdout is TRUE)
#' @param stderrFile Character; file path where to save standard error (if captureStderr is TRUE)
#'
#' @return Integer value representing the exit status of the bcftools command
#'   (typically 0 for success, non-zero for errors)
#'
#' @examples
#' \dontrun{
#' # View a VCF file header
#' vcfFile <- system.file("exdata", "imputed.gt.vcf.gz", package = "RBCFLib")
#' BCFToolsRun(c("view", "-h", vcfFile))
#'
#' # Query specific fields from a VCF file
#' BCFToolsRun(c("query", "-f", "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n", vcfFile))
#' }
#'
#' @export
BCFToolsRun <- function(args,
                        captureStdout = FALSE,
                        captureStderr = FALSE,
                        stdoutFile = tempfile("bcftools_stdout_"),
                        stderrFile = tempfile("bcftools_stderr_")) {
    # Input validation
    if (!is.character(args)) {
        stop("'args' must be a character vector")
    }

    if (!is.logical(captureStdout) || length(captureStdout) != 1) {
        stop("'captureStdout' must be a logical value")
    }

    if (!is.logical(captureStderr) || length(captureStderr) != 1) {
        stop("'captureStderr' must be a logical value")
    }

    if (!is.character(stdoutFile) || length(stdoutFile) != 1) {
        stop("'stdoutFile' must be a character string")
    }

    if (!is.character(stderrFile) || length(stderrFile) != 1) {
        stop("'stderrFile' must be a character string")
    }

    # Call the C function
    status <- .Call(
        RC_bcftools_run,
        args,
        captureStdout,
        captureStderr,
        stdoutFile,
        stderrFile
    )

    # Return the status code
    return(status)
}
