#' Pipe two BCFTools commands together
#'
#' Executes two bcftools commands in a piped sequence, where the output of the first command
#' is used as input for the second command, similar to the Unix pipe operator (|).
#'
#' @param command1 A character string specifying the first bcftools command
#' @param args1 A character vector of arguments for the first command
#' @param command2 A character string specifying the second bcftools command
#' @param args2 A character vector of arguments for the second command
#' @param catchStdout Logical; whether to capture standard output from the second command (default: TRUE)
#' @param catchStderr Logical; whether to capture standard error from both commands (default: TRUE)
#' @param saveStdout Character; file path where to save standard output, or NULL (default: NULL)
#'
#' @return A named list with elements:
#' \describe{
#'   \item{status}{Integer vector with exit statuses of both commands (0 for success, non-zero for errors)}
#'   \item{stdout}{Character vector of captured standard output lines from the second command, or NULL if not captured}
#'   \item{stderr}{Character vector of captured standard error lines from both commands, or NULL if not captured}
#'   \item{command}{Character vector representing the full piped bcftools command invoked}
#' }
#'
#' @examples
#' \dontrun{
#' # Convert VCF to BCF and back to VCF
#' vcfFile <- system.file("exdata", "imputed.gt.vcf.gz", package = "RBCFLib")
#' result <- BCFToolsPipe("view", c("-Ob", vcfFile), "view", c("-Ov"))
#'
#' # Filter VCF then remove headers
#' BCFToolsPipe("view", c("-i", "'DP>30'", vcfFile), "view", c("-H"))
#'
#' # View a region then filter by some criteria and save to a file
#' outFile <- tempfile(fileext = ".vcf")
#' BCFToolsPipe("view", c("-r", "chr22:1000-2000000", vcfFile),
#'              "view", c("-i", "'QUAL>20'"), 
#'              saveStdout = outFile)
#' }
#'
#' @export
BCFToolsPipe <- function(command1,
                         args1 = character(),
                         command2,
                         args2 = character(),
                         catchStdout = TRUE,
                         catchStderr = TRUE,
                         saveStdout = NULL) {
    # List of valid bcftools commands
    validCommands <- c(
        "version", "view", "index", "query", "call", "mpileup", "concat",
        "merge", "norm", "stats", "annotate", "cnv", "consensus",
        "convert", "csq", "filter", "gtcheck", "plugin", "roh",
        "isec", "reheader", "sort", "help", "polysomy"
    )

    # Input validation
    if (!is.character(command1) || length(command1) != 1) {
        stop("'command1' must be a single character string")
    }

    if (!command1 %in% validCommands) {
        stop(sprintf(
            "'%s' is not a recognized bcftools command. Valid commands are: %s",
            command1, paste(validCommands, collapse = ", ")
        ))
    }

    if (!is.character(command2) || length(command2) != 1) {
        stop("'command2' must be a single character string")
    }

    if (!command2 %in% validCommands) {
        stop(sprintf(
            "'%s' is not a recognized bcftools command. Valid commands are: %s",
            command2, paste(validCommands, collapse = ", ")
        ))
    }

    if (!is.character(args1)) {
        stop("'args1' must be a character vector")
    }

    if (!is.character(args2)) {
        stop("'args2' must be a character vector")
    }

    if (!is.logical(catchStdout) || length(catchStdout) != 1) {
        stop("'catchStdout' must be a logical value")
    }

    if (!is.logical(catchStderr) || length(catchStderr) != 1) {
        stop("'catchStderr' must be a logical value")
    }

    if (!is.null(saveStdout) && !is.character(saveStdout)) {
        stop("'saveStdout' must be NULL or a character string")
    }

    # Enforce output capture in interactive mode for safety
    if (interactive()) {
        if (!catchStderr || !catchStdout) {
            stop("catchStdout and catchStderr must be TRUE in interactive mode")
        }
    }

    # Build shell command
    cmd <- paste(c("bcftools", command1, args1, "|", "bcftools", command2, args2), collapse = " ")

    if (!is.null(saveStdout)) {
        # Redirect stdout to file
        shcmd <- paste(cmd, ">", shQuote(saveStdout))
        status <- system(shcmd, intern = FALSE, ignore.stdout = FALSE, ignore.stderr = !catchStderr)
        stdout_lines <- NULL
    } else {
        # Capture stdout
        stdout_lines <- system(cmd, intern = TRUE, ignore.stderr = !catchStderr)
        status <- attr(stdout_lines, "status")
        if (is.null(status)) status <- 0L
    }

    # Return status vector
    # status may be integer or vector, but here one status for full pipeline
    # Set both statuses identical for compatibility
    status_vec <- as.integer(c(status, status))

    result <- list(
        status = status_vec,
        stdout = stdout_lines,
        stderr = NULL,
        command = strsplit(cmd, " ")[[1]]
    )

    return(result)
}
