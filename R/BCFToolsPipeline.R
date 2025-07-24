#' Create a pipeline of multiple BCFTools commands
#'
#' Executes a sequence of bcftools commands in a pipeline, where the output of each command
#' is used as input for the next command, similar to chaining multiple Unix pipe operators (|).
#'
#' @param ... A series of command specifications, where each specification consists of:
#'   - A character string specifying the bcftools command
#'   - A character vector of arguments for the command
#' @param catchStdout Logical; whether to capture standard output from the last command (default: TRUE)
#' @param catchStderr Logical; whether to capture standard error from all commands (default: TRUE)
#' @param saveStdout Character; file path where to save standard output, or NULL (default: NULL)
#'
#' @return A named list with elements:
#' \describe{
#'   \item{status}{Integer vector with exit statuses of all commands (0 for success, non-zero for errors)}
#'   \item{stdout}{Character vector of captured standard output lines from the last command, or NULL if not captured}
#'   \item{stderr}{Character vector of captured standard error lines from all commands, or NULL if not captured}
#'   \item{command}{Character vector representing the full piped bcftools command sequence invoked}
#' }
#'
#' @examples
#' \dontrun{
#' # Pipeline of three commands: view, filter, annotate
#' vcfFile <- system.file("exdata", "imputed.gt.vcf.gz", package = "RBCFLib")
#' result <- BCFToolsPipeline(
#'   "view", c("-r", "chr1:1000-2000", vcfFile), 
#'   "view", c("-i", "QUAL>20"),
#'   "annotate", c("-x", "INFO")
#' )
#'
#' # Filter, sort and index, saving to a file
#' outFile <- tempfile(fileext = ".vcf")
#' BCFToolsPipeline(
#'   "view", c("-i", "'DP>30'", vcfFile),
#'   "sort", character(),
#'   "view", c("-Oz"),
#'   saveStdout = outFile
#' )
#' }
#'
#' @export
BCFToolsPipeline <- function(...,
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

    # Collect arguments
    args <- list(...)
    if (length(args) < 2 || length(args) %% 2 != 0) {
        stop("Arguments must be pairs of command and argument vectors")
    }

    n_commands <- length(args) / 2
    commands <- vector("list", n_commands)
    command_args <- vector("list", n_commands)

    # Process each command-args pair
    for (i in 1:n_commands) {
        cmd_idx <- (i - 1) * 2 + 1
        args_idx <- cmd_idx + 1
        
        # Get command
        command <- args[[cmd_idx]]
        if (!is.character(command) || length(command) != 1) {
            stop(sprintf("Command %d must be a single character string", i))
        }
        
        if (!command %in% validCommands) {
            stop(sprintf(
                "'%s' is not a recognized bcftools command. Valid commands are: %s",
                command, paste(validCommands, collapse = ", ")
            ))
        }
        
        # Get arguments
        cmd_args <- args[[args_idx]]
        if (!is.character(cmd_args)) {
            stop(sprintf("Arguments for command %d must be a character vector", i))
        }
        
        # Store command and args
        commands[[i]] <- command
        command_args[[i]] <- cmd_args
    }

    # Validate output parameters
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

    # Create temporary files for stdout/stderr capture if needed
    stdout_file <- if (is.null(saveStdout)) tempfile() else saveStdout
    stderr_file <- tempfile()

    # Call the C function
    status <- .Call(RC_bcftools_pipeline, 
        as.list(commands),
        command_args,
        as.integer(n_commands),
        catchStdout,
        catchStderr,
        stdout_file,
        stderr_file,
        PACKAGE = "RBCFLib"
    )

    # Read captured output if needed
    stdout_lines <- NULL
    stderr_lines <- NULL
    
    if (catchStdout && is.null(saveStdout)) {
        if (file.exists(stdout_file) && file.info(stdout_file)$size > 0) {
            stdout_lines <- readLines(stdout_file)
            file.remove(stdout_file)
        }
    }
    
    if (catchStderr) {
        if (file.exists(stderr_file) && file.info(stderr_file)$size > 0) {
            stderr_lines <- readLines(stderr_file)
            file.remove(stderr_file)
        }
    }

    # Build the result
    result <- list(
        status = status,
        stdout = stdout_lines,
        stderr = stderr_lines,
        command = attr(status, "command")
    )

    return(result)
}
