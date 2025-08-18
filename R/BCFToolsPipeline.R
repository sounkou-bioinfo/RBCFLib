#' Create a pipeline of multiple BCFTools commands
#'
#' Executes a sequence of bcftools commands in a pipeline, where the output of each command
#' is used as input for the next command, similar to chaining multiple Unix pipe operators (|).
#'
#' The function automatically validates that -o/--output options are only used in the last
#' command of the pipeline and only for commands that support output redirection.
#'
#' @param ... A series of command specifications, where each specification consists of:
#'   - A character string specifying the bcftools command
#'   - A character vector of arguments for the command
#' @param catchStdout Logical; whether to capture standard output from the last command (default: TRUE)
#' @param catchStderr Logical; whether to capture standard error from all commands (default: TRUE)
#' @param saveStdout Character; file path where to save standard output, or NULL (default: NULL)
#'
#' @details
#' **Output Handling:**
#' - Only the last command in the pipeline can contain -o/--output/--output-file arguments
#' - Commands that don't support -o/--output options:  index, roh, stats
#' - If -o/--output is used in any non-final command, an error will be thrown
#' - If -o/--output is used with unsupported commands, an error will be thrown
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
#' # Filter, sort and output to file (only last command has -o)
#' outFile <- tempfile(fileext = ".vcf.gz")
#' BCFToolsPipeline(
#'   "view", c("-i", "'DP>30'", vcfFile),
#'   "sort", character(),
#'   "view", c("-Oz", "-o", outFile)
#' )
#'
#' # INVALID: -o in non-final command (will throw error)
#' # BCFToolsPipeline(
#' #   "view", c("-o", "temp.vcf", vcfFile),  # ERROR: -o not allowed here
#' #   "sort", character()
#' # )
#'
#' # INVALID: -o with unsupported command (will throw error)
#' # BCFToolsPipeline(
#' #   "view", c(vcfFile),
#' #   "stats", c("-o", "stats.txt")  # ERROR: stats doesn't support -o
#' # )
#' }
#'
#' @export
BCFToolsPipeline <- function(
    ...,
    catchStdout = TRUE,
    catchStderr = TRUE,
    saveStdout = NULL) {
  # List of valid bcftools commands
  validCommands <- c(
    "version",
    "view",
    "index",
    "query",
    "call",
    "mpileup",
    "concat",
    "merge",
    "norm",
    "stats",
    "annotate",
    "cnv",
    "consensus",
    "convert",
    "csq",
    "filter",
    "gtcheck",
    "plugin",
    "roh",
    "isec",
    "reheader",
    "sort",
    "head",
    "help",
    "+GTisec", "+GTsubset", "+ad-bias", "+add-variantkey", "+af-dist", "+allele-length", "+blup", "+check-ploidy", "+check-sparsity", "+color-chrs", "+contrast", "+counts", "+dosage", "+fill-AN-AC", "+fill-from-fasta", "+fill-tags", "+fixploidy", "+fixref", "+frameshifts", "+guess-ploidy", "+gvcfz", "+impute-info", "+indel-stats", "+isecGT", "+liftover", "+mendelian2", "+metal", "+missing2ref", "+munge", "+parental-origin", "+prune", "+remove-overlaps", "+scatter", "+score", "+setGT", "+smpl-stats", "+split", "+split-vep", "+tag2tag", "+trio-dnm2", "+trio-stats", "+trio-switch-rate", "+variant-distance", "+variantkey-hex", "+vcf2table", "+vrfs"
  )

  # Commands that don't support standard output redirection with -o option
  # Based on the pysam implementation and BCFToolsRun.R
  # head is for samtools, maybe we will wrap it too
  EXCLUDED_COMMANDS <- c("head", "index", "roh", "stats", "+guess-ploidy")

  # Collect arguments
  # TODO this is brittle, we should make a proper DSL
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
        command,
        paste(validCommands, collapse = ", ")
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

  # Validate -o/--output arguments in pipeline
  # Function to check if arguments contain output option
  has_output_option <- function(args) {
    for (arg in args) {
      if (
        arg == "-o" ||
          arg == "--output" ||
          arg == "--output-file" ||
          startsWith(arg, "-o=") ||
          startsWith(arg, "--output=") ||
          startsWith(arg, "--output-file=")
      ) {
        return(TRUE)
      }
    }
    return(FALSE)
  }

  # Check each command for output options
  for (i in 1:n_commands) {
    cmd <- commands[[i]]
    args <- command_args[[i]]

    if (has_output_option(args)) {
      # Only the last command can have output option
      if (i < n_commands) {
        stop(sprintf(
          "Command %d ('%s') contains -o/--output/--output-file option, but only the last command in a pipeline can specify output",
          i,
          cmd
        ))
      }

      # Check if the command supports output option
      if (cmd %in% EXCLUDED_COMMANDS) {
        stop(sprintf(
          "Command '%s' does not support -o/--output/--output-file option. Commands that don't support output redirection: %s",
          cmd,
          paste(EXCLUDED_COMMANDS, collapse = ", ")
        ))
      }
    }
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
  status <- .Call(
    RC_bcftools_pipeline,
    as.list(commands),
    command_args,
    as.integer(n_commands),
    catchStdout,
    catchStderr,
    stdout_file,
    stderr_file,
    PACKAGE = "RBCFLib"
  )

  # Function to collect output from file (handles both text and binary)
  collect_output <- function(file_path) {
    if (!file.exists(file_path) || file.info(file_path)$size == 0) {
      return(character(0))
    }

    result <- tryCatch(
      {
        con <- file(file_path, "r", encoding = "UTF-8")
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

  # Read captured output if needed
  stdout_lines <- NULL
  stderr_lines <- NULL

  # Only read stdout content if catchStdout is TRUE AND saveStdout is NULL
  # This matches BCFToolsRun behavior
  if (catchStdout && is.null(saveStdout)) {
    stdout_lines <- collect_output(stdout_file)
    file.remove(stdout_file)
  }

  if (catchStderr) {
    stderr_lines <- collect_output(stderr_file)
    file.remove(stderr_file)
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
