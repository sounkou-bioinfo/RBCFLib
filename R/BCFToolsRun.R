#' Run BCFTools commands
#'
#' Executes bcftools commands from within R using the embedded bcftools library.
#' This function allows running any bcftools command and its options directly from R.
#' It automatically handles output redirection for commands that accept -o/--output options.
#'
#' @param command A character string specifying the bcftools command to run
#'   (e.g., "view", "query", "call", etc.)
#' @param args A character vector of arguments to pass to the bcftools command
#' @param catchStdout Logical; whether to capture standard output (default: TRUE)
#' @param catchStderr Logical; whether to capture standard error (default: TRUE)
#' @param saveStdout Character; file path where to save standard output, or NULL (default: NULL)
#' @param isUsage Logical; if TRUE, will not redirect output for commands with explicit output options (default: FALSE)
#'
#' @return A named list with elements:
#' \describe{
#'   \item{status}{Integer exit status of the bcftools command (0 for success, non-zero for errors)}
#'   \item{stdout}{Character vector of captured standard output lines, or NULL if not captured}
#'   \item{stderr}{Character vector of captured standard error lines, or NULL if not captured}
#'   \item{command}{Character vector representing the full bcftools command invoked}
#' }
#'
#' @examples
#' \dontrun{
#' # View a VCF file header
#' vcfFile <- system.file("exdata", "imputed.gt.vcf.gz", package = "RBCFLib")
#' BCFToolsRun("view", c("-h", vcfFile))
#'
#' # Query specific fields from a VCF file
#' BCFToolsRun("query", c("-f", "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n", vcfFile))
#'
#' # Save output directly to a file
#' outFile <- tempfile(fileext = ".vcf")
#' BCFToolsRun("view", c("-h", vcfFile), saveStdout = outFile)
#' }
#'
#' @export
BCFToolsRun <- function(
  command,
  args = character(),
  catchStdout = TRUE,
  catchStderr = TRUE,
  saveStdout = NULL,
  isUsage = FALSE
) {
  # Create temporary files for stderr (and stdout if needed)
  stderrFile <- tempfile("bcftools_stderr_")
  stdoutFile <- if (is.null(saveStdout)) tempfile("bcftools_stdout_") else
    saveStdout

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
    "polysomy",
    # plugin commands
    "+GTisec","+GTsubset","+ad-bias","+add-variantkey","+af-dist","+allele-length","+blup","+check-ploidy","+check-sparsity","+color-chrs","+contrast","+counts","+dosage","+fill-AN-AC","+fill-from-fasta","+fill-tags","+fixploidy","+fixref","+frameshifts","+guess-ploidy","+gvcfz","+impute-info","+indel-stats","+isecGT","+liftover","+mendelian2","+metal","+missing2ref","+munge","+parental-origin","+prune","+remove-overlaps","+scatter","+score","+setGT","+smpl-stats","+split","+split-vep","+tag2tag","+trio-dnm2","+trio-stats","+trio-switch-rate","+variant-distance","+variantkey-hex","+vcf2table","+vrfs"
  )

  # Input validation
  if (!is.character(command) || length(command) != 1) {
    stop("'command' must be a single character string")
  }

  if (!command %in% validCommands) {
    stop(sprintf(
      "'%s' is not a recognized bcftools command. Valid commands are: %s",
      command,
      paste(validCommands, collapse = ", ")
    ))
  }

  if (!is.character(args)) {
    stop("'args' must be a character vector")
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

  if (!is.logical(isUsage) || length(isUsage) != 1) {
    stop("'isUsage' must be a logical value")
  }

  # Enforce output capture in interactive mode for safety
  if (interactive()) {
    if (!catchStderr || !catchStdout) {
      stop("catchStdout and catchStderr must be TRUE in interactive mode")
    }
  }

  # Check for index command file existence
  if (command == "index" && length(args) > 0) {
    # Skip option arguments
    ARGUMENTS <- c(
      "-m",
      "--min-shift",
      "-o",
      "--output",
      "--output-file",
      "-@",
      "--threads"
    )
    skip_next <- FALSE
    for (arg in args) {
      if (skip_next) {
        skip_next <- FALSE
        next
      }
      if (startsWith(arg, "-")) {
        # Skip next argument for e.g. '--min-shift' '12' or '-m' '12' but not '-m12'
        if (arg %in% ARGUMENTS) {
          skip_next <- TRUE
        }
        next
      }
      if (!file.exists(arg)) {
        stop(sprintf("No such file or directory: '%s'", arg))
      } else {
        break
      }
    }
  }

  # Make a copy of args that we'll potentially modify
  working_args <- list(args)

  # Handle automatic output redirection similar to pysam dispatcher
  if (catchStdout) {
    # Commands that don't support standard output redirection with -o option
    # Based on the pysam implementation
    EXCLUDED_COMMANDS <- c("head", "index", "roh", "stats")

    if (!is.null(saveStdout)) {
      # User explicitly requested to save stdout to a file
      # Nothing special needed here, we'll use the file as provided
    } else if (catchStdout && !isUsage) {
      stdout_option <- NULL

      # In bcftools, most methods accept -o for output redirection
      if (!(command %in% EXCLUDED_COMMANDS)) {
        stdout_option <- "-o"

        # Check if the option is already present
        has_output_option <- FALSE
        for (i in seq_along(working_args[[1]])) {
          arg <- working_args[[1]][i]
          if (
            arg == "-o" ||
              arg == "--output" ||
              startsWith(arg, "-o=") ||
              startsWith(arg, "--output=")
          ) {
            has_output_option <- TRUE
            break
          }
        }

        # Add output option if not already present
        if (!has_output_option && !is.null(stdout_option)) {
          working_args[[1]] <- c(working_args[[1]], stdout_option, stdoutFile)
        }
      }
    }
  }

  # Call the C function using the unified pipeline interface
  pipeline_result <- .Call(
    RC_bcftools_pipeline,
    list(command),           # Single command wrapped in list
    list(working_args[[1]]), # Single args list wrapped in list  
    1L,                      # Number of commands = 1
    catchStdout,
    catchStderr,
    stdoutFile,
    stderrFile
  )
  
  # Extract the single exit code and command attribute
  res_int <- pipeline_result[1]
  attr(res_int, "command") <- attr(pipeline_result, "command")

  # Extract command attribute
  cmd <- attr(res_int, "command") # Function to collect output from file
  collect_output <- function(file_path) {
    if (!file.exists(file_path)) {
      return(character(0))
    }

    result <- tryCatch(
      {
        readLines(file_path)
      },
      error = function(e) {
        # Handle binary output by reading as raw if text reading fails
        readBin(file_path, what = "raw", n = file.info(file_path)$size)
      }
    )

    # Clean up temporary files, but preserve user-requested output files
    if (!is.null(saveStdout) && file_path != saveStdout) {
      tryCatch(file.remove(file_path), error = function(e) NULL)
    } else if (is.null(saveStdout)) {
      tryCatch(file.remove(file_path), error = function(e) NULL)
    }

    return(result)
  }

  # Read captured output
  stdout_lines <- if (catchStdout && is.null(saveStdout)) {
    collect_output(stdoutFile)
  } else {
    NULL
  }

  stderr_lines <- if (catchStderr) {
    collect_output(stderrFile)
  } else {
    NULL
  }

  # Build result list
  result <- list(
    status = as.integer(res_int),
    stdout = stdout_lines,
    stderr = stderr_lines,
    command = cmd
  )
  return(result)
}
