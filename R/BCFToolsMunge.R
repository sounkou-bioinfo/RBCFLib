#' Run BCFTools Munge Command
#'
#' Executes bcftools munge command from within R using the embedded bcftools library.
#' This function is specifically designed for the munge command which converts GWAS
#' summary statistics to VCF format.
#'
#' @param input_file Path to input GWAS summary statistics TSV file
#' @param fasta_ref Path to reference sequence in FASTA format
#' @param output File path for output VCF (if NULL, will use a temporary file)
#' @param columns Column header preset (one of: PLINK, PLINK2, REGENIE, SAIGE, BOLT, METAL, PGS, SSF)
#' @param columns_file Path to tab-delimited file with column headers mapping
#' @param fai Path to reference sequence .fai index
#' @param set_cache_size Cache size in bytes for FASTA parsing
#' @param iffy_tag FILTER annotation tag for when reference allele is uncertain (default: "IFFY")
#' @param mismatch_tag FILTER annotation tag for when reference doesn't match any allele (default: "REF_MISMATCH")
#' @param sample_name Sample name for the phenotype (default: "SAMPLE")
#' @param ns Number of samples in the study
#' @param nc Number of cases in the study (for case-control studies)
#' @param args Additional command-line arguments as character vector
#' @param catchStdout Logical; whether to capture standard output (default: TRUE)
#' @param catchStderr Logical; whether to capture standard error (default: TRUE)
#' @param saveStdout Character; file path where to save standard output, or NULL (default: NULL)
#' @param isUsage Logical; if TRUE, will not redirect output (default: FALSE)
#'
#' @return A named list with elements:
#' \describe{
#'   \item{status}{Integer exit status of the bcftools munge command (0 for success, non-zero for errors)}
#'   \item{stdout}{Character vector of captured standard output lines, or NULL if not captured}
#'   \item{stderr}{Character vector of captured standard error lines, or NULL if not captured}
#'   \item{command}{Character vector representing the full bcftools munge command invoked}
#' }
#'
#' @examples
#' \dontrun{
#' # Convert GWAS summary statistics to VCF with basic parameters
#' inFile <- system.file("exdata", "test_plink.tsv", package = "RBCFLib")
#' refFile <- system.file("exdata", "Test.fa", package = "RBCFLib")
#' outFile <- tempfile(fileext = ".vcf")
#'
#' # Using the main parameters
#' result <- BCFToolsMunge(input_file = inFile, fasta_ref = refFile, output = outFile)
#'
#' # With column preset
#' result <- BCFToolsMunge(
#'     input_file = inFile,
#'     fasta_ref = refFile,
#'     output = outFile,
#'     columns = "PLINK"
#' )
#'
#' # With custom column mappings
#' colFile <- system.file("exdata", "colheaders.tsv", package = "RBCFLib")
#' result <- BCFToolsMunge(
#'     input_file = inFile,
#'     fasta_ref = refFile,
#'     output = outFile,
#'     columns_file = colFile
#' )
#' }
#'
#' @export
BCFToolsMunge <- function(input_file = NULL,
                          fasta_ref = NULL,
                          output = NULL,
                          columns = NULL,
                          columns_file = NULL,
                          fai = NULL,
                          set_cache_size = NULL,
                          iffy_tag = NULL,
                          mismatch_tag = NULL,
                          sample_name = NULL,
                          ns = NULL,
                          nc = NULL,
                          args = character(),
                          catchStdout = TRUE,
                          catchStderr = TRUE,
                          saveStdout = NULL,
                          isUsage = FALSE) {
    # Create temporary files for stderr (and stdout if needed)
    stderrFile <- tempfile("bcftools_munge_stderr_")
    stdoutFile <- if (is.null(saveStdout)) tempfile("bcftools_munge_stdout_") else saveStdout

    # Build command-line arguments from function parameters
    all_args <- character() # Start with empty arguments array, C function will add bcftools and +munge

    # Required parameters - fasta_ref must be first for option parsing to work correctly
    if (!is.null(fasta_ref)) {
        all_args <- c(all_args, "-f", fasta_ref)
    }

    # Optional parameters
    if (!is.null(columns)) {
        all_args <- c(all_args, "-c", columns)
    }

    if (!is.null(output)) {
        all_args <- c(all_args, "--output", output)
    }

    if (!is.null(columns_file)) {
        all_args <- c(all_args, "-C", columns_file)
    }

    if (!is.null(fai)) {
        all_args <- c(all_args, "--fai", fai)
    }

    if (!is.null(set_cache_size)) {
        all_args <- c(all_args, "--set-cache-size", as.character(set_cache_size))
    }

    if (!is.null(iffy_tag)) {
        all_args <- c(all_args, "--iffy-tag", iffy_tag)
    }

    if (!is.null(mismatch_tag)) {
        all_args <- c(all_args, "--mismatch-tag", mismatch_tag)
    }

    if (!is.null(sample_name)) {
        all_args <- c(all_args, "-s", sample_name)
    }

    if (!is.null(ns)) {
        all_args <- c(all_args, "--ns", as.character(ns))
    }

    if (!is.null(nc)) {
        all_args <- c(all_args, "--nc", as.character(nc))
    }

    # Add any additional arguments
    if (length(args) > 0) {
        all_args <- c(all_args, args)
    }

    # Input validation (most validations can stay here)
    if (!is.logical(catchStdout) || length(catchStdout) != 1) {
        stop("CatchStdout must be a logical value")
    }

    if (!is.logical(catchStderr) || length(catchStderr) != 1) {
        stop("CatchStderr must be a logical value")
    }

    if (!is.null(saveStdout) && !is.character(saveStdout)) {
        stop("SaveStdout must be NULL or a character string")
    }

    if (!is.logical(isUsage) || length(isUsage) != 1) {
        stop("IsUsage must be a logical value")
    }

    # Enforce output capture in interactive mode for safety
    if (interactive()) {
        if (!catchStderr || !catchStdout) {
            stop("catchStdout and catchStderr must be TRUE in interactive mode")
        }
    }

    # Prepare the final set of arguments for the C call
    # Start with all_args which includes user-specified options but not the input_file yet.
    final_call_args <- all_args

    # Handle automatic output redirection if output is not specified by user and we are capturing stdout
    if (catchStdout && !isUsage && is.null(saveStdout) && is.null(output)) {
        # If no output file is specified via the 'output' parameter,
        # and we are capturing stdout (and not saving it to a user file via saveStdout),
        # then bcftools munge needs an explicit --output. We use the temp stdoutFile.
        # The 'output' variable is also updated here to reflect this choice for consistency.
        output <- stdoutFile
        final_call_args <- c(final_call_args, "--output", stdoutFile)
    }

    # Now, add input file as the very last positional argument
    if (!is.null(input_file)) {
        final_call_args <- c(final_call_args, input_file)
    }

    # Validate the type of the final arguments list
    if (!is.character(final_call_args)) {
        stop("Arguments must be a character vector")
    }

    # Make a list for .Call as it expects a list or pairlist for the arguments parameter
    working_args <- list(final_call_args)

    # Call the C function, which returns an integer with 'command' attribute
    # Pass the arguments directly to the bcftools munge plugin
    res_int <- .Call(
        RC_bcftools_munge,
        working_args[[1]], # This now contains the correctly ordered final_call_args
        catchStdout,
        catchStderr,
        stdoutFile,
        stderrFile,
        isUsage
    )

    # Extract command attribute
    cmd <- attr(res_int, "command")

    # Function to collect output from file
    collect_output <- function(file_path_arg) { # Renamed parameter to avoid conflict
        if (!file.exists(file_path_arg)) {
            return(character(0))
        }

        result <- tryCatch(
            {
                readLines(file_path_arg)
            },
            error = function(e) {
                # Handle binary output by reading as raw if text reading fails
                readBin(file_path_arg, what = "raw", n = file.info(file_path_arg)$size)
            }
        )

        # Clean up temporary files, but preserve user-requested output files
        if (!is.null(saveStdout) && file_path_arg != saveStdout && file_path_arg != output) { # also check against output
            tryCatch(file.remove(file_path_arg), error = function(e) NULL)
        } else if (is.null(saveStdout) && file_path_arg != output) { # also check against output
            tryCatch(file.remove(file_path_arg), error = function(e) NULL)
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

    # Clean up the main output file if it was a temporary one and not saved via saveStdout
    if (!is.null(output) && output == stdoutFile && is.null(saveStdout)) {
        if (file.exists(output)) {
            tryCatch(file.remove(output), error = function(e) NULL)
        }
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
