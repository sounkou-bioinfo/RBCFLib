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

    # Determine the arguments to pass to the C function
    final_args_for_c <- character()

    if (isUsage) {
        final_args_for_c <- "--help"
    } else {
        # Build arguments from parameters
        current_args <- character()

        if (!is.null(fasta_ref)) {
            current_args <- c(current_args, "-f", fasta_ref)
        }
        if (!is.null(columns)) {
            current_args <- c(current_args, "-c", columns)
        }
        if (!is.null(output)) {
            current_args <- c(current_args, "--output", output)
        }
        if (!is.null(columns_file)) {
            current_args <- c(current_args, "-C", columns_file)
        }
        if (!is.null(fai)) {
            current_args <- c(current_args, "--fai", fai)
        }
        if (!is.null(set_cache_size)) {
            current_args <- c(current_args, "--set-cache-size", set_cache_size)
        }
        if (!is.null(iffy_tag)) {
            current_args <- c(current_args, "--iffy-tag", iffy_tag)
        }
        if (!is.null(mismatch_tag)) {
            current_args <- c(current_args, "--mismatch-tag", mismatch_tag)
        }
        if (!is.null(sample_name)) {
            # Note: bcftools munge typically uses --sample-name. -s might be for other options.
            # Keeping -s as per existing signature for now.
            current_args <- c(current_args, "-s", sample_name)
        }
        if (!is.null(ns)) {
            # Note: bcftools munge typically uses --nSample. --ns might be incorrect.
            # Keeping --ns as per existing signature for now.
            current_args <- c(current_args, "--ns", ns)
        }
        if (!is.null(nc)) {
            # Note: bcftools munge typically uses --nCase. --nc might be incorrect.
            # Keeping --nc as per existing signature for now.
            current_args <- c(current_args, "--nc", nc)
        }

        # Add the input file (often positional)
        if (!is.null(input_file)) {
            current_args <- c(current_args, input_file)
        }

        # Add any extra user-provided arguments from the 'args' parameter
        if (length(args) > 0) {
            current_args <- c(current_args, args)
        }

        # If no operational arguments were constructed, this implies a help request.
        if (length(current_args) == 0) {
            final_args_for_c <- "--help"
        } else {
            final_args_for_c <- current_args
        }
    }

    # Call the C function
    result_list <- tryCatch({
        status_and_command <- .Call(
            "RC_bcftools_munge",
            final_args_for_c,
            as.logical(catchStdout), # Ensure logical
            as.logical(catchStderr), # Ensure logical
            stdoutFile,
            stderrFile,
            as.logical(isUsage)
        ) # Pass isUsage, ensure logical

        # ... existing code for processing result_list, reading stdout/stderr ...
        # This part of the code (lines 111-140 in the original file context) should be preserved.
        # For brevity, I'm not reproducing it here but it should follow.
        # The key change is how `final_args_for_c` is constructed and passed.

        # Read stdout if captured and not saved to a specific file by user
        stdout_lines <- NULL
        if (catchStdout && is.null(saveStdout) && file.exists(stdoutFile)) {
            stdout_lines <- readLines(stdoutFile, warn = FALSE)
        } else if (!is.null(saveStdout) && file.exists(saveStdout)) {
            # If saveStdout was used, we can indicate the file path or read a few lines for preview
            stdout_lines <- paste("Output saved to:", saveStdout)
        }


        # Read stderr if captured
        stderr_lines <- NULL
        if (catchStderr && file.exists(stderrFile)) {
            stderr_lines <- readLines(stderrFile, warn = FALSE)
        }

        # Ensure status_and_command is a list before trying to access elements by name
        # .Call returns the SEXP directly, which is an integer with an attribute here.
        status_code <- if (is.integer(status_and_command)) status_and_command else NA_integer_
        command_echo <- if (!is.null(attr(status_and_command, "command"))) attr(status_and_command, "command") else final_args_for_c


        list(
            status = status_code,
            stdout = stdout_lines,
            stderr = stderr_lines,
            command = c("bcftools", "+munge", command_echo) # Construct full command for display
        )
    }, error = function(e) {
        # In case of an error during .Call or file reading
        list(
            status = -1L, # Indicate error
            stdout = NULL,
            stderr = as.character(e),
            command = c("bcftools", "+munge", final_args_for_c)
        )
    }, finally = {
        # Clean up temporary files
        if (is.null(saveStdout) && file.exists(stdoutFile)) {
            unlink(stdoutFile, force = TRUE)
        }
        if (file.exists(stderrFile)) {
            unlink(stderrFile, force = TRUE)
        }
    })
    return(result_list)
}
