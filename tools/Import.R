#!/usr/bin/Rscript
# ref : https://github.com/pysam-developers/pysam/tree/master/import
# Import.R - Process bcftools C source files for R package integration
#
# This script:
# 1. Locates all C files in the bcftools source directory
# 2. Creates modified versions with RBCFLIB extensions
# 3. Adds bcftools.RBCFLIB.h include and modifies I/O operations
#

# Get the directory of the currently running script
get_script_path <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    script_path <- NULL
    for (arg in args) {
        if (grepl("--file=", arg)) {
            script_path <- sub("--file=", "", arg)
            break
        }
    }

    if (is.null(script_path)) {
        # Fallback to current working directory if script path cannot be determined
        script_path <- getwd()
    } else {
        script_path <- dirname(normalizePath(script_path))
    }

    return(script_path)
}

Import <- function(bcftools_src_dir, dest_dir = NULL) {
    # Default to the parent directory if dest_dir not provided
    if (is.null(dest_dir)) {
        dest_dir <- dirname(bcftools_src_dir)
    }

    # Ensure directories exist
    stopifnot(dir.exists(bcftools_src_dir))
    if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE)

    # Define files to exclude
    exclude_patterns <- c("test/", "version.c", "plugins/", "RBCFLIB\\.c$")

    # Find all .c files in source directory
    c_files <- list.files(bcftools_src_dir,
        pattern = "\\.c$",
        recursive = TRUE,
        full.names = TRUE
    )

    # Filter out excluded files
    for (pattern in exclude_patterns) {
        c_files <- c_files[!grepl(pattern, c_files)]
    }

    message(sprintf("Found %d C files to process", length(c_files)))

    # Keep track of the RBCFLIB files we create
    rbcflib_files <- character()

    # Process each file
    modified_count <- 0
    for (src_file in c_files) {
        # Determine output file path
        rel_path <- sub(paste0("^", gsub("([.\\\\+])", "\\\\\\1", bcftools_src_dir)), "", src_file)
        rel_dir <- dirname(rel_path)

        # Create destination directory if it doesn't exist
        dest_subdir <- file.path(dest_dir, rel_dir)
        if (!dir.exists(dest_subdir)) dir.create(dest_subdir, recursive = TRUE)

        # Set destination filename
        basename_no_ext <- tools::file_path_sans_ext(basename(src_file))
        dest_file <- file.path(dest_subdir, paste0(basename_no_ext, ".RBCFLIB.c"))

        # Add to our list of RBCFLIB files
        rbcflib_files <- c(rbcflib_files, paste0(basename_no_ext, ".RBCFLIB.o"))

        # Read source file
        lines <- readLines(src_file, warn = FALSE, encoding = "UTF-8")

        # Modify the content
        modified_content <- ProcessCFile(lines, basename_no_ext)

        # Write to destination
        writeLines(modified_content, dest_file, useBytes = TRUE)
        modified_count <- modified_count + 1
    }

    message(sprintf("Processed %d files", modified_count))

    # Get path to script directory where our header files are located
    script_dir <- get_script_path()

    # Copy header and implementation files from script directory to dest_dir
    file.copy(file.path(script_dir, "bcftools.RBCFLIB.h"),
        file.path(dest_dir, "bcftools.RBCFLIB.h"),
        overwrite = TRUE
    )

    file.copy(file.path(script_dir, "bcftools.RBCFLIB.c"),
        file.path(dest_dir, "bcftools.RBCFLIB.c"),
        overwrite = TRUE
    )

    # Add bcftools.RBCFLIB.o to our list
    rbcflib_files <- c("main.RBCFLIB.o", "bcftools.RBCFLIB.o", rbcflib_files)

    # Remove "main.RBCFLIB.o" from its original position if it's duplicated
    rbcflib_files <- unique(rbcflib_files)

    # Process the Makefile
    makefile_path <- file.path(bcftools_src_dir, "Makefile")
    if (file.exists(makefile_path)) {
        dest_makefile <- file.path(dest_dir, "Makefile.RBCFLIB")
        ProcessMakefile(makefile_path, dest_makefile, rbcflib_files)
        message(sprintf("Created modified Makefile at %s", dest_makefile))
    } else {
        warning("Could not find Makefile in bcftools source directory")
    }

    return(invisible(modified_count))
}

# Process a C file to make it compatible with R
ProcessCFile <- function(lines, basename) {
    # Add include for our header at the top
    # Find first include line
    first_include <- grep("^\\s*#\\s*include", lines)[1]

    if (!is.na(first_include)) {
        lines <- c(
            lines[1:(first_include - 1)],
            "#include \"bcftools.RBCFLIB.h\"",
            lines[first_include:length(lines)]
        )
    } else {
        # If no include found, add it after any initial comments
        lines <- c("#include \"bcftools.RBCFLIB.h\"", lines)
    }

    # Replace main function if present
    main_pattern <- "\\bint\\s+main\\s*\\("
    main_lines <- grep(main_pattern, lines)

    if (length(main_lines) > 0) {
        for (line_idx in main_lines) {
            lines[line_idx] <- gsub(main_pattern, "int bcftools_main(", lines[line_idx])
        }
    }

    # Replace standard I/O operations
    io_replacements <- list(
        "\\bstderr\\b" = "bcftools_stderr",
        "\\bstdout\\b" = "bcftools_stdout",
        "\\bexit\\s*\\(" = "bcftools_exit(",
        "\\bprintf\\s*\\(" = "fprintf(bcftools_stdout, ",
        "\\bputs\\s*\\(" = "bcftools_puts(",
        "putchar\\s*\\(([^)]+)\\)" = "fputc(\\1, bcftools_stdout)"
    )

    for (pattern in names(io_replacements)) {
        for (i in 1:length(lines)) {
            # Skip lines with preprocessor directives to avoid changing macros
            if (!grepl("^\\s*#", lines[i])) {
                lines[i] <- gsub(pattern, io_replacements[[pattern]], lines[i])
            }
        }
    }

    return(lines)
}

# Process the Makefile to adapt it for RBCFLIB
ProcessMakefile <- function(src_makefile, dest_makefile, rbcflib_files = NULL) {
    # Read the original Makefile
    makefile_content <- readLines(src_makefile, warn = FALSE, encoding = "UTF-8")

    # Define the RCOBJECTS value
    if (!is.null(rbcflib_files) && length(rbcflib_files) > 0) {
        # Use the actual list of RBCFLIB files we generated
        rcobjects_value <- paste(rbcflib_files, collapse = " ")
    } else {
        # Fallback to the pattern approach if no list is provided
        rcobjects_value <- "main.RBCFLIB.o bcftools.RBCFLIB.o $(patsubst %.o,%.RBCFLIB.o,$(filter-out main.o version.o, $(OBJS)))"
    }

    # Replace @RCOBJECTS@ with the actual value
    for (i in seq_along(makefile_content)) {
        makefile_content[i] <- gsub("@RCOBJECTS@", rcobjects_value, makefile_content[i])
    }

    # Write the modified Makefile
    writeLines(makefile_content, dest_makefile, useBytes = TRUE)

    return(TRUE)
}

# If run from command line
if (!interactive()) {
    args <- commandArgs(TRUE)
    if (length(args) < 1) {
        cat("Usage: Rscript Import.R path/to/bcftools-version [destination-dir]\n")
        quit(status = 1)
    }

    src_dir <- args[1]
    dest_dir <- if (length(args) > 1) args[2] else NULL

    Import(src_dir, dest_dir)
}
