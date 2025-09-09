# Helper script to extract clean compiler paths for CMake
# Based on the pattern from configure.R examples

get_clean_compiler_paths <- function() {
    # Get R's compiler configuration
    CC_RAW <- system2(
        "${R_HOME}/bin/R",
        c("CMD", "config", "CC"),
        stdout = TRUE
    )
    CXX_RAW <- system2(
        "${R_HOME}/bin/R",
        c("CMD", "config", "CXX"),
        stdout = TRUE
    )

    # Parse CC configuration
    CC_ARGS <- strsplit(CC_RAW, " ")[[1]]

    # Handle ccache if present
    uses_ccache <- FALSE
    if (grepl("ccache", CC_ARGS[1])) {
        uses_ccache <- TRUE
        CC <- paste(CC_ARGS[-1], collapse = " ")
    } else {
        CC <- CC_ARGS[1]
    }

    # Extract just the compiler executable (first part before any flags)
    CC_COMPILER <- strsplit(CC, " ")[[1]][1]

    # Get full path
    CC_FULL <- normalizePath(Sys.which(CC_COMPILER), winslash = "/")

    # Handle CXX similarly if needed
    CXX_ARGS <- strsplit(CXX_RAW, " ")[[1]]
    if (grepl("ccache", CXX_ARGS[1])) {
        CXX <- paste(CXX_ARGS[-1], collapse = " ")
    } else {
        CXX <- CXX_ARGS[1]
    }
    CXX_COMPILER <- strsplit(CXX, " ")[[1]][1]
    CXX_FULL <- normalizePath(Sys.which(CXX_COMPILER), winslash = "/")

    return(list(
        CC_FULL = CC_FULL,
        CXX_FULL = CXX_FULL,
        uses_ccache = uses_ccache
    ))
}

# Export for use in shell scripts
if (length(commandArgs(trailingOnly = TRUE)) > 0) {
    paths <- get_clean_compiler_paths()
    cat(sprintf("CC_FULL=%s\n", paths$CC_FULL))
    cat(sprintf("CXX_FULL=%s\n", paths$CXX_FULL))
}
