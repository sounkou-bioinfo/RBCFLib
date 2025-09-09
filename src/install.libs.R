# ref : h/t @Zilong-li https://github.com/Zilong-Li/vcfppR/blob/main/src/install.libs.R
## Install shared objects and static libraries
files <- Sys.glob(paste0("*", SHLIB_EXT))

## Also copy the bcftools static library
files <- c(
  files,
  list.files("bcftools-1.22/", pattern = "\\.a$", full.names = TRUE)
)

## Create destination directory
dest <- file.path(R_PACKAGE_DIR, paste0("libs", R_ARCH))
dir.create(dest, recursive = TRUE, showWarnings = FALSE)

## Copy the shared objects and static libraries
file.copy(files, dest, overwrite = TRUE)

## Copy bcftools header files
headers <- c(
  list.files("bcftools-1.22/", pattern = "\\.h$", full.names = TRUE)
)
headers_dest <- file.path(R_PACKAGE_DIR, "include/bcftools")
dir.create(headers_dest, recursive = TRUE, showWarnings = FALSE)
file.copy(headers, headers_dest, overwrite = TRUE)

# copy htslib header giles and libraries
htslib_headers <- c(
  list.files(
    "bcftools-1.22/htslib-1.22/",
    pattern = "\\.h$",
    full.names = TRUE,
    recursive = TRUE
  )
)
htslib_headers_dest <- file.path(
  R_PACKAGE_DIR,
  "include/htslib",
  basename(dirname(htslib_headers))
)
htslib_headers_dest <- sub("htslib-1.22", "", htslib_headers_dest)
lapply(htslib_headers_dest, function(x) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
})
lapply(seq_along(htslib_headers), function(x) {
  file.copy(htslib_headers[x], htslib_headers_dest[x], overwrite = TRUE)
})

htslib_libs <- c(
  list.files(
    "htslib-1.22/",
    pattern = paste0(".*", SHLIB_EXT),
    full.names = TRUE
  )
)

htslib_libs <- c(
  htslib_libs,
  list.files("htslib-1.22/", pattern = "\\.a$", full.names = TRUE)
)

htslib_libs_dest <- file.path(R_PACKAGE_DIR, paste0("libs", R_ARCH))
dir.create(htslib_libs_dest, recursive = TRUE, showWarnings = FALSE)
file.copy(htslib_libs, htslib_libs_dest, overwrite = TRUE)

## Copy bcftools executable to inst/bin
bcftools_exe <- "bcftools-1.22/bcftools"
if (WINDOWS) {
  bcftools_exe <- paste0(bcftools_exe, ".exe")
}
if (file.exists(bcftools_exe)) {
  bin_dest <- file.path(R_PACKAGE_DIR, "bin")
  dir.create(bin_dest, recursive = TRUE, showWarnings = FALSE)
  file.copy(bcftools_exe, file.path(bin_dest, "bcftools"), overwrite = TRUE)
  Sys.chmod(file.path(bin_dest, "bcftools"), mode = "0755")
}

## Copy bcftools plugins to inst/bin/plugins (not bcftools/plugins since bcftools is a file)
plugins_dir <- "bcftools-1.22/plugins"
if (dir.exists(plugins_dir)) {
  plugins_dest <- file.path(R_PACKAGE_DIR, "bin", "plugins")
  dir.create(plugins_dest, recursive = TRUE, showWarnings = FALSE)

  # Copy all shared library files (plugins)
  plugin_files <- list.files(
    plugins_dir,
    pattern = paste0(SHLIB_EXT, "$"),
    full.names = TRUE
  )
  if (length(plugin_files) > 0) {
    # Copy each plugin file individually
    for (plugin_file in plugin_files) {
      plugin_name <- basename(plugin_file)
      dest_file <- file.path(plugins_dest, plugin_name)
      file.copy(plugin_file, dest_file, overwrite = TRUE)
      Sys.chmod(dest_file, mode = "0755")
    }
  }
}

# remove  libexec directory if it exists
libexec_dir <- file.path(R_PACKAGE_DIR, "libexec")
if (dir.exists(libexec_dir)) {
  unlink(libexec_dir, recursive = TRUE, force = TRUE)
}

## Copy SuiteSparse libraries and headers if they exist
suitesparse_src <- "SuiteSparse/install"
if (dir.exists(suitesparse_src)) {
  # Copy SuiteSparse static libraries
  suitesparse_libs <- list.files(
    file.path(suitesparse_src, "lib"),
    pattern = "\\.a$",
    full.names = TRUE
  )
  if (length(suitesparse_libs) > 0) {
    file.copy(suitesparse_libs, dest, overwrite = TRUE)
  }

  # Copy SuiteSparse headers
  suitesparse_headers_src <- file.path(suitesparse_src, "include")
  if (dir.exists(suitesparse_headers_src)) {
    suitesparse_headers_dest <- file.path(
      R_PACKAGE_DIR,
      "include",
      "suitesparse"
    )
    dir.create(suitesparse_headers_dest, recursive = TRUE, showWarnings = FALSE)

    # Copy all header files
    suitesparse_headers <- list.files(
      suitesparse_headers_src,
      pattern = "\\.h$",
      full.names = TRUE,
      recursive = TRUE
    )

    # Maintain directory structure
    for (header in suitesparse_headers) {
      rel_path <- sub(paste0(suitesparse_headers_src, "/"), "", header)
      dest_file <- file.path(suitesparse_headers_dest, rel_path)
      dest_dir <- dirname(dest_file)
      dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
      file.copy(header, dest_file, overwrite = TRUE)
    }
  }
}
