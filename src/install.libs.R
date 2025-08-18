# ref : h/t @Zilong-li https://github.com/Zilong-Li/vcfppR/blob/main/src/install.libs.R
## Install shared objects and static libraries
files <- Sys.glob(paste0("*", SHLIB_EXT))

## Also copy the bcftools static library
files <- c(files, dir("bcftools-1.22/", pattern = "a$", full.names = TRUE))

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
  list.files("htslib-1.22/", pattern = "\\.h$", full.names = TRUE)
)
htslib_headers_dest <- file.path(R_PACKAGE_DIR, "include/htslib")

dir.create(htslib_headers_dest, recursive = TRUE, showWarnings = FALSE)
file.copy(htslib_headers, htslib_headers_dest, overwrite = TRUE)
htslib_libs <- c(
  list.files(
    "htslib-1.22/",
    pattern = paste0(".*", SHLIB_EXT),
    full.names = TRUE
  )
)
htslib_libs <- c(
  htslib_libs,
  dir("htslib-1.22/", pattern = "a$", full.names = TRUE)
)
htslib_libs_dest <- file.path(R_PACKAGE_DIR, paste0("libs", R_ARCH))
dir.create(htslib_libs_dest, recursive = TRUE, showWarnings = FALSE)
file.copy(htslib_libs, htslib_libs_dest, overwrite = TRUE)

## Copy bcftools executable to inst/bin
bcftools_exe <- "bcftools-1.22/bcftools"
if (file.exists(bcftools_exe)) {
  bin_dest <- file.path(R_PACKAGE_DIR, "bin")
  dir.create(bin_dest, recursive = TRUE, showWarnings = FALSE)
  file.copy(bcftools_exe, file.path(bin_dest, "bcftools"), overwrite = TRUE)
  # Make executable on Unix systems
  if (.Platform$OS.type == "unix") {
    Sys.chmod(file.path(bin_dest, "bcftools"), mode = "0755")
  }
}

## Copy bcftools plugins to inst/bin/plugins (not bcftools/plugins since bcftools is a file)
plugins_dir <- "bcftools-1.22/plugins"
if (dir.exists(plugins_dir)) {
  plugins_dest <- file.path(R_PACKAGE_DIR, "bin", "plugins")
  dir.create(plugins_dest, recursive = TRUE, showWarnings = FALSE)
  
  # Copy all .so files (plugins)
  plugin_files <- list.files(plugins_dir, pattern = "\\.so$", full.names = TRUE)
  if (length(plugin_files) > 0) {
    # Copy each plugin file individually
    for (plugin_file in plugin_files) {
      plugin_name <- basename(plugin_file)
      dest_file <- file.path(plugins_dest, plugin_name)
      file.copy(plugin_file, dest_file, overwrite = TRUE)
    }
    # Make plugins executable on Unix systems
    if (.Platform$OS.type == "unix") {
      for (plugin in list.files(plugins_dest, pattern = "\\.so$", full.names = TRUE)) {
        Sys.chmod(plugin, mode = "0755")
      }
    }
  }
}

# remove  libexec directory if it exists
libexec_dir <- file.path(R_PACKAGE_DIR, "libexec")
if (dir.exists(libexec_dir)) {
  unlink(libexec_dir, recursive = TRUE, force = TRUE)
}
