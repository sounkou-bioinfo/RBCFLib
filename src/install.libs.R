# ref : h/t @Zilong-li https://github.com/Zilong-Li/vcfppR/blob/main/src/install.libs.R
## Install shared objects and static libraries
files <- Sys.glob(paste0("*", SHLIB_EXT))

## Also copy the bcftools static library
files <- c(files, dir("bcftools-1.21/", pattern = "a$", full.names = TRUE))

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
    list.files("htslib-1.22/", pattern = paste0(".*", SHLIB_EXT), full.names = TRUE)
)
htslib_libs <- c(htslib_libs, dir("htslib-1.22/", pattern = "a$", full.names = TRUE))
htslib_libs_dest <- file.path(R_PACKAGE_DIR, paste0("libs", R_ARCH))
dir.create(htslib_libs_dest, recursive = TRUE, showWarnings = FALSE)
file.copy(htslib_libs, htslib_libs_dest, overwrite = TRUE)
