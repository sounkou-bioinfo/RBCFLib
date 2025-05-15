# Tests for IndexVcf using tinytest

library(tinytest)

# Test that the function exists
expect_true(exists("indexVcf"))

# Locate the example VCF file (installed or development)
vcf_file <- system.file("exdata", "imputed.gt.vcf.gz", package = "RBCFLib")
if (vcf_file == "" || !file.exists(vcf_file)) {
    # fallback to development source file
    vcf_file <- file.path("inst", "exdata", "imputed.gt.vcf.gz")
}
# Skip if test file not found
exit_if_not(file.exists(vcf_file))

# Prepare temporary file for indexing
temp_file <- tempfile(fileext = ".vcf.gz")
file.copy(vcf_file, temp_file)

# Test creating a CSI index without error
result <- tryCatch(
    indexVcf(temp_file, index_type = "csi", force = TRUE),
    error = function(e) e
)
expect_true(!inherits(result, "error"))
# Check if index was created
expect_true(file.exists(paste0(temp_file, ".csi")))

# Clean up
if (file.exists(temp_file)) file.remove(temp_file)
if (file.exists(paste0(temp_file, ".csi"))) file.remove(paste0(temp_file, ".csi"))

# Test input validation
expect_error(indexVcf("nonexistent_file.vcf.gz"), "File does not exist")
expect_error(indexVcf(vcf_file, index_type = "invalid"), "Invalid index type")
