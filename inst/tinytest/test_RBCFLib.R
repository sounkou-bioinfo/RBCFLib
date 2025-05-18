# Main test file for RBCFLib
library(tinytest)

# Basic sanity test
expect_equal(1 + 1, 2)

# Run module-specific tests
if (requireNamespace("RBCFLib", quietly = TRUE)) {
    # Check if BCFToolsRun test file exists and source it
    bcftools_run_test <- system.file("tinytest", "testBCFToolsRun.R", package = "RBCFLib")
    if (file.exists(bcftools_run_test)) {
        source(bcftools_run_test)
    }

    # Check if BCFToolsMunge test file exists and source it
    bcftools_munge_test <- system.file("tinytest", "testBCFToolsMunge.R", package = "RBCFLib")
    if (file.exists(bcftools_munge_test)) {
        source(bcftools_munge_test)
    }
}
