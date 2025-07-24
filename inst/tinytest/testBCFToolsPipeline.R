# Test for BCFToolsPipeline functionality
library(tinytest)
library(RBCFLib)

# Test file path
vcf_file <- system.file("exdata", "imputed.gt.vcf.gz", package = "RBCFLib")
expect_true(file.exists(vcf_file), info = "Test VCF file exists")

# Test 1: Basic pipeline functionality - view and filter
result <- BCFToolsPipeline("view", c(vcf_file), "view", c("-H"))

# Check status codes - both should be 0 for success
expect_identical(result$status[1], 0L, info = "First command (view) executed successfully")
expect_identical(result$status[2], 0L, info = "Second command (view -H) executed successfully")

# Check that output is captured correctly
expect_true(is.character(result$stdout), info = "Stdout captured as character vector")

# Check command attribute
expect_true(any(grepl("bcftools.*view.*\\|.*view", paste(result$command, collapse = " "))),
    info = "Command string shows pipe operation"
)

# Test 2: Simple two-command pipeline
result2 <- BCFToolsPipeline(
    "view", c("-h", vcf_file),
    "view", c("-H")
)

# Make sure the command ran without errors
expect_identical(result2$status[1], 0L, info = "First command (view -h) executed successfully")
expect_identical(result2$status[2], 0L, info = "Second command (view -H) executed successfully")

# Test 3: Error handling - invalid command
expect_error(BCFToolsPipeline("invalid_cmd", character(0), "view", c("-h")),
    pattern = "not a recognized bcftools command",
    info = "Error on invalid first command"
)

expect_error(BCFToolsPipeline("view", c("-h"), "invalid_cmd", character(0)),
    pattern = "not a recognized bcftools command",
    info = "Error on invalid second command"
)

# Test 4: Test robust -o/--output validation - only last command can have -o
expect_error(
    BCFToolsPipeline(
        "view", c("-o", "temp.vcf", vcf_file),
        "view", c("-H")
    ),
    pattern = "contains -o.*option.*only the last command",
    info = "Error when -o used in non-final command"
)

# Test 5: Error handling for -o with excluded commands
expect_error(
    BCFToolsPipeline(
        "view", c(vcf_file),
        "stats", c("-o", "stats.txt")
    ),
    pattern = "does not support -o.*option",
    info = "Error when -o used with excluded command (stats)"
)

# Test 6: Error handling for --output in non-final command
expect_error(
    BCFToolsPipeline(
        "view", c("--output", "temp.vcf", vcf_file),
        "view", c("-H")
    ),
    pattern = "contains.*--output.*option.*only the last command",
    info = "Error when --output used in non-final command"
)
