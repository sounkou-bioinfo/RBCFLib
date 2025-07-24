# Test for BCFToolsPipe functionality
library(tinytest)
library(RBCFLib)

# Test file path
vcf_file <- system.file("exdata", "imputed.gt.vcf.gz", package = "RBCFLib")
expect_true(file.exists(vcf_file), info = "Test VCF file exists")

# Test 1: Basic piping functionality - view and filter
result <- BCFToolsPipe("view", c(vcf_file), "view", c("-H"))

# Check status codes - both should be 0 for success
expect_identical(result$status[1], 0L, info = "First command (view) executed successfully")
expect_identical(result$status[2], 0L, info = "Second command (view -H) executed successfully")

# Check that output is captured correctly
expect_true(is.character(result$stdout), info = "Stdout captured as character vector")

# Check command attribute
expect_true(any(grepl("bcftools.*view.*\\|.*view", paste(result$command, collapse=" "))), 
            info = "Command string shows pipe operation")

# Test 2: Pipe to a file
output_file <- tempfile(fileext = ".vcf")
on.exit(unlink(output_file), add = TRUE)

result2 <- BCFToolsPipe("view", c(vcf_file),
                       "view", c("-h"),
                       saveStdout = output_file)

# Make sure the command ran without errors
expect_identical(result2$status[1], 0L, info = "First command (view) executed successfully")
expect_identical(result2$status[2], 0L, info = "Second command (view) executed successfully")

# Check that the file was created
expect_true(file.exists(output_file), info = "Output file was created")

# Check file content
file_lines <- readLines(output_file)
expect_true(length(file_lines) > 0, info = "Output file contains content")
expect_true(any(grepl("^#", file_lines)), info = "Output file contains VCF header lines")

# Test 3: Error handling - invalid command
expect_error(BCFToolsPipe("invalid_cmd", character(0), "view", c("-h")),
             pattern = "not a recognized bcftools command",
             info = "Error on invalid first command")

expect_error(BCFToolsPipe("view", c("-h"), "invalid_cmd", character(0)),
             pattern = "not a recognized bcftools command",
             info = "Error on invalid second command")

# Test 4: More practical piping example - view to BCF then convert back to VCF
result3 <- BCFToolsPipe("view", c("-Ob", vcf_file),
                       "view", c("-Ov"))

expect_identical(result3$status[1], 0L, info = "First view command executed successfully")
expect_identical(result3$status[2], 0L, info = "Second view command executed successfully")
expect_true(length(result3$stdout) > 0, info = "Output contains content")
expect_true(any(grepl("^##", result3$stdout)), info = "Output contains VCF format")
