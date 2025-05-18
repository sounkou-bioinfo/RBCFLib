#
# Tests for BCFToolsMunge function
#

library(tinytest)
library(RBCFLib)

# Setup test files
inputFile <- system.file("exdata", "test_plink.tsv", package = "RBCFLib")
fastaRef <- system.file("exdata", "Test.fa", package = "RBCFLib")
colHeaders <- system.file("exdata", "colheaders.tsv", package = "RBCFLib")
outputFile <- tempfile(fileext = ".bcf")
outputFileVCF <- tempfile(fileext = ".vcf")

# Check if input files exist
expect_true(file.exists(inputFile), "Test input file doesn't exist")
expect_true(file.exists(fastaRef), "Test reference FASTA doesn't exist")
expect_true(file.exists(colHeaders), "Test column headers file doesn't exist")

# Test 1: Basic functionality with PLINK preset
test_basic <- BCFToolsMunge(
    InputFileName = inputFile,
    Columns = "PLINK",
    FastaRef = fastaRef,
    OutputFile = outputFile,
    OutputType = "b"
)

# Check return structure
expect_true(is.list(test_basic), "BCFToolsMunge should return a list")
expect_true(
    all(c("status", "stdout", "stderr", "command") %in% names(test_basic)),
    "BCFToolsMunge return list missing expected elements"
)
expect_identical(as.integer(test_basic$status), 0L, "BCFToolsMunge should exit with status 0")

# Check that output file was created
expect_true(file.exists(outputFile), "Output BCF file was not created")

# Test 1.5: Add debug print to help diagnose the issue
cat("Before Test 2, checking fastaRef path: ", fastaRef, "\n")
if (file.exists(fastaRef)) {
    cat("Reference FASTA exists\n")
} else {
    cat("Reference FASTA does not exist\n")
}
if (file.exists(paste0(fastaRef, ".fai"))) {
    cat("Reference FASTA index exists\n")
} else {
    cat("Reference FASTA index does not exist\n")
}

# Test 2: Test with custom column headers
test_custom <- BCFToolsMunge(
    InputFileName = inputFile,
    ColumnsFile = colHeaders,
    FastaRef = fastaRef,
    FaiFile = paste0(fastaRef, ".fai"),
    SampleName = "TESTSAMPLE",
    OutputFile = outputFileVCF,
    OutputType = "v"
)

expect_identical(as.integer(test_custom$status), 0L, "BCFToolsMunge with custom columns should exit with status 0")
expect_true(file.exists(outputFileVCF), "Output VCF file was not created")

# Test 3: Test error conditions
# Missing required arguments
expect_error(
    BCFToolsMunge(),
    "InputFileName is required"
)

expect_error(
    BCFToolsMunge(InputFileName = inputFile),
    "Either Columns or ColumnsFile must be provided"
)

expect_error(
    BCFToolsMunge(InputFileName = inputFile, Columns = "PLINK", ColumnsFile = colHeaders),
    "Only one of Columns or ColumnsFile should be provided, not both"
)

expect_error(
    BCFToolsMunge(InputFileName = inputFile, Columns = "PLINK"),
    "Either FastaRef or FaiFile must be provided"
)

# Test 4: Test with various options
test_options <- BCFToolsMunge(
    InputFileName = inputFile,
    Columns = "PLINK",
    FastaRef = fastaRef,
    SampleName = "TEST_GWAS",
    NumSamples = 1000,
    NumCases = 500,
    EffSampleSize = 950,
    NoVersion = TRUE,
    OutputFile = tempfile(fileext = ".vcf.gz"),
    OutputType = "z",
    WriteIndex = TRUE
)

expect_identical(as.integer(test_options$status), 0L, "BCFToolsMunge with extra options should exit with status 0")

# Test 5: Test with non-existent input file
nonexistentFile <- tempfile(fileext = ".tsv")
expect_error(
    BCFToolsMunge(
        InputFileName = nonexistentFile,
        Columns = "PLINK",
        FastaRef = fastaRef
    ),
    "Could not open .* No such file or directory"
)

# Clean up temp files
if (file.exists(outputFile)) file.remove(outputFile)
if (file.exists(outputFileVCF)) file.remove(outputFileVCF)

# Report successful completion
cat("All BCFToolsMunge tests completed\n")
