# Tests for updated BCFToolsRun function and helper functions

library(tinytest)
library(RBCFLib)

# locate example VCF
vcfFile <- system.file("exdata", "imputed.gt.vcf.gz", package = "RBCFLib")
expect_true(file.exists(vcfFile), "Sample VCF file not found")

# Test the updated BCFToolsRun
# Help command with new function signature
out <- BCFToolsRun(
  "version",
  character(),
  catchStdout = TRUE,
  catchStderr = TRUE
)
expect_equal(
  out$status,
  0L,
  info = "BCFToolsRun with updated signature works correctly"
)
expect_true(is.character(out$command), info = "Command returned")
expect_true(length(out$stdout) > 0, info = "Stdout captured")

# Test command validation
expect_error(
  BCFToolsRun("nonexistent"),
  "not a recognized bcftools command",
  info = "Invalid command validation works"
)

# Test file validation
nonexistentFile <- "/path/to/nonexistent/file.vcf"
expect_error(
  BCFToolsRun("index", nonexistentFile),
  "No such file or directory",
  info = "File validation for index command works"
)

# Test direct BCFToolsRun calls instead of helper functions
if (file.exists(vcfFile)) {
  # View header only (directly using BCFToolsRun)
  out <- BCFToolsRun("view", c("-h", vcfFile))
  expect_equal(out$status, 0L, info = "BCFToolsRun view with -h option works")
  expect_true(
    any(grepl("##fileformat=VCF", out$stdout)),
    info = "Header contains expected VCF format line"
  )

  # Get stats directly using BCFToolsRun
  out <- BCFToolsRun("stats", vcfFile)
  expect_equal(out$status, 0L, info = "BCFToolsRun stats command works")
  expect_true(
    any(grepl("summary numbers", out$stdout, ignore.case = TRUE)),
    info = "Stats contains summary information"
  )

  # Test query directly using BCFToolsRun
  out <- BCFToolsRun(
    "query",
    c("-f", "%CHROM\\t%POS\\t%REF\\t%ALT\\n", "-H", vcfFile)
  )
  expect_equal(out$status, 0L, info = "BCFToolsRun query command works")
  expect_true(
    any(grepl("CHROM", out$stdout)),
    info = "Query output contains column headers"
  )

  # We're not testing helper functions directly anymore
  # Just using BCFToolsRun directly
}

# Test output file handling
outFile <- tempfile(fileext = ".txt")
out <- BCFToolsRun("view", c("-h", vcfFile, "-o", outFile))
expect_equal(out$status, 0L, info = "Output file explicit specification works")
expect_true(file.exists(outFile), info = "Output file was created")
expect_true(file.info(outFile)$size > 0, info = "Output file contains data")

# Test isUsage parameter but also ensure we redirect to a file
outFile2 <- tempfile(fileext = ".txt")
out <- BCFToolsRun("view", c("-h", vcfFile, "-o", outFile2), isUsage = TRUE)
expect_equal(out$status, 0L, info = "isUsage parameter works")
