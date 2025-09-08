# Tinytest for VBI index and query
library(tinytest)
library(RBCFLib)
# get script path
exdata <- system.file("exdata", package = "RBCFLib")
vcf <- file.path(exdata, "imputed.gt.vcf.gz")
# vcf <- "../1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.bcf"
# vcf <- "../data/clinvar_20250504.vcf.gz"
# vcf <- "../concat.bcf"
vbi <- paste0(vcf, ".vbi")
force_index <- FALSE
# Indexing
if (!file.exists(vbi) || force_index) {
  timing <- system.time(res_idx <- VBIIndex(vcf, vbi, Threads = 1))
  expect_true(file.exists(vbi))
  cat(sprintf("[Info] VBI indexing took %g sec\n", timing[3]))
} else {
  message("VBI index file found, skipping indexing step.")
}


# Load the VBI index once as an external pointer
cat("\n[Benchmark] Loading VBI index as external pointer...\n")
system.time(vbi_ptr <- VBIIndexLoad(vbi))
# Print the first 5 lines of the VBI index (for debug/coverage)
system.time(expect_silent(VBIPrintIndex(vbi_ptr, 5)))
VBIPrintIndex(vbi_ptr, 100) |> cat()
tim <- system.time(rmem <- VBIIndexMemoryUsage(vbi_ptr))
cat(sprintf(
  "[Info] vbi_index_t C-level memory usage: %.2f MB\n",
  rmem[["vbi_index_t_bytes"]] / 1e6
))
cat(sprintf(
  "[Info] cgranges_t  C-level memory usage: %.2f MB\n",
  rmem[["cgranges_t_bytes"]] / 1e6
))
cat(sprintf(
  "[Info] Total       C-level memory usage: %.2f MB\n",
  sum(unlist(rmem)) / 1e6
))
cat("timing for memory usage: ", tim[3])
cat("Extracting ranges from VBI index...\n")
et <- system.time({
  ranges <- VBIExtractRanges(vbi_ptr)
})
cat("Extracting ranges took: ", et[3], " sec\n")
print(length(ranges[[1]]))
# Query by index range
cat("\n[Benchmark] Querying by index range...\n")
tim <- system.time(hits2 <- VBIQueryByIndices(vcf, vbi_ptr, 5, 11))
expect_true(is.data.frame(hits2))
expect_true(all(c('chrom', 'pos', 'ref', 'alt', 'index') %in% names(hits2)))
expect_true(nrow(hits2) == 7)
cat(sprintf("Retrieved %d records in %g sec\n", nrow(hits2), tim[3]))
# second timing
system.time(hits2b <- VBIQueryByIndices(vcf, vbi_ptr, 5, 11))
stopifnot(identical(nrow(hits2), nrow(hits2b)))

cat("[Benchmark] Querying region 5x (linear scan)...\n")
nrange <- min(25, length(ranges[[1]]))
nsamples <- 2
region_str <- vector("list", nsamples)
region_str_cgranges <- vector("list", nsamples)
for (i in 1:nsamples) {
  thisSample <- sample(seq_along(ranges[[1]]), size = nrange)
  region_str[[i]] <- sample(
    paste0(
      ranges$chrom[thisSample],
      ":",
      ranges$start[thisSample],
      "-",
      ranges$end[thisSample]
    ),
    size = nrange
  ) |>
    paste0(collapse = ",")
  region_str_cgranges[[i]] <- sample(
    paste0(
      ranges$chrom[thisSample],
      ":",
      ranges$start[thisSample] - 1,
      "-",
      ranges$end[thisSample] + 1
    ),
    size = nrange
  ) |>
    paste0(collapse = ",")
}
tm1 <- system.time({
  for (i in 1:nsamples) {
    df <- VBIQueryRange(vcf, vbi_ptr, region_str[[i]])
    expect_true(is.data.frame(df))
    cat("rows=", nrow(df), "\n")
  }
})
cat("[Benchmark] Querying region 5x (cgranges)...\n")
tm2 <- system.time({
  for (i in 1:nsamples) {
    df <- VBIQueryRegionCGRanges(vcf, vbi_ptr, region_str_cgranges[[i]])
    expect_true(is.data.frame(df))
    cat("rows=", nrow(df), "\n")
  }
})
cat(sprintf("[Benchmark] Linear scan: %g sec\n", tm1[3]))
cat(sprintf("[Benchmark] cgranges:   %g sec\n", tm2[3]))
