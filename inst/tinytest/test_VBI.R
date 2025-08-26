# Tinytest for VBI index and query
library(tinytest)
library(RBCFLib)
# get script path
exdata <- system.file("exdata", package = "RBCFLib")
vcf <- file.path(exdata, "imputed.gt.vcf.gz")
vcf <- "../1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.bcf"
#vcf <- "../data/clinvar_20250504.vcf.gz"
vbi <- paste0(vcf, ".vbi")

# Indexing
if (!file.exists(vbi)) {
  system.time(res_idx <- VBIIndex(vcf, vbi, Threads = 1))
  expect_true(file.exists(vbi))
} else {
  message("VBI index file found, skipping indexing step.")
}

# Load the VBI index once as an external pointer
cat("\n[Benchmark] Loading VBI index as external pointer...\n")
system.time(vbi_ptr <- VBIIndexLoad(vbi))
# Print the first 5 lines of the VBI index (for debug/coverage)
system.time(expect_silent(VBIPrintIndex(vbi_ptr, 5)))
VBIPrintIndex(vbi_ptr, 100) |> cat()
rmem <- VBIIndexMemoryUsage(vbi_ptr)
cat(sprintf(
  "[Info] vbi_index_t C-level memory usage: %.2f MB\n",
  rmem[["vbi_index_t_bytes"]] / 1e6
))
cat(sprintf(
  "[Info] cgranges_t  C-level memory usage: %.2f MB\n",
  rmem[["cgranges_t_bytes"]] / 1e6
))
ranges <- VBIExtractRanges(vbi_ptr)
print(length(ranges[[1]]))
# Compare outputs for the same region
region_str <- paste0(ranges$chrom, ":", ranges$start, "-", ranges$end)
region_str_cgranges <- paste0(
  ranges$chrom,
  ":",
  ranges$start - 1,
  "-",
  ranges$end + 1
)
# Query by region
print(head(region_str))
system.time(
  hits <- VBIQueryRange(vcf, vbi_ptr, region_str)
)
print(length(hits))
expect_true(is.character(hits))

# Query by index range
tim <- system.time(hits2 <- VBIQueryByIndices(vcf, vbi_ptr, 554, 10000))
expect_true(is.character(hits2))
cat(sprintf("Retrieved %d records in %g sec\n", length(hits2), tim[3]))
tim <- system.time(hits2 <- VBIQueryByIndices(vcf, vbi_ptr, 554, 10000))
cat(sprintf("Retrieved %d records in %g sec\n", length(hits2), tim[3]))


cat("[Benchmark] Querying region 100x (linear scan)...\n")
nrange <- min(100, length(ranges[[1]]))
tm1 <- system.time({
  for (i in 1:20) {
    VBIQueryRange(
      vcf,
      vbi_ptr,
      sample(region_str, size = nrange) |> paste0(collapse = ",")
    )
  }
})
cat("[Benchmark] Querying region 100x (cgranges)...\n")
tm2 <- system.time({
  for (i in 1:20) {
    VBIQueryRegionCGRanges(
      vcf,
      vbi_ptr,
      sample(region_str_cgranges, size = nrange) |> paste0(collapse = ",")
    )
  }
})
cat(sprintf("[Benchmark] Linear scan: %g sec\n", tm1[3]))
cat(sprintf("[Benchmark] cgranges:   %g sec\n", tm2[3]))
