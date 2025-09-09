# Tinytest for VBI index and query
library(tinytest)
library(RBCFLib)
# get script path
exdata <- system.file("exdata", package = "RBCFLib")
vcf <- file.path(exdata, "imputed.gt.vcf.gz")
# vcf <- "../1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.bcf"
# vcf <- "../clinvar_20250504.vcf.gz"
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
# First load VCF with VBI integration for the new API
vcf_ctx <- VCFLoad(vcf, vbi)
tim <- system.time(hits2 <- VBIQueryByIndices(vcf_ctx, 5, 11))
expect_true(is.data.frame(hits2))
expect_true(all(c('chrom', 'pos', 'ref', 'alt', 'index') %in% names(hits2)))
expect_true(nrow(hits2) == 7)
cat(sprintf("Retrieved %d records in %g sec\n", nrow(hits2), tim[3]))
# second timing
system.time(hits2b <- VBIQueryByIndices(vcf_ctx, 5, 11))
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
    df <- VBIQueryRange(vcf_ctx, region_str[[i]])
    expect_true(is.data.frame(df))
    cat("rows=", nrow(df), "\n")
  }
})
cat("[Benchmark] Querying region 5x (cgranges)...\n")
tm2 <- system.time({
  for (i in 1:nsamples) {
    df <- VBIQueryRegionCGRanges(vcf_ctx, region_str_cgranges[[i]])
    expect_true(is.data.frame(df))
    cat("rows=", nrow(df), "\n")
  }
})
cat(sprintf("[Benchmark] Linear scan: %g sec\n", tm1[3]))
cat(sprintf("[Benchmark] cgranges:   %g sec\n", tm2[3]))

# Test new VBI functions: VCFLoad and VBIQueryRegion
cat("\n[Test] Testing VCFLoad and VBIQueryRegion...\n")

# Load VCF with VBI integration
cat("Loading VCF with VBI integration...\n")
tm_load <- system.time(vcf_obj <- VCFLoad(vcf, vbi))
expect_true(!is.null(vcf_obj))
cat(sprintf("VCFLoad took: %g sec\n", tm_load[3]))

# Test basic VBIQueryRegion
cat("Testing basic VBIQueryRegion...\n")
tm_query1 <- system.time(
  hits_basic <- VBIQueryRegion(vcf_obj, "chr21:5030082-5030356")
)
expect_true(is.data.frame(hits_basic))
expect_true(all(
  c('chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'n_allele') %in%
    names(hits_basic)
))
cat(sprintf(
  "Basic query retrieved %d records in %g sec\n",
  nrow(hits_basic),
  tm_query1[3]
))
print(head(hits_basic))

# Test VBIQueryRegion with INFO fields
cat("Testing VBIQueryRegion with INFO fields...\n")
tm_query2 <- system.time(
  hits_with_info <- VBIQueryRegion(
    vcf_obj,
    "chr21:5030082-5030356",
    include_info = TRUE
  )
)
expect_true(is.data.frame(hits_with_info))
expect_true(all(
  c('chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'n_allele') %in%
    names(hits_with_info)
))
# Should have more columns when including INFO
expect_true(ncol(hits_with_info) >= ncol(hits_basic))
cat(sprintf(
  "INFO query retrieved %d records in %g sec\n",
  nrow(hits_with_info),
  tm_query2[3]
))
cat("Columns with INFO:", colnames(hits_with_info), "\n")
print(head(hits_with_info))

# Test with a region that should return no results
cat("Testing empty region query...\n")
empty_hits <- VBIQueryRegion(vcf_obj, "chr21:1-10")
expect_true(is.data.frame(empty_hits))
expect_true(nrow(empty_hits) == 0)
cat("Empty query returned 0 records as expected\n")

cat("[Test] VCFLoad and VBIQueryRegion tests completed successfully!\n")

# Test genotype support
cat("\n[Test] Testing VBI with genotype data...\n")

# Test VBIQueryRegion with genotypes
cat("Testing VBIQueryRegion with genotypes...\n")
tm_gt <- system.time(
  hits_with_gt <- VBIQueryRegion(
    vcf_obj,
    "chr21:5030082-5030356",
    include_genotypes = TRUE
  )
)
expect_true(is.data.frame(hits_with_gt))
expect_true("GT" %in% names(hits_with_gt))
cat(sprintf(
  "Genotype query retrieved %d records in %g sec\n",
  nrow(hits_with_gt),
  tm_gt[3]
))
cat("Columns with genotypes:", colnames(hits_with_gt), "\n")
if (nrow(hits_with_gt) > 0) {
  cat("Sample GT data (first variant):\n")
  print(head(hits_with_gt[1, "GT", drop = FALSE]))
}

# Test VBIQueryRegion with all fields
cat("Testing VBIQueryRegion with all fields...\n")
tm_all <- system.time(
  hits_all <- VBIQueryRegion(
    vcf_obj,
    "chr21:5030082-5030356",
    include_info = TRUE,
    include_format = TRUE,
    include_genotypes = TRUE
  )
)
expect_true(is.data.frame(hits_all))
expect_true(all(c("INFO", "FORMAT_IDS", "GT") %in% names(hits_all)))
cat(sprintf(
  "Full query retrieved %d records in %g sec\n",
  nrow(hits_all),
  tm_all[3]
))
cat("All columns:", colnames(hits_all), "\n")

cat("[Test] VBI genotype tests completed successfully!\n")
