# Tinytest for VBI index and query
library(tinytest)
library(RBCFLib)
# get script path
exdata <- system.file("exdata", package = "RBCFLib")
vcf <- file.path(exdata, "imputed.gt.vcf.gz")
vcf <- "../1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
vbi <- tempfile(fileext = ".vbi")
vbi <- "test.vbi"

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
VBIPrintIndex(vbi_ptr, 5) |> cat()

# Query by region
system.time(
    hits <- VBIQueryRange(vcf, vbi_ptr, rep("chr21:5030082-5030082", 100))
)
expect_true(is.character(hits))

# Query by index range
system.time(hits2 <- VBIQueryIndex(vcf, vbi_ptr, 100, 100000))
expect_true(is.character(hits2))
print(length(hits2))

# Benchmark: compare linear scan vs cgranges region query
region_str <- "chr21:5030082-5030082"
cat("[Benchmark] Querying region 100x (linear scan)...\n")
tm1 <- system.time({
    for (i in 1:100) {
        VBIQueryRange(vcf, vbi_ptr, region_str)
    }
})
cat("[Benchmark] Querying region 100x (cgranges)...\n")
tm2 <- system.time({
    for (i in 1:100) {
        VBIQueryRegionCGRanges(vbi_ptr, region_str)
    }
})
cat(sprintf("[Benchmark] Linear scan: %g sec\n", tm1[3]))
cat(sprintf("[Benchmark] cgranges:   %g sec\n", tm2[3]))
