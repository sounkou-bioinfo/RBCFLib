# Tinytest for VBI index and query
library(tinytest)
library(RBCFLib)
# get script path
exdata <- system.file("exdata", package = "RBCFLib")
vcf <- file.path(exdata, "imputed.gt.vcf.gz")
vcf <- "../1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
vbi <- tempfile(fileext = ".vbi")


# Indexing
system.time(res_idx <- VBIIndex(vcf, vbi, Threads = 2))
expect_true(file.exists(vbi))

# Print the first 5 lines of the VBI index (for debug/coverage)
system.time(expect_silent(VBIPrintIndex(vbi, 5)))
VBIPrintIndex(vbi, 5) |> cat()
# Query by region
system.time(hits <- VBIQueryRange(vcf, vbi, "chr21:5030082-5030082"))
expect_true(is.character(hits))

# Query by index range
system.time(hits2 <- VBIQueryIndex(vcf, vbi, 1, 2))
expect_true(is.character(hits2))
