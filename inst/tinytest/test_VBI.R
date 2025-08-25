# Tinytest for VBI index and query
library(tinytest)
library(RBCFLib)
exdata <- system.file("exdata", package = "RBCFLib")
vcf <- file.path(exdata, "imputed.gt.vcf.gz")
vbi <- tempfile(fileext = ".vbi")

# Indexing
res_idx <- VBIIndex(vcf, vbi)
expect_true(file.exists(vbi))

# Query by region
hits <- VBIQueryRange(vcf, vbi, "chr21:5030082-5030082")
print(hits)
expect_true(is.character(hits))

# Query by index range
hits2 <- VBIQueryIndex(vcf, vbi, 1, 2)
expect_true(is.character(hits2))
cat("Hits2:", hits2, "\n")
