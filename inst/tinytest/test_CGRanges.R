# Tinytest for cgranges R binding
library(tinytest)
library(RBCFLib)

# Create cgranges object
cr <- CGRangesCreate()

# Add intervals
CGRangesAdd(cr, "chr1", 10, 20, 1)
CGRangesAdd(cr, "chr1", 15, 25, 2)
CGRangesAdd(cr, "chr1", 30, 40, 3)
CGRangesAdd(cr, "chr2", 5, 15, 4)
CGRangesIndex(cr)


# Overlap query: should hit intervals 1 and 2
hits <- CGRangesOverlap(cr, "chr1", 18, 22)
print(hits)
expect_true(all(hits[[1]] %in% c(1, 2)))
expect_true(length(hits[[1]]) == 2)

# Overlap query: should hit interval 3 only
hits2 <- CGRangesOverlap(cr, "chr1", 35, 36)
expect_true(all(hits2 == 3))

# Overlap query: should hit interval 4 only
hits3 <- CGRangesOverlap(cr, "chr2", 10, 12)
expect_true(all(hits3 == 4))

# Overlap query: should hit nothing
hits4 <- CGRangesOverlap(cr, "chr2", 20, 25)
print(hits4)
expect_true(length(hits4[[1]]) == 0)

# Clean up
CGRangesDestroy(cr)

# Additional tests for CGRanges and VBI R bindings
# CGRangesExtractByIndex test
cr <- CGRangesCreate()
CGRangesAdd(cr, "chr1", 10, 20, 1)
CGRangesAdd(cr, "chr1", 15, 25, 2)
CGRangesAdd(cr, "chr2", 5, 15, 3)
CGRangesIndex(cr)

# Extract by index (should match what was added)
cat('CGRangesExtractByIndex(cr, 1):\n')
print(CGRangesExtractByIndex(cr, 1))
cat('CGRangesExtractByIndex(cr, 2):\n')
print(CGRangesExtractByIndex(cr, 2))
cat('CGRangesExtractByIndex(cr, 3):\n')
print(CGRangesExtractByIndex(cr, 3))
res <- CGRangesExtractByIndex(cr, c(1, 2, 3))
print(res$chrom)
expect_equal(is.na(res$chrom), c(FALSE, FALSE, TRUE))
expect_equal(res$start, c(10, 15, 5))
expect_equal(res$end, c(20, 25, 15))
expect_equal(res$label, c(2, 3, 4)) # label is 1-based in R

# CGRangesOverlapVec test (vectorized)
hits <- CGRangesOverlapVec(cr, c("chr1", "chr2"), c(12, 6), c(18, 10))
expect_equal(length(hits), 2)
expect_true(all(sapply(hits, is.integer)))
expect_true(1 %in% hits[[1]]) # chr1:12-18 overlaps first interval
expect_true(3 %in% hits[[2]]) # chr2:6-10 overlaps third interval

CGRangesDestroy(cr)
