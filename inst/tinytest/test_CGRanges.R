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
expect_true(all(hits %in% c(1, 2)))
expect_true(length(hits) == 2)

# Overlap query: should hit interval 3 only
hits2 <- CGRangesOverlap(cr, "chr1", 35, 36)
expect_true(all(hits2 == 3))

# Overlap query: should hit interval 4 only
hits3 <- CGRangesOverlap(cr, "chr2", 10, 12)
expect_true(all(hits3 == 4))

# Overlap query: should hit nothing
hits4 <- CGRangesOverlap(cr, "chr2", 20, 25)
expect_true(length(hits4) == 0)

# Clean up
CGRangesDestroy(cr)
