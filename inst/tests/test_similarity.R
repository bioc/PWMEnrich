library(PWMEnrich)
library(testthat)

# test function for motif similarity
motifs = readMotifs(system.file(package="Biomotifs", dir="extdata", file="jaspar-insecta.jaspar"), remove.acc=TRUE)

# check if two things are numerically equal, e.g. to precision of 1e-8
numEqual = function(x, y, prec=1e-8){
	all(abs(x-y) < prec)
}

test_that("motifSimilarity", {
	expect_equal(motifSimilarity(motifs$tin, motifs$tin), 1)
	expect_true(numEqual(motifSimilarity(motifs$tin, motifs$vnd), 0.8785692, 1e-6))
})
