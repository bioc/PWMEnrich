library(PWMEnrich)
library(testthat)

# test for cases when the input is somehow invalid although it looks like a sequence

motifs.denovo = readMotifs(system.file(package="PWMEnrich", dir="extdata", file="example.transfac"), remove.acc=TRUE)

inseq = list(DNAString("AGTACCGATGACCGATGACC"), DNAString("NNNNNNNNNNNNNNNNNNN"))

# this should produce an error on the second sequence
test_that("N as input", {
	expect_error(motifScores(inseq, motifs.denovo))

	expect_error(motifEnrichment(inseq, motifs.denovo))
})
