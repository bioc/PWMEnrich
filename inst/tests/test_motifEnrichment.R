library(PWMEnrich)
library(testthat)

# check if two things are numerically equal, e.g. to precision of 1e-8
numEqual = function(x, y, prec=1e-8){
	all(abs(x-y) < prec)
}

motifs.denovo = readMotifs(system.file(package="PWMEnrich", dir="extdata", file="example.transfac"), remove.acc=TRUE)

# test that motif reading in works properly
test_that("readMotifs", {
	expect_equal(length(motifs.denovo), 2)
	expect_equal(names(motifs.denovo), c("tin_like_motif", "gata_like_motif"))	
})

# create a GATA motif ACGT 
gata = rbind("A"=c(0, 10, 0, 10), "C"= c(0, 0, 0, 0), "G"=c(10, 0, 0, 0), "T"=c(0, 0, 10, 0) )
gata = t(apply(gata, 1, as.integer))

# GATA PWM and example sequence
gata.pwm = PFMtoPWM(gata)
s = DNAString("AAAAGATAAA")
s1 = DNAString("AAAAAAAAAGATA")

sequences = list(s, s1)

scores = motifScores(s, gata.pwm, raw.scores=T)
scores.count = motifScores(s, gata.pwm, cutoff=log2(exp(5)))

scores2 = motifScores(list(s,s1), gata.pwm)

test_that("motifScores works properly", {
	expect_true(numEqual( scores[[1]][5,1], 2^(gata.pwm$pwm["G",1]*4))) #  hand-calculate PWM score
	expect_equal(nrow(scores[[1]]), 20)
	
	# one significant motif hit
	expect_true( all( scores.count == matrix(1) ) )
	
	expect_equal(dim(scores2), c(2,1))
})

# now do something else
res.affinity = motifEnrichment(s, gata.pwm)
res.cutoff = motifEnrichment(s, gata.pwm, cutoff=log2(exp(5)), score="cutoff")

# check affinity and cutoff scores
test_that("motifEnrichment for raw affinity and cutoff", {
	# check if the affinity is right
	expect_true( numEqual( mean(scores[[1]], na.rm=TRUE), res.affinity$group.nobg[1] ))
	expect_equal( res.cutoff$group.nobg[1], 1)
})

###
# make backgrounds and check scores

# load the testing bg.seq object
load(system.file(package="PWMEnrich", dir="extdata", file="bg.seq-test.RData"))

## make diferent background distributions and use them to scan

bg.logn = makePWMLognBackground(bg.seq, gata.pwm, bg.len=1000)
bg.z5 = makePWMCutoffBackground(bg.seq, gata.pwm, cutoff=log2(exp(5)))

# check if the Z-score calculation works, on background the z-score should be 0 for the group
res.bg.z5 = motifEnrichment(bg.seq, bg.z5)
res.z5 = motifEnrichment(sequences, bg.z5)
test_that("motifEnrichment Z-score", {
	expect_true(numEqual(res.bg.z5$group.bg,0))
	# both sequence have one binding site
	expect_equal(as.vector(res.z5$sequence.nobg), c(1,1))
	# both have 1 binding site, but first sequence is longer
	expect_true(res.z5$sequence.bg[1,1] > res.z5$sequence.bg[2,1]) 
})

# check lognormal results 
res.logn = motifEnrichment(sequences, bg.logn)

# calculate lognormal p-value by hand
mx = bg.logn$bg.mean
sx = bg.logn$bg.sd * sqrt(1000-3) / sqrt(10-3) # first sequence is 10-4+1 = 7bp long

ml = 2*log(mx) - 0.5*log(mx^2+sx^2)
sl = sqrt(-2*log(mx) + log(mx^2+sx^2))

# compare lognormal P-value 
test_that("motifEnrichment P-value", {
	expect_true(numEqual(res.logn$sequence.bg[1,1], plnorm(res.logn$sequence.nobg[1,], meanlog=ml, sdlog=sl, lower.tail=FALSE)))
})

## other backgrounds
test_that("other background produce no error", {
	bg.p = makePWMEmpiricalBackground(bg.seq, gata.pwm)
	
	expect_equal(nrow(bg.p$bg.fwd), 10000)
	expect_equal(nrow(bg.p$bg.rev), 10000)

	bg.1e2 = makePWMPvalCutoffBackground(bg.p, 1e-2)
	expect_equal(length(bg.1e2$bg.cutoff), 1)
	expect_equal(length(bg.1e2$bg.P), 1)
	
	res.bg.1e2 = motifEnrichment(bg.seq, gata.pwm, cutoff=bg.1e2$bg.cutoff, score="cutoff")
})

bg.gev = makePWMGEVBackground(bg.seq, gata.pwm, bg.len=seq(200, 1000, 200))
test_that("GEV", {
	expect_equal(length(bg.gev$bg.loc), 1)
	expect_equal(class(bg.gev$bg.loc[[1]]), "lm")
	expect_equal(class(bg.gev$bg.scale[[1]]), "lm")		
	expect_equal(class(bg.gev$bg.shape[[1]]), "lm")		
})

####
## Differential enrichment

diff.logn = motifDiffEnrichment(s, s1, bg.logn)
diff.z5 = motifDiffEnrichment(s, s1, bg.z5)	
test_that("motifDiffEnrichment", {	
	expect_true(diff.logn$group.bg > 0)
	expect_equal(diff.z5$group.nobg, 0)
	expect_true(diff.z5$group.bg > 0)
})

### test priors calculation
prior = makePriors(list(DNAString(c("AAAACCCCCC"))), 1)
test_that("makePriors", {
	expect_equal(prior, c("A"=0.2, "C"=0.3, "G"=0.3, "T"=0.2))
})


#### Test MotifEnrichmentResults

pfm.all = list(tin_like_motif=motifs.denovo[[1]], gata_like_motif=motifs.denovo[[2]], gata=gata)
motifs.all = PFMtoPWM(pfm.all, name=c("tin", "GATA", "GATA"))

test_that("PFMtoPWM names", {
	expect_equal(sapply(motifs.all, function(x) x$name), c(tin_like_motif="tin", gata_like_motif="GATA", gata="GATA"))
})

# setup stuff for motifEnrichment

bg.logn = makePWMLognBackground(bg.seq, motifs.all)
bg.z5 = makePWMCutoffBackground(bg.seq, motifs.all, cutoff=log2(exp(4)))

sequences = list(DNAString("AAATTTAAGATAAAATTGCGT"), DNAString("AAAAAGATAAAAAAAA"))

res.logn = motifEnrichment(sequences, bg.logn)
res.z5 = motifEnrichment(sequences, bg.z5)

### test new methods for motifEnrichmentResults

test_that("MotifEnrichmentResults methods for logn", {
	expect_equal(names(motifRankingForGroup(res.logn)), c("GATA", "GATA", "tin"))
	expect_equal(names(motifRankingForGroup(res.logn, id=TRUE)), c("gata", "gata_like_motif", "tin_like_motif"))
	expect_equal(motifRankingForGroup(res.logn, rank=TRUE), c("tin"=3, "GATA"=2, "GATA"=1))
	expect_equal(motifRankingForGroup(res.logn, order=TRUE, id=TRUE), c("gata"=3, "gata_like_motif"=2, "tin_like_motif"=1))
	expect_equal(names(motifRankingForGroup(res.logn, unique=TRUE)), c("GATA", "tin"))
	expect_equal(motifRankingForGroup(res.logn, rank=TRUE, unique=TRUE), c("tin"=2, "GATA"=1))
	expect_error(motifRankingForGroup(res.logn, order=TRUE, unique=TRUE))
})

test_that("MotifEnrichmentResults methods for z5", {
	expect_equal(names(motifRankingForGroup(res.z5)), c("GATA", "GATA", "tin"))
	expect_equal(names(motifRankingForGroup(res.z5, id=TRUE)), c("gata_like_motif", "gata", "tin_like_motif"))
	expect_equal(motifRankingForGroup(res.z5, rank=TRUE), c("tin"=3, "GATA"=1, "GATA"=2))
	expect_equal(motifRankingForGroup(res.z5, order=TRUE, id=TRUE), c("gata_like_motif"=2, "gata"=3, "tin_like_motif"=1))
	expect_equal(names(motifRankingForGroup(res.z5, unique=TRUE)), c("GATA", "tin"))
	expect_equal(motifRankingForGroup(res.z5, rank=TRUE, unique=TRUE), c("tin"=2, "GATA"=1))
	expect_error(motifRankingForGroup(res.z5, order=TRUE, unique=TRUE))
})

### test the new parameter to PWMUnscaled
gatan = gata/10
gatan.pwm = PWMUnscaled(gatan, seq.count=10)
gatan.pwm2 = PFMtoPWM(gatan, seq.count=10)
test_that("seq.count parameter", {
	expect_equal(gata.pwm$pfm, gatan.pwm$pfm)
	expect_equal(gata.pwm$pwm, gatan.pwm$pwm)
	expect_equal(gata.pwm$pfm, gatan.pwm2$pfm)
	expect_equal(gata.pwm$pwm, gatan.pwm2$pwm)
})

