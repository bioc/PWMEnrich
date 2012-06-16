library(PWMEnrich)
library(testthat)
library(gtools)

# test cloverScore
scores = cbind(c(1, 5, 10, 2), c(5, 4, 3, 1))
colnames(scores) = c("m1", "m2")

clover = PWMEnrich:::cloverScore(scores, lr3=TRUE)

clover.mean = colMeans(clover)

# calculate manually using all subset combinations
s = scores[,1]
lr3 = rep(0, 4)
for(i in 1:4){
	comb = combinations(4, i)
	
	lr3[i] = mean(apply(comb, 1, function(x) prod(s[x])))
}

# clover scores should match
test_that("Clover score", {
	expect_equal(lr3, clover[,1])
})


