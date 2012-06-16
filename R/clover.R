# Calculate the clover score

#' Calculate the Clover score using the recursive formula from Frith et al
#'
#' @param scores a matrix of average odds scores, where columns are motifs, and rows sequences
#' @param lr3 if to return a matrix of LR3 scores, where columns correpond to motifs, and rows to subset sizes
#' @param verbose if to produce verbose output of progress
#' @return the LR4 score, which is the mean of LR3 scores over subset sizes
cloverScore = function(scores, lr3=FALSE, verbose=FALSE){
	if(is.vector(scores))
		return(scores)

	# LR3 scores
	clover = matrix(0, nrow=nrow(scores), ncol=ncol(scores))
	colnames(clover) = colnames(scores)
	rownames(clover) = rownames(scores)
	
	# iterate over motifs
	for(m.inx in 1:ncol(scores)){
		if(verbose)
			message(paste("Calculating Clover score for motif", m.inx, "/", ncol(scores)))
		s = scores[,m.inx]
		
		N = length(s)
		A = matrix(NA, ncol=N+1, nrow=N+1)
		
		# load initial values
		A[1,] = 1
		for(i in 2:ncol(A)){
			A[i,i-1] = 0
		}
		
		# calculate the A matrix from the paper (Frith et al, NAR, 2004)
		for(i in 2:(N+1)){
			for(j in 2:(N+1)){
				if(j >= i){
					# i,j in their formula
					di = i-1
					dj = j-1
					# i,j is still used to subscript A
					A[i,j] = ( di * s[dj] * A[i-1,j-1] + (dj - di) * A[i,j-1] ) / dj	
				}
			}
		}
		clover[,m.inx] = A[2:(N+1),N+1]
		 
	}
	
	# final score (LR4) is the mean of LR3 scores
	if(lr3)
		clover
	else
		colMeans(clover)
}

#' Calculate the Clover P-value as described in the Clover paper
#'
#' This function only take one background sequence as input, it also just calculates the P-value
#' so it is more efficient. 
#'
#' @param scores the affinity scores for individual sequences
#' @param seq.len lengths of sequences
#' @param pwm.len lengths of PWMs
#' @param bg.fwd the raw score of forward strand
#' @param bg.rev the raw scores of reverse strand
#' @param B the number of random replicates
#' @param verbose if to give verbose progress reports
#' @param clover the clover scores if already calculated
#' @return P-value
cloverPvalue1seq = function(scores, seq.len, pwm.len, bg.fwd, bg.rev, B=1000, verbose=TRUE, clover=NULL){

	if(is.vector(scores))
		scores = matrix(scores, nrow=1, dimnames=list(NULL, names(scores)))

	# length of sequences
	bg.len = nrow(bg.fwd)
	
	if(max(seq.len) > bg.len){
		stop("The maximal string length in 'sequences' is greater than the background length")
	}
	
	# original scores for the set of sequences
	if(is.null(clover))
		original.clover = cloverScore(scores)
	else
		original.clover = clover
	
	# the empricial p-value
	clover.p = rep(0, ncol(scores))
	names(clover.p) = colnames(scores)
	
	# do average over strands
	bg.res = (bg.fwd + bg.rev) / 2
		
	for(i in 1:B){
		if(verbose)
			message("Random sample ", i, " / ", B)
		
		# sample background sequences with the same length distribution
		new.scores = matrix(0, nrow=nrow(scores), ncol=ncol(scores))
		colnames(new.scores) = colnames(scores)
		
		for(j in 1:nrow(scores)){
			# choose a subsequence with matching length
			start.range = 1:(bg.len - seq.len[j])
			start = sample(start.range, 1)
			
			# fetch the scores
			new.scores[j,] = colMeans(bg.res[start:(start+seq.len[j]-1),], na.rm=TRUE)
			
		}

		# calculate scores and clover
		new.clover = cloverScore(new.scores)
		
		clover.p = clover.p + ((new.clover > original.clover) + 0)
		
	}
	
	return(clover.p / B)
}


