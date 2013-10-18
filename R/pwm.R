
# Globals from Biostrings
DNA_BASES = c("A", "C", "G", "T")
DNA_ALPHABET = c("A", "C", "G", "T", "M", "R", "W", "S", "Y", "K", "V", "H", "D", "B", "N", "-", "+")

#' Input parameter normalization function for PWMUnscaled
#'
#' This function is from Biostrings package
#'
#' @param prior.params Typical 'prior.params' vector: c(A=0.25, C=0.25, G=0.25, T=0.25)
.normargPriorParams <- function(prior.params)
{
    if (!is.numeric(prior.params))
        stop("'prior.params' must be a numeric vector")
    if (length(prior.params) != length(DNA_BASES) ||
        !setequal(names(prior.params), DNA_BASES))
        stop("'prior.params' elements must be named A, C, G and T")
    ## Re-order the elements.
    prior.params <- prior.params[DNA_BASES]
    if (any(is.na(prior.params)) || any(prior.params < 0))
        stop("'prior.params' contains NAs and/or negative values")
    prior.params
}

#' Input parameter normalization for PWMUnscaled
#'
#' This function is from Biostrings package.
#' A Position Frequency Matrix (PFM) is also represented as an ordinary
#' matrix. Unlike a PWM, it must be of type integer (it will typically be
#' the result of consensusMatrix()).
#'
#' @param x a frequency matrix
.normargPfm <- function(x)
{
    if (!is.matrix(x) || !is.integer(x))
        stop("invalid PFM 'x': not an integer matrix")
    ## Check the row names.
    if (is.null(rownames(x)))
        stop("invalid PFM 'x': no row names")
    if (!all(rownames(x) %in% DNA_ALPHABET))
        stop("invalid PFM 'x': row names must be in 'DNA_ALPHABET'")
    if (!all(DNA_BASES %in% rownames(x)))
        stop("invalid PFM 'x': row names must contain A, C, G and T")
    if (any(duplicated(rownames(x))))
        stop("invalid PFM 'x': duplicated row names")
    ## Check the nb of cols.
    if (ncol(x) == 0L)
        stop("invalid PFM 'x': no columns")
    ## Check the values.
    if (any(is.na(x)) || any(x < 0L))
        stop("invalid PFM 'x': values cannot be NA or negative")
    if (any(x[!(rownames(x) %in% DNA_BASES), ] != 0L))
        stop("invalid PFM 'x': IUPAC ambiguity letters are represented")
    x <- x[DNA_BASES, , drop=FALSE]
    x
}

#' The PWM function from Biostrings without unit scaling
#'
#' By default the Biostrings package scales the log-odds score so it is within 0 and 1. In this function
#' we take a more traditional approach with no unit scaling and offer unit scaling as an additional parameter.
#'
#' See ?PWM from Biostrings for more information on input arguments. 
#'
#' @title Create a PWM from PFM
#' @param x the integer count matrix representing the motif, rows as nucleotides
#' @param id a systematic ID given to this PWM, could include the source, version, etc
#' @param name the name of the transcription factor (TF) to which the PWM corresponds to
#' @param type the type of PWM calculation, either as log2-odds, or posterior probability (frequency matrix)
#' @param prior.params the pseudocounts for each of the nucleotides 
#' @param pseudo.count the pseudo-count values if different from priors
#' @param unit.scale if to unit.scale the pwm (default is no unit scaling)
#' @param seq.count if x is a normalised PFM (i.e. with probabilities instead of sequence counts), then this sequence count
#'                  will be used to convert \code{x} into a count matrix
#' @return a new PWM object representing the PWM 
#' @export
#' @examples
#'
#' if(require("PWMEnrich.Dmelanogaster.background")){
#'    data(MotifDb.Dmel.PFM)
#'
#'    PWMUnscaled(MotifDb.Dmel.PFM$ttk, id="ttk-JASPAR", name="ttk") # make a PWM with uniform background
#'    PWMUnscaled(MotifDb.Dmel.PFM$ttk, id="ttk-JASPAR", name="ttk", prior.params=c("A"=0.2, "C"=0.3, "G"=0.3, "T"=0.2)) # custom background
#'
#'    prior = getBackgroundFrequencies("dm3", quick=TRUE) # get background for drosophila (quick mode on a reduced dataset)
#'    PWMUnscaled(MotifDb.Dmel.PFM$ttk, id="ttk-JASPAR", name="ttk", prior.params=prior) # convert using genomic background
#' }
#'
PWMUnscaled = function(x, id="", name="", type=c("log2probratio", "prob"), prior.params=c(A=0.25, C=0.25, G=0.25, T=0.25), pseudo.count=prior.params, 
	unit.scale=FALSE, seq.count=NULL){
	
	# convert to PFM if needed
	if(!is.null(seq.count)){
		x = apply(round(x * seq.count), 1:2, as.integer)
	}
	
	# match input params
    x <- .normargPfm(x)
    nseq <- colSums(x)
    type <- match.arg(type)
    prior.params <- .normargPriorParams(prior.params)
    
    # total prior added
    priorN <- sum(prior.params)
    
    # posterior probabilities with unequal number of counts per row
    postProbs <- divideRows(x + prior.params, nseq + priorN)
    rownames(postProbs) = rownames(x)
    colnames(postProbs) = colnames(x)
    
    # calculate logodds, or return probabilities
    if (type == "log2probratio") {
        if (any(prior.params == 0)) 
            stop("infinite values in PWM due to 0's in 'prior.params'")
        prior.probs <- prior.params/priorN
        ans <- log2(postProbs/prior.probs)
    }
    else {
        ans <- postProbs
    }
    
    # if to apply unit.scale
    if(unit.scale)
	   	ans = unitScale(ans)
	   	
	# return an object 
	new("PWM", id=id, name=name, pfm=x, prior.params=prior.params, pwm=ans)
}

#' Scan the whole sequence on both strands
#'
#' The whole sequence is scanned with a PWM and scores returned beginning at each position. Partial motif
#' matches are not done, thus the last #[length of motif]-1 scores are NA. 
#' 
#' The function returns either an odds average (*not* log-odds average), maximal score on each strand,
#' or scores on both strands. 
#'
#' The function by default returns the score in log2 following the package \code{Biostrings}. 
#'
#' @param pwm PWM object
#' @param dna a DNAString or other sequence from Biostrings
#' @param pwm.rev the reverse complement for a pwm (if it is already pre-computed)
#' @param odds.score if to return raw scores in odds (not logodds) space
#' @param both.strands if to return results on both strands
#' @param strand.fun which function to use to summarise values over two strands (default is "mean")
#' @return a vector representing scores starting at each position, or a matrix with score in the two strands
#' @export
#' @examples
#' 
#' if(require("PWMEnrich.Dmelanogaster.background")){
#'    data(MotifDb.Dmel)
#'
#'   scanWithPWM(MotifDb.Dmel$ttk, DNAString("CGTAGGATAAAGTAACT")) # odds average over the two strands expressed as log2-odds
#'   scanWithPWM(MotifDb.Dmel$ttk, DNAString("CGTAGGATAAAGTAACT"), both.strands=TRUE) # log2-odds scores on both strands
#' }
#' 
scanWithPWM = function(pwm, dna, pwm.rev=NULL, odds.score=FALSE, both.strands=FALSE, strand.fun="mean"){
	if(is.character(dna))
		dna = DNAString(dna)

	if(!(class(dna) %in% c("DNAString", "DNAStringSet")))
		stop("The input sequence needs to be either of type DNAString or DNAStringSet")
	

	if(is.null(pwm.rev))
		pwm.rev = reverseComplement(pwm)
		
	# extract only the PWM matrices
	pwm = pwm$pwm
	pwm.rev = pwm.rev$pwm
	
	if(length(dna) < ncol(pwm)){
		stop("DNA sequence needs to be at least as long as the PWM")
	}

	fwd.motif = PWMscoreStartingAt(pwm, dna, starting.at=1:(length(dna)-ncol(pwm)+1))
	back.motif = PWMscoreStartingAt(pwm.rev, dna, starting.at=1:(length(dna)-ncol(pwm)+1))
	
	if(both.strands){
		res = cbind(fwd.motif, back.motif)
		colnames(res) = c("fwd", "rev")
		res = rbind(res, matrix(NA, nrow=(ncol(pwm)-1), ncol=2))
	} else {
		if(strand.fun == "max"){
			res = apply(cbind(fwd.motif, back.motif), 1, max)
		} else if(strand.fun == "mean"){
			res = 2^cbind(fwd.motif, back.motif)
			res = log2(rowMeans(res))
		} else {
			stop("Unknown strand function", strand.fun, ". Valid values are: mean, max")
		}
		# pad out the result so it is of same length as input sequence
		res = c(res, rep(NA, ncol(pwm)-1))
	}
	
	if(odds.score)
		return(2^res)
	else
		return(res)
}

#' Convert frequencies into motifs using PWMUnscaled
#'
#' @param motifs a list of motifs represented as matrices of frequencies (PFM)
#' @param id the set of IDs for the motifs (defaults to names of the 'motifs' list)
#' @param name the set of names for the motifs (defaults to names of the 'motifs' list)
#' @param seq.count if frequencies in the motifs are normalized to 1, provides a vector of sequence counts (e.g. for MotifDb motifs)
#' @param ... other parameters to PWMUnscaled
#'
#' @export
#' @examples
#'
#' if(require("PWMEnrich.Dmelanogaster.background")){
#'    data(MotifDb.Dmel.PFM)
#'
#'    PFMtoPWM(MotifDb.Dmel.PFM) # convert to PWM with uniform background
#'
#'    prior = getBackgroundFrequencies("dm3", quick=TRUE) # get background for drosophila (quick mode on a reduced dataset)
#'    PFMtoPWM(MotifDb.Dmel.PFM, prior.params=prior) # convert with genomic background 
#' }
PFMtoPWM = function(motifs, id=names(motifs), name=names(motifs), seq.count=NULL, ...){
	if(!is.list(motifs)){
		was.list = FALSE
		motifs = list(motifs)
	} else {
		was.list = TRUE
	}
	
	if(is.null(id))
		id = rep("", length(motifs))
	if(is.null(name))
		name = rep("", length(motifs))
	
	if(length(id) != length(motifs))
		stop("The number of IDs (parameter 'id') need to be the same as number of motifs (parameter 'motifs')")

	if(length(name) != length(motifs))
		stop("The number of names (parameter 'name') need to be the same as number of motifs (parameter 'motifs')")
		
	if(!is.null(seq.count) && length(seq.count) != length(motifs)){
		stop("The 'seq.count' vector needs to be of the same length as the list of motifs (parameters 'motifs')")
	}
	
	# call PWMUnscaled
	res = list()
	for(i in 1:length(motifs)){
		if(is.null(seq.count)){
			res[[i]] = PWMUnscaled(motifs[[i]], id=id[i], name=name[i], ...)
		} else {
			res[[i]] = PWMUnscaled(motifs[[i]], id=id[i], name=name[i], seq.count=seq.count[i], ...)
		}
	}
	
	names(res) = names(motifs)
	
	# convert back to single object if that's how the input was
	if(!was.list){
		return(res[[1]])
	} else {
		return(res)
	}

}

#' Normalizes the motifs input argument for multiple functions
#'
#' @param motifs a list of motifs either as frequency matrices (PFM) or as PWM objects. If PFMs are specified
#'               they are converted to PWMs using uniform background. 
.inputParamMotifs = function(motifs){
	# check motifs format and convert to PWM
	if(!is.list(motifs))
		motifs = list(motifs)
		
	if(is.matrix(motifs[[1]]))
		pwms = PFMtoPWM(motifs)
	else if(class(motifs[[1]]) == "PWM")
		pwms = motifs
	else
		stop("motifs need to be either frequency matrices or PWM objects")
		
	return(pwms)
}

#' Normalize the sequences input argument
#'
#' @param sequences a set of sequences to be scanned, a list of DNAString or other scannable objects
.inputParamSequences = function(sequences){
	if(is.character(sequences)){
		sequences = readDNAStringSet(sequences)
	}
	# make sure sequences are in the right formar
	if(!is.list(sequences) & class(sequences) != "DNAStringSet")
		sequences = list(sequences)
	
	if(is.list(sequences) & length(sequences)>0){
		sequences = lapply(sequences, function(s){
			if(is.character(s))
				DNAString(s)
			else
				s
		})
	}
		
	#if(class(sequences) == "DNAStringSet")
	#	sequences = DNAStringSetToList(sequences)
		
	return(sequences)
}

#' Motif affinity of number of hits over a threshold
#'
#' Scan a number of sequences either to find overall affinity, or a number of hits over a score threshold. 
#'
#' @param sequences a set of sequences to be scanned, a list of DNAString or other scannable objects
#' @param motifs a list of motifs either as frequency matrices (PFM) or as PWM objects. If PFMs are specified
#'               they are converted to PWMs using uniform background. 
#' @param raw.scores if to return raw scores (odds) for each position in the sequence. Note that scores for forward and reverse 
#'                   strand are concatenated into a single long vector of scores (twice the length of the sequence)
#' @param verbose if to print verbose output
#' @param cutoff if not NULL, will count number of matches with score above value specified (instead of returning the average affinity). 
#'               Can either be one value, or a vector of values for each of the motifs. 
#' @return if raw.scores=FALSE, returns a matrix of mean scores (after cutoff if any), where columns are motifs. 
#'         The returned values are either mean odd scores (not log-odd), or number of hits above a threshold;
#'         otherwise if raw.scores=TRUE, returns a list of raw score values (before cutoff)
#' @export
#' @examples
#'
#' if(require("PWMEnrich.Dmelanogaster.background")){
#'    data(MotifDb.Dmel)
#'
#'    affinity = motifScores(DNAString("CGTAGGATAAAGTAACTAGTTGATGATGAAAG"), MotifDb.Dmel) # affinity scores
#'    counts = motifScores(DNAString("CGTAGGATAAAGTAACTAGTTGATGATGAAAG"), MotifDb.Dmel, cutoff=log2(exp(4))) # motif hit count with Patser score of 4
#'    print(affinity)
#'    print(counts)
#'
#'    # scanning multiple sequences
#'    sequences = list(DNAString("CGTAGGATAAAGTAACTAGTTGATGATGAAAG"), DNAString("TGAGACGAAGGGGATGAGATGCGGAAGAGTGAAA"))
#'    affinity2 = motifScores(sequences, MotifDb.Dmel)
#'    print(affinity2)
#' }
motifScores = function(sequences, motifs, raw.scores=FALSE, verbose=TRUE, cutoff=NULL){	
	if(!is.null(.PWMEnrich.Options[["useBigMemory"]]) && .PWMEnrich.Options[["useBigMemory"]])
		return(motifScoresBigMemory(sequences, motifs, raw.scores, verbose, cutoff))

	# check motifs format and convert to PWM
	pwms = .inputParamMotifs(motifs)
	sequences = .inputParamSequences(sequences)
		
	pwms.rev = lapply(pwms, reverseComplement)
	
	# always unfold cutoff into a vector of values
	if(!is.null(cutoff)){
		if(length(cutoff) == 1){
			cutoff = rep(cutoff, length(pwms))
		} else if(length(cutoff) != length(pwms)) {
			stop("The lengths of cutoff and pwms do not match. Either provide one values for cutoff, or a vector of values, one for each PWM")
		} 
	}
	
	# inner loop function to use for parallel processing
	motifScoresLoop = function(i){
		if(verbose)
			message(paste("Scanning sequence", i, "/", length(sequences)))
		s = sequences[[i]]
		
		# if using raw scores, record matrix of values;
		# otherwise a single final value for each pwm
		if(raw.scores){
			r = matrix(0, nrow=length(s)*2, ncol=length(pwms))
			colnames(r) = names(pwms)
		} else {
			r = rep(0, length(pwms))
			names(r) = names(pwms)
		}
		
		for(j in 1:length(pwms)){
			# if we are interested in counts, use max over two strands
			
			if(raw.scores){
				# record raw scores that are averages over strands
				r[,j] = as.vector(scanWithPWM(pwms[[j]], s, pwms.rev[[j]], odds.score=TRUE, both.strands=TRUE))
			} else if(!is.null(cutoff)){
				# count number of hits on both strands
				r.pwm = scanWithPWM(pwms[[j]], s, pwms.rev[[j]], odds.score=FALSE, both.strands=TRUE)
				r[j] = sum(r.pwm >= cutoff[j], na.rm=TRUE)
			} else {
				# do average over both strands and also over all values
				r[j] = mean(scanWithPWM(pwms[[j]], s, pwms.rev[[j]], odds.score=TRUE), na.rm=TRUE)
			}
		}
		
		r
	}
	
	####################################### END OF INNER LOOP ################

	# either do it parallel or serial
	if(!is.null(.PWMEnrich.Options[["numCores"]])){
		cat("Parallel scanning with", .PWMEnrich.Options[["numCores"]], "cores\n")
		# do it in parallel
		res = mclapply(1:length(sequences), motifScoresLoop, mc.cores = .PWMEnrich.Options[["numCores"]])
		
		if(is.list(res)){
			if( any(sapply(res, is.null)) ){
				stop("Parallel scanning failed for some sequences. This could be due to a number of reasons including not enough memory.")
			}
		}
	} else {
		# do it serial
		res = lapply(1:length(sequences), motifScoresLoop)
	}
	
	## return either raw scores, or counts or means
	if(raw.scores){
		names(res) = names(sequences)
		return(res)
	} else {
		if(length(pwms) == 1){
			return( matrix(sapply(res, identity), ncol=1, dimnames=list(NULL, names(pwms))) )
		} else {
			r = t(sapply(res, identity))
			rownames(r) = names(sequences)
			colnames(r) = names(pwms)
			return( r )
		}
	}
}

#' This is a memory intensive version of motifScore() which is about 2 times faster
#'
#' The parameters and functionality are the same as \code{\link{motifScores}}. Please refer to documentation of this function
#' for detailed explanation of functionality. 
#'
#' This function is not meant to be called directly, but is indirectly called by motifScores() once a global parameters useBigMemory is set. 
#'
#' @param sequences set of input sequences
#' @param motifs set of input PWMs or PFMs
#' @param raw.scores if to return scores for each base-pair
#' @param verbose if to produce verbose output
#' @param cutoff the cutoff for calling binding sites (in base 2 log). 
#' @param seq.all already concatenated sequences if already available (used to internally speed up things)
#'
#' @seealso \code{\link{motifScores}}
motifScoresBigMemory = function(sequences, motifs, raw.scores=FALSE, verbose=TRUE, cutoff=NULL, seq.all=NULL){	
	# check motifs format and convert to PWM
	pwms = .inputParamMotifs(motifs)
	sequences = .inputParamSequences(sequences)
		
	pwms.rev = lapply(pwms, reverseComplement)
	pwms.len = sapply(pwms, length)
	min.pwm.len = min(pwms.len)
	
	# always unfold cutoff into a vector of values
	if(!is.null(cutoff)){
		if(length(cutoff) == 1){
			cutoff = rep(cutoff, length(pwms))
		} else if(length(cutoff) != length(pwms)) {
			stop("The lengths of 'cutoff' and 'pwms' do not match. Either provide one values for 'cutoff' or a vector of values, one for each PWM")
		} 
		
		# cutoff in odds space
		cutoff.2 = 2^cutoff
	}
	
	# work on a single sequence that is all concatenated together because it's faster
	if(is.null(seq.all))
		seq.all = DNAString(concatenateSequences(sequences))
	seq.len = sapply(sequences, length)
	
	# a common error with shorter sequences
	shorter.seq = which(seq.len < min.pwm.len)
	if(length(shorter.seq)>0){
		stop(paste(length(shorter.seq), "sequence(s) have length shorter than the shortest PWM (", min.pwm.len, "), please identify and remove these sequences."))
	}
	
	# a grouping data frame
	seq.len.sum = cumsum(seq.len)
	seq.group = data.frame("from"=rep(0, length(sequences)), "to"=0)
	seq.group$from = c(1, 1+seq.len.sum[-length(seq.len.sum)])
	seq.group$to = seq.len.sum
	
	# inner loop function to use for parallel processing
	#
	# @param motif.start the index of the first motif
	# @param motif.end the index of the last motif 
	motifScoresLoop = function(input.params){
		motif.start = input.params[[1]]
		motif.end = input.params[[2]]
	
		if(verbose){
			if(motif.start == motif.end){
				message(paste("Starting scanning with motif", names(pwms)[motif.start]))
			} else {
				message(paste("Starting scanning with motifs from", motif.start, "to", motif.end))
			}
		}
		# local copy of the sequences and motifs used in the function
		s = seq.all
		s.group = seq.group
		pwms.l = pwms[motif.start:motif.end]
		pwms.rev.l = pwms.rev[motif.start:motif.end]
		pwms.len.l = pwms.len[motif.start:motif.end]
		if(!is.null(cutoff))
			cutoff.2.l = cutoff.2[motif.start:motif.end]
		
		# create the output structure
		res = list()
		for(i in 1:nrow(s.group)){
			if(raw.scores){			
				r = matrix(0, nrow=seq.len[i]*2, ncol=length(pwms.l))
				colnames(r) = names(pwms.l)
				res[[i]] = r		
			} else {
				r = rep(0, length(pwms.l))
				names(r) = names(pwms.l)		
				res[[i]] = r
			}
		}
		
		### do the scanning, convert back to individual sequences and do the post-processing
		for(j in 1:length(pwms.l)){
			if(verbose){
				message(paste("Scanning all sequences with motif", j, "/", length(pwms.l)))
			}
			raw = as.vector(scanWithPWM(pwms.l[[j]], s, pwms.rev.l[[j]], odds.score=TRUE, both.strands=TRUE))
			
			### split the raw results by sequence
			for(i in 1:nrow(s.group)){
				# selector for this sequence
				from = seq.group$from[i]
				to = seq.group$to[i]
				from.rev = length(raw)/2 + seq.group$from[i]
				to.rev = length(raw)/2 + seq.group$to[i]
				sel = c(from:to, from.rev:to.rev)
			
				r = raw[sel]
				
				# fwd strand
				r[((seq.len[i]-pwms.len.l[j])+2):seq.len[i]] = NA
				# rev strand
				r[((2*seq.len[i]-pwms.len.l[j])+2):(2*seq.len[i])] = NA

				# save the values				
				if(raw.scores){
					res[[i]][,j] = r
				} else if(!is.null(cutoff)){
					# count number of hits on both strands 
					res[[i]][j] = sum(r >= cutoff.2.l[j], na.rm=TRUE)
				} else {
					# MeanAffinity
					res[[i]][j] = mean(r, na.rm=TRUE)
				}
				
			}
		}
				
		res
	}
	######################################## END OF INNER LOOP ######################################

	# either do it parallel or serial
	if(!is.null(.PWMEnrich.Options[["numCores"]]) && length(pwms) >= .PWMEnrich.Options[["numCores"]]){
		cores = .PWMEnrich.Options[["numCores"]]
		sel = round(seq(0, length(pwms), length.out=cores+1))
		start = sel[1:cores]+1
		end = sel[2:(cores+1)]
		input = lapply(1:cores, function(i) list(start[i], end[i]))
		# do it in parallel
		res.parallel = mclapply(input, motifScoresLoop, mc.cores = cores)
		
		if(is.list(res.parallel)){
			if( any(sapply(res.parallel, is.null)) ){
				stop("Parallel scanning failed for some sequences. This could be due to a number of reasons including not enough memory.")
			}
		}	

		# concatenate results for different motifs		
		if(raw.scores)
			res = lapply(1:length(sequences), function(i) do.call("cbind", lapply(res.parallel, function(x) x[[i]])))	
		else
			res = lapply(1:length(sequences), function(i) do.call("c", lapply(res.parallel, function(x) x[[i]])))
	} else {
		# do it serial
		res = motifScoresLoop(list(1, length(pwms)))
	}
	
	
	
	## return either raw scores, or counts or means
	if(raw.scores){
		names(res) = names(sequences)
		return(res)
	} else {
		if(length(pwms) == 1){
			return( matrix(sapply(res, identity), ncol=1, dimnames=list(names(sequences), names(pwms))) )
		} else {
			r = t(sapply(res, identity))
			rownames(r) = names(sequences)
			colnames(r) = names(pwms)
			return( r )
		}
	}
}


#' Information content for a PWM or PFM
#'
#' @param motif a matrix of frequencies, or a PWM object
#' @param prior.params the prior parameters to use when a matrix is given (ignored if motif is already a PWM)
#' @param bycol if to return values separately for each column
#' @return information content in bits (i.e. log2)
#' @export
#' @examples
#'
#' if(require("PWMEnrich.Dmelanogaster.background")){
#'    data(MotifDb.Dmel)
#'    data(MotifDb.Dmel.PFM)
#'
#'    motifIC(MotifDb.Dmel$ttk) # the nucleotide distribution is taken from the PWM (in this case genomic background)
#'    motifIC(MotifDb.Dmel.PFM$ttk) # information content with default uniform background because the input is a matrix, not PWM object
#' }
motifIC = function(motif, prior.params=c(A=0.25, C=0.25, G=0.25, T=0.25), bycol=FALSE){
	if(class(motif) == "PWM"){
		p = PWMUnscaled(motif$pfm, type="prob", prior.params=motif$prior.params)
		bg = motif$prior.params
	} else {
		p = PWMUnscaled(motif, type="prob", prior.params=prior.params)
		bg = prior.params
	}
		
	bg = bg/sum(bg) # normalize to 1
	
	cols = p$pwm * log2(p$pwm/bg)
	if(bycol)
		return(colSums(cols))
	else
		return(sum(cols))
}

#' Calculate motif enrichment using one of available scoring algorithms and background corrections. 
#'
#' This function provides and interface to all algorithms available in PWMEnrich
#' to find motif enrichment in a single or a group of sequences with/without
#' background correction. 
#'
#' Since for all algorithms the first step involves calculating raw scores without background correction, the output
#' always contains the scores without background correction together with (optional) background-corrected
#' scores. 
#'
#' Unless otherwise specified the scores are returned both separately for each sequence (without/with background) and
#' for the whole group of sequences (without/with background). 
#'
#' To use a background correction you need to supply a set of PWMs with precompiled background distribution parameters
#' (see function \code{\link{makeBackground}}). When such an object is supplied as the \code{pwm} parameter, the scoring 
#' scheme and background correction are automatically determined. 
#'
#' There are additional packages with already pre-computed background (e.g. see package \code{PWMEnrich.Dmelanogaster.background}).
#'
#' Please refer to (Stojnic & Adryan, 2012) for more details on the algorithms. 
#'
#' @title Motif enrichment
#' @param sequences the sequences to be scanned for enrichment. Can be either a single sequence 
#'        (an object of class DNAString), or a list of DNAString objects, or a DNAStringSet object.
#' @param pwms this parameter can take multiple values depending on the scoring scheme and background correction used.
#'             When the \code{method} parameter is set to "autodetect", the following default algorithms
#'             are going to be used: 
#'        \itemize{
#'           \item if \code{pwms} is a list containing either frequency matrices or a list of PWM objects then
#'                 the "affinity" algorithm is selected. If frequency matrices are given, they are converted 
#'                 to PWMs using uniform background. For best performance, convert frequency matrices to PWMs
#'                 before calling this function using realistic genomic background. 
#'           \item Otherwise, appropriate scoring scheme and background correction are selected based on the
#'                 class of the object (see below). 
#'        }
#' @param score this parameter determines which scoring scheme to use. Following scheme as available:
#'        \itemize{
#'           \item "autodetect" - default value. Scoring method is determined based
#'                 on the type of \code{pwms} parameter. 
#'           \item "affinity" - use threshold-free affinity score. The \code{pwms}
#'                 parameter can either be a list of frequency matrices, \code{PWM} objects, or a 
#'                 \code{PWMLognBackground} object. 
#'           \item "cutoff" - use number of motif hits above a score cutoff. The \code{pwms}
#'                 parameter can either be a list of frequency matrices, \code{PWM} objects, or a 
#'                 \code{PWMCutoffBackground} object.
#'           \item "clover" - use the Clover algorithm (Frith et al, 2004). The Clover score of a single 
#'                 sequence is identical to the affinity score, while for a group of sequences is an
#'                 average of products of affinities over all sequence subsets. 
#'         }
#'            
#' @param bg this parameter determines how the raw score is compared to the background distribution. 
#'        \itemize{
#'           \item "autodetect" - default value. Background correction is determined based on the type
#'                  of the \code{pwms} parameter. 
#'           \item "logn" - use a lognormal distribution background pre-computed for a set of PWMs.
#'                 This requires \code{pwms} to be of class \code{PWMLognBackground}.
#'           \item "z" - use a z-score for the number of significant motif hits compared to background number of hits.
#'                 This requires \code{pwms} to be of class \code{PWMCutoffBackground}.
#'           \item "pval" - use empirical P-value based on a set of background sequences. This requires
#'                 \code{pwms} to be of class \code{PWMEmpiricalBackground}. Note that PWMEmpiricalBackground
#'                 objects tend to be very large so that the empirical P-value can be calculated in reasonable time.
#'           \item "ms" - shuffle columns of motif matrices and use that as basis for P-value calculation. Note that
#'                  since the sequences need to rescanned with all of the new shuffled motifs this can be very slow. 
#'                  Also, this also works only no *individual* sequences, not groups. 
#'           \item "none" - no background correction
#'         }
#' @param cutoff the score cutoff for a significant motif hit if scoring scheme "cutoff" is selected. 
#' @param verbose if to print verbose output
#' @param motif.shuffles number of times to shuffle motifs if using "ms" background correction
#' @param B number of replicates when calculating empirical P-value
#' @param group.only if to return statistics only for the group of sequences, not individual sequences. In the case of
#'                   empirical background the P-values for individual sequences are not calculated (thus saving time), for other 
#'                   backgrounds they are calculated but not returned.  
#' @return a MotifEnrichmentResults object containing a subset following elements:
#'         \itemize{
#'           \item "score" - scoring scheme used
#'           \item "bg" - background correction used
#'           \item "params" - any additional parameters
#'			 \item "sequences" - the set of sequences used
#'           \item "pwms" - the set of pwms used
#'           \item "sequence.nobg" - per-sequence scores without any background correction. 
#'                   For "affinity" and "clover" a matrix of mean affinity scores; for
#'                   "cutoff" number of significant hits above a cutoff 
#'           \item "sequence.bg" - per-sequence scores after background correction. For "logn" and "pval" the P-value (smaller is better);
#'                                 for "z" and "ms" background corrections the z-scores (bigger is better). 
#'           \item "group.nobg" - aggregate scores for the whole group of sequences without background correction. For "affinity"
#'                                and "clover" the mean affinity over all sequences in the set; for "cutoff" the total number of hits in all
#'                                sequences.
#'           \item "group.bg" - aggregate scores for the whole group of sequences with background correction. For "logn" and "pval",
#'                              the P-value for the whole group (smaller is better); for "z" and "ms" the z-score for the whole set (bigger is better).
#'           \item "sequence.norm" - (only for "logn") the length-normalized scores for each of the sequences. Currently only implemented
#'                                   for "logn", where it returns the values normalized from LogN(0,1) distribution
#'           \item "group.norm" - (only for "logn") similar to sequence.norm, but for the whole group of sequences
#'         }
#' @export
#' @references \itemize{ 
#'  \item R. Stojnic & B. Adryan: Identification of functional DNA motifs using a binding affinity lognormal background distribution, submitted.
#' 	\item MC Frith et al: Detection of functional DNA motifs via statistical over-representation, Nucleid Acid Research (2004).
#' }
#' @examples 
#'
#' if(require("PWMEnrich.Dmelanogaster.background")){
#'    ###
#'    # load the pre-compiled lognormal background
#'    data(PWMLogn.dm3.MotifDb.Dmel)
#'
#'    # scan two sequences for motif enrichment
#'    sequences = list(DNAString("GAAGTATCAAGTGACCAGTAGATTGAAGTAGACCAGTC"), DNAString("AGGTAGATAGAACAGTAGGCAATGGGGGAAATTGAGAGTC"))
#'    res = motifEnrichment(sequences, PWMLogn.dm3.MotifDb.Dmel)
#'
#'    # most enriched in both sequences (lognormal background P-value)
#'    head(motifRankingForGroup(res))
#'
#'    # most enriched in both sequences (raw affinity, no background)
#'    head(motifRankingForGroup(res, bg=FALSE))
#'
#'    # most enriched in the first sequence (lognormal background P-value)
#'    head(motifRankingForSequence(res, 1))
#'
#'    # most enriched in the first sequence (raw affinity, no background)
#'    head(motifRankingForSequence(res, 1, bg=FALSE))
#'
#'    ###
#'    # Load the pre-compiled background for hit-based motif counts with cutoff of P-value = 0.001 
#'    data(PWMPvalueCutoff1e3.dm3.MotifDb.Dmel)
#'
#'    res.count = motifEnrichment(sequences, PWMPvalueCutoff1e3.dm3.MotifDb.Dmel)
#'
#'    # Enrichment in the whole group, z-score for the number of motif hits
#'    head(motifRankingForGroup(res))
#'
#'    # First sequence, sorted by number of motif hits with P-value < 0.001
#'    head(motifRankingForSequence(res, 1, bg=FALSE))
#'    
#' }
motifEnrichment = function(sequences, pwms, score="autodetect", bg="autodetect", cutoff=NULL, verbose=TRUE, motif.shuffles=30, B=1000,
	group.only=FALSE){
	
	# detect scoring method
	if(score == "autodetect"){
		if(class(pwms) == "PWMLognBackground"){
			score = "affinity"
		} else if(class(pwms) == "PWMCutoffBackground"){
			score = "cutoff"
		} else {
			score = "affinity"
		}
	}
	
	if(score == "cutoff" & is.null(cutoff))
		cutoff = pwms@bg.cutoff
	
	# the PWM* object
	pwmobj = NULL
	
	# detect background correction
	if(!(bg %in% c("none", "ms"))){
		if(class(pwms) == "PWMLognBackground"){
			bg = "logn"
			pwmobj = pwms
			pwms = pwmobj@pwms
		} else if(class(pwms) == "PWMCutoffBackground"){
			bg = "z"
			pwmobj = pwms
			pwms = pwmobj@pwms
		} else if(class(pwms) == "PWMEmpiricalBackground"){
			bg = "pval"
			pwmobj = pwms
			pwms = pwmobj@pwms
		} else if(class(pwms) == "PWMGEVBackground"){
			bg = "gev"
			pwmobj = pwms
			pwms = pwmobj@pwms
		} else {
			bg = "none"
		}
	}

	if(!is.list(pwms))
		pwms = list(pwms)
	
	sequences = .inputParamSequences(sequences)
	
	# return list
	res = list(score=score, bg=bg, pwms=pwms, sequences=sequences, params=list())
	
	# apply no-bg scoring
	if(score == "affinity"){
		# lengths of sequences and PWMs
		seq.len = sapply(sequences, length)
		pwm.len = sapply(pwms, length)
		
		res$sequence.nobg = motifScores(sequences, pwms, verbose=verbose)
		res$group.nobg = affinitySequenceSet(res$sequence.nobg, seq.len, pwm.len)
	} else if(score == "clover"){
		res$sequence.nobg = motifScores(sequences, pwms, verbose=verbose)
		res$group.nobg = cloverScore(res$sequence.nobg, verbose=verbose)
	} else if(score == "cutoff"){
		res$params = list(cutoff=cutoff)
		res$sequence.nobg = motifScores(sequences, pwms, cutoff=cutoff, verbose=verbose)
		res$group.nobg = colSums(res$sequence.nobg)
	} else {
		stop(paste("Unknown scoring algorithm: '", score, "'. Please select one of: 'affinity', 'cutoff', 'clover'", sep=""))
	}
	
	# apply background correction if needed
	if(bg == "none"){
		res$sequence.bg = NULL
		res$group.bg = NULL
	} else if(bg == "logn"){
		# lengths of sequences and PWMs
		seq.len = sapply(sequences, length)
		pwm.len = sapply(pwms, length)
		# do it per sequence
		res$sequence.bg = t(sapply(1:length(seq.len), 
			function(i) logNormPval(res$sequence.nobg[i,], seq.len[i], pwm.len, pwmobj@bg.mean, pwmobj@bg.sd, pwmobj@bg.len)))
		colnames(res$sequence.bg) = names(pwms)
		# convert logn P-values into normalized observations and run Clover on that
		res$sequence.norm = apply(res$sequence.bg, 1:2, qlnorm, lower.tail=FALSE)
		
		if(score == "affinity"){
			# and for the group
			#### these two lines are the old implementation!!!
			# res$group.bg = logNormPvalSequenceSet(res$sequence.nobg, seq.len, pwm.len, pwmobj@bg.mean, pwmobj@bg.sd, pwmobj@bg.len)
			# res$group.norm = sapply(res$group.bg, qlnorm, lower.tail=FALSE)
			#####
			res$group.bg = colMeans(res$sequence.norm)
			res$group.norm = res$group.bg
			
		} else if(score == "clover"){
			# now run Clover on normalized scores
			res$group.bg = cloverScore(res$sequence.norm, verbose=verbose)
		}
	} else if(bg == "z"){
		# lengths of sequences and PWMs
		seq.len = sapply(sequences, length)
		pwm.len = sapply(pwms, length)
		res$sequence.bg = cutoffZscore(res$sequence.nobg, seq.len, pwm.len, pwmobj@bg.P)
		res$group.bg = cutoffZscoreSequenceSet(res$sequence.nobg, seq.len, pwm.len, pwmobj@bg.P)
	} else if(bg == "pval"){
		seq.len = sapply(sequences, length)
		pwm.len = sapply(pwms, length)
		
		# if to use cutoff
		usecutoff = NULL
		if(score == "cutoff")
			usecutoff = cutoff
		
		if(!group.only){
			res$sequence.bg = t(sapply(1:length(seq.len), function(i)
				 empiricalPvalue(res$sequence.nobg[i,], seq.len[i], pwm.len, pwmobj@bg.fwd, pwmobj@bg.rev, cutoff=usecutoff, B=B, verbose=verbose)))
		}
		
		# do empirical P-value
		if(score == "clover")
			res$group.bg = cloverPvalue1seq(res$sequence.nobg, seq.len, pwm.len, pwmobj@bg.fwd, pwmobj@bg.rev, B=B, verbose=verbose, clover=res$group.nobg)
		else	 
			res$group.bg = empiricalPvalueSequenceSet(res$sequence.nobg, seq.len, pwm.len, pwmobj@bg.fwd, pwmobj@bg.rev, cutoff=usecutoff, B=B, verbose=verbose)
	} else if(bg == "ms"){
		# if to use cutoff
		usecutoff = NULL
		if(score == "cutoff")
			usecutoff = cutoff
		res$sequence.bg = matrixShuffleZscorePerSequence(res$sequence.nobg, sequences, pwms, cutoff=usecutoff, B=motif.shuffles)
	} else if(bg == "gev"){
		seq.len = sapply(sequences, length)
		pwm.len = sapply(pwms, length)
		
		res$sequence.bg = gevPerSequence(res$sequence.nobg, seq.len, pwm.len, pwmobj@bg.loc, pwmobj@bg.scale, pwmobj@bg.shape)
		res$group.bg = NULL
	} else {
		stop(paste("Uknown background correction algorithm: '", bg, "', Please select one of: 'none', 'logn', 'z', 'pval', 'ms'", sep=""))
	}
	
	# only store the values for the group, not individual sequences!
	if(group.only){
		seq.n = grep("^sequence[.]", names(res))
		if(length(seq.n)>0){
			res = res[-seq.n]
		}
	}
	
	# restore the sequence names
	if("sequence.nobg" %in% names(res)){
		rownames(res$sequence.nobg) = names(sequences)
	}
	if("sequence.bg" %in% names(res)){
		rownames(res$sequence.bg) = names(sequences)
	}
	
	return(new("MotifEnrichmentResults", res=res))
		
}

#' Obtain z-score for motif column shuffling
#'
#' All PWMs are shuffled at the same time. This function would be too slow to produce
#' empirical P-values, thus we return a z-score from a small number of shuffles. 
#'
#' The z-scores are calculated for each sequence individually. 
#' 
#' @param scores a set of already calculated scores
#' @param sequences either one sequence or a list/set of sequences (objects of type DNAString or DNAStringSet)
#' @param pwms a list of PWMs
#' @param cutoff if NULL, will use affinity, otherwise will use number of hits over this log2 odds cutoff
#' @param B number of replicates, i.e. PWM column shuffles
matrixShuffleZscorePerSequence = function(scores, sequences, pwms, cutoff=NULL, B=30){
	# scores from different replicates for each of the sequences
	scores.matrix = array(0, dim=c(nrow(scores), ncol(scores), B))
	dimnames(scores.matrix) = list(rownames(scores), colnames(scores), NULL)
	
	# perfom shuffles	
	for(i in 1:B){
		cat("Perfoming shuffle", i, "/", B, "\n")
		pwms.shuffled = pwms
		for(j in 1:length(pwms)){
			p = pwms[[j]]
			pwms.shuffled[[j]]@pwm = p@pwm[, sample(ncol(p@pwm))]
		}
		
		scores.suffled = motifScores(sequences, pwms.shuffled, cutoff=cutoff, verbose=FALSE)
		scores.matrix[,,i] = scores.suffled
		
	}
	
	# work out mean and sd
	shuffled.mean = apply(scores.matrix, c(1,2), mean)
	shuffled.sd = apply(scores.matrix, c(1,2), sd)
	
	return( (scores - shuffled.mean) / shuffled.sd )
}


#' Calculate the P-value from lognormal distribution with background of equal length
#' 
#' @param scores affinity scores for the PWMs, can contain scores for more than one sequence (as rows), P-values are extracted separately
#' @param seq.len the length distribution of the sequences
#' @param pwm.len the leggths of PWMs
#' @param bg.mean the mean values from the background for PWMs
#' @param bg.sd the sd values from the background
#' @param bg.len the length distribution of the background (we currently support only constant length)
logNormPval = function(scores, seq.len, pwm.len, bg.mean, bg.sd, bg.len){
	if(is.vector(scores)){
		scores = matrix(scores, nrow=1, dimnames=list(NULL, names(scores)))
	}
	if( length(seq.len) != nrow(scores) ){
		stop("The sequence length distribution and number of rows for scores doesn't match")
	}
	if( !all(colnames(scores) == names(bg.mean)) ){
		stop("The motifs column names in 'scores' and 'bg.mean' do not match")
	}
	
	# calculate the expected sd values for the length distribution for each sequence
	seq.sd = matrix(0, nrow=nrow(scores), ncol=ncol(scores))

	for(i in 1:nrow(seq.sd)){
		seq.sd[i,] = bg.sd / sqrt( (seq.len[i]-pwm.len+1) / bg.len)
	}
	
	# calculate p-values
	seq.p = matrix(0, nrow=nrow(scores), ncol=ncol(scores))
	colnames(seq.p) = colnames(scores)
	rownames(seq.p) = rownames(scores)
	
	for(i in 1:nrow(seq.sd)){
		for(j in 1:ncol(seq.sd)){
			# calculate the mean/sd parameters of the lognormal distribution
			mx = bg.mean[j]
			sx = seq.sd[i,j]
			
			ml = 2*log(mx) - 0.5*log(mx^2+sx^2)
			sl = sqrt(-2*log(mx) + log(mx^2+sx^2))
			
			seq.p[i,j] = plnorm(scores[i,j], meanlog=ml, sdlog=sl, lower.tail=FALSE)
		}
	}
	
	seq.p
}

#' Lognormal P-value for a set of sequences
#'
#' @param scores a matrix of per-sequence affinity scores
#' @param seq.len lengths of sequences
#' @param pwm.len lengths of pwms
#' @param bg.mean mean background at length of bg.len
#' @param bg.sd standard deviation of background at length of bg.len
#' @param bg.len the length for which mean and sd are calculated
#' @return P-value
logNormPvalSequenceSet = function(scores, seq.len, pwm.len, bg.mean, bg.sd, bg.len){
	# total score for all the sequences
	s = structure(rep(0, ncol(scores)), names=colnames(scores))
	res = structure(rep(0, ncol(scores)), names=colnames(scores))
	for(i in 1:ncol(scores)){
		# lengths of sequences for this pwm, all of them are shorter by pwm.len-1
		seq.len.pwm = seq.len - pwm.len[i] + 1
		
		# sum over sequences then average by total length
		s[i] = sum(scores[,i] * seq.len.pwm) / sum(seq.len.pwm)
		res[i] = logNormPval(s[i], sum(seq.len.pwm)+pwm.len[i]-1, pwm.len[i], bg.mean[i], bg.sd[i], bg.len[i])
	}
	
 	return(res)
}


#' Replace all infinite values by 0
#'
#' @param x a vector of values
keepFinite = function(x) {
	if(any(!is.finite(x)))
		x[!is.finite(x)] = 0
	x
}

#' Z-score calculation for cutoff hits
#'
#' The Z-score is calculated separately for each sequence
#'
#' @param scores the hit counts for the sequences
#' @param seq.len the length distribution of sequences
#' @param pwm.len the length distribution of the PWMs
#' @param bg.P background probabilities of observing a motif hit at nucleotide resolution 
#'             (scaled to sequence length, not 2 * length)
#' @return Z-score
cutoffZscore = function(scores, seq.len, pwm.len, bg.P){
	if(is.vector(scores))
		scores = matrix(scores, nrow=1)

	# actual probability is 1/2 because we have two strands
	bg.P = bg.P / 2
	
	z = matrix(0, nrow=nrow(scores), ncol=ncol(scores))
	rownames(z) = rownames(scores)
	colnames(z) = colnames(scores)
	for(i in 1:nrow(scores)){
		# actualy length is 2x because we have two strands
		P.n = 2 * (seq.len[i] - pwm.len + 1)
		P.mean = bg.P * P.n
		P.sd = sqrt(P.n * bg.P * (1-bg.P))
		# so this would be the poisson approximation
		# z[i,] = 1 / ppois(scores[i,], P.mean, lower.tail=FALSE)
		z[i,] = keepFinite( (scores[i,] - P.mean) / P.sd )
	}
	
	return(z)
}

#' Z-score calculation for cutoff hits for group of sequences
#'
#' The Z-score is calculated as if the sequence came for one very long sequence
#'
#' @param scores the hit counts for the sequences
#' @param seq.len the length distribution of sequences
#' @param pwm.len the length distribution of the PWMs
#' @param bg.P background probabilities of observing a motif hit at nucleotide resolution
#' @return Z-score
cutoffZscoreSequenceSet = function(scores, seq.len, pwm.len, bg.P){	
	total.hits = colSums(scores)
	total.len = sum(seq.len) - length(seq.len)*(pwm.len-1) + pwm.len - 1
	
	# do it seperate for each PWM to make sure the lengths are right
	res = sapply(1:length(total.hits), function(i) cutoffZscore(total.hits[i], total.len[i], pwm.len[i], bg.P[i]))
	names(res) = names(total.hits)
	
	return(res)
}

#' Calculate the empirical P-value by affinity of cutoff. 
#'
#' This is the new backend function for empirical P-values for either affinity or cutoff. 
#' The function only works on single sequences. 
#'
#' @param scores the scores obtained for the sequence
#' @param seq.len the length of the sequence, if a single value will take a single sequence
#'                of given length. If a vector of values, will take sequences of given lengths
#'                and joint them together
#' @param pwm.len the lengths of PWMs
#' @param bg.fwd raw odds scores for the forward strand of background
#' @param bg.rev raw odds scores for the reverse strand of background
#' @param cutoff if not NULL, will use hit count above this cutoff. The cutoff should be specified in log2. 
#' @param B the number of random replicates
#' @param verbose if to give verbose progress reports
#' @param exact.length if to take into consideration that the actual sequence lengths differ for different PWMs.
#'                     For very long sequences (i.e. seq.len >> pwm.len) this make very little difference, however
#'                     the run time with exact.length is much longer. 
empiricalPvalue = function(scores, seq.len, pwm.len, bg.fwd, bg.rev, cutoff=NULL, B=10000, verbose=FALSE, exact.length=FALSE){
	# length of sequences
	bg.len = nrow(bg.fwd)
	
	if(max(seq.len) > bg.len){
		stop("The length of sequence greater than background!")
	} 
	
	if(max(seq.len) * B > bg.len){
		warning("The length of the sequence multiplied by number of runs (B) is greater than background length. This might lead to unreliable P-value estimates.")
	}
		
	# the empricial p-value
	score.p = rep(0, length(scores))
		
	for(i in 1:B){
		if(verbose)
			message("Random sample ", i, " / ", B)
		
		# calculate the final score as sum over different seq.len sequences
		final.score = 0
		for(k in 1:length(seq.len)){
			# choose a subsequence with matching length
			start.range = 1:(bg.len - seq.len[k])
			start = sample(start.range, 1)
		
			# fetch the scores from the two strands
			score.fwd = bg.fwd[start:(start+seq.len[k]-1),]
			score.rev = bg.rev[start:(start+seq.len[k]-1),]

			# simulate that we don't have scores for last pwm.len+1 nucleotides
			if(exact.length){
				for(j in 1:length(pwm.len)){
					score.fwd[(nrow(score.fwd)-pwm.len[j]+1):nrow(score.fwd),j] = NA
					score.rev[(nrow(score.rev)-pwm.len[j]+1):nrow(score.rev),j] = NA
				}
			}
			
			# append both strands to one long string		
			score.both = rbind(score.fwd, score.rev)
		
			if(is.null(cutoff)){
				# affinity algorithm
				new.score = colSums(score.both, na.rm=TRUE)
			} else {
				# number of hits above threshold
				new.score = colSums(score.both >= 2^cutoff, na.rm=TRUE)
			}
			final.score = final.score + new.score
		}
		
		# to obtain the mean we divide by twice the length (because it's both strands)
		final.score = final.score / (2 * sum(seq.len)) 
		
		score.p = score.p + ((final.score > scores) + 0)
		
	}
	
	return(score.p / B)
}

#' Empirical P-value for a set of sequences
#'
#' Calculate empirical P-value for a set of sequences, using either affinity or cutoff. When cutoff is used, the score
#' is a number of motif hits above a certain log-odds cutoff. 
#'
#' @param scores a matrix of scores, rows for sequences, columns for PWMs
#' @param seq.len the lengths of sequences
#' @param pwm.len the lengths of PWMs
#' @param bg.fwd raw odds scores for the forward strand of background
#' @param bg.rev raw odds scores for the reverse strand of background
#' @param cutoff if not NULL, will use hit count above this cutoff. The cutoff should be specified in log2. 
#' @param B the number of random replicates
#' @param verbose if to give verbose progress reports
empiricalPvalueSequenceSet = function(scores, seq.len, pwm.len, bg.fwd, bg.rev, cutoff=NULL, B=10000, verbose=FALSE){
	# affinity
	if(is.null(cutoff)){
		# total score for all the sequences
		s = structure(rep(0, ncol(scores)), names=colnames(scores))
		res = structure(rep(0, ncol(scores)), names=colnames(scores))
		exact.lengths = structure(rep(0, ncol(scores)), names=colnames(scores))
		for(i in 1:ncol(scores)){
			# lengths of sequences for this pwm, all of them are shorter by pwm.len-1
			seq.len.pwm = seq.len - pwm.len[i] + 1
			exact.lengths[i] = sum( seq.len.pwm )
		
			# sum over sequences then average by total length
			s[i] = sum(scores[,i] * seq.len.pwm) / sum(seq.len.pwm)			
		}
		
		# approximate sequence length without all the last positions being matched
		total.len = seq.len - round(mean(pwm.len))
		
		# run it on mean sequence length
		res = empiricalPvalue(s, total.len, pwm.len, bg.fwd, bg.rev, cutoff, B, verbose, exact.length=FALSE)
	
	 	return(res)
	} else {
		total.hits = colSums(scores)
		total.len = seq.len - round(mean(pwm.len)) #sum(seq.len) - length(seq.len)*pwm.len
	
		# run ti on mean sequence lengths (with pwm.len)
		res = empiricalPvalue(total.hits, total.len, pwm.len, 
					bg.fwd, bg.rev, cutoff, B, verbose, exact.length=FALSE)
					
		return(res)
	}
}

#' Apply GEV background normalization per every sequence
#' 
#' @param scores affinity scores for the PWMs, can contain scores for more than one sequence (as rows), P-values are extracted separately
#' @param seq.len the length distribution of the sequences
#' @param pwm.len the lengths of PWMs
#' @param bg.loc list of linear regression for location parameter
#' @param bg.scale list of linear regression for scale parameter
#' @param bg.shape list of linear regression for shape parameter
gevPerSequence = function(scores, seq.len, pwm.len, bg.loc, bg.scale, bg.shape){
	if(is.vector(scores)){
		scores = matrix(scores, nrow=1, dimnames=list(NULL, names(scores)))
	}
	if( length(seq.len) != nrow(scores) ){
		stop("The sequence length distribution and number of rows for scores doesn't match")
	}
	
	# calculate the expected sd values for the length distribution for each sequence
	seq.sd = matrix(0, nrow=nrow(scores), ncol=ncol(scores))
	
	log.seqlen = data.frame("log.len"=log(seq.len))
	
	# predict loc, scale, shape for each PWM
	loc = matrix(0, nrow=nrow(scores), ncol=ncol(scores))
	scale = matrix(0, nrow=nrow(scores), ncol=ncol(scores))
	shape = matrix(0, nrow=nrow(scores), ncol=ncol(scores))
	for(i in 1:ncol(scores)){
		loc[,i] = as.vector(predict.lm(bg.loc[[i]], newdata=log.seqlen))
		scale[,i] = as.vector(predict.lm(bg.scale[[i]], newdata=log.seqlen))
		shape[,i] = as.vector(predict.lm(bg.shape[[i]], newdata=log.seqlen))
	}
	
	## work out the P-values
	seq.p = matrix(0, nrow=nrow(scores), ncol=ncol(scores))
	colnames(seq.p) = colnames(scores)
	rownames(seq.p) = rownames(scores)
	
	for(i in 1:nrow(seq.p)){
		for(j in 1:ncol(seq.p)){
			# if scale is <0 the GEV approximation is no longer valid
			if(scale[i,j] < 0)
				seq.p[i,j] = NA
			else
				seq.p[i,j] = pgev(log(scores[i,j]), loc=loc[i,j], scale=scale[i,j], shape=shape[i,j], lower.tail=FALSE)
		}
	}
	
	seq.p
}

#' Calculate total affinity over a set of sequences
#'
#' @param scores affinity scores for individual sequences
#' @param seq.len lengths of sequences
#' @param pwm.len lengths of PWMs
affinitySequenceSet = function(scores, seq.len, pwm.len){
	if(is.vector(scores))
		scores = matrix(scores, nrow=1, dimnames=list(NULL, names(scores)))

	# final scores for the group
	final = structure(rep(0, ncol(scores)), names=colnames(scores))
	
	for(i in 1:ncol(scores)){
		# actual length of scanned sequence
		seq.len.pwm = seq.len - pwm.len[i] + 1
		
		final[i] = sum(scores[,i] * seq.len.pwm) / sum(seq.len.pwm)
	}
	
	return(final)
}

#' Calculate Recovery-AUC for motifs ranked according to some scoring scheme
#'
#' Note that this function asssumes that smaller values are better!
#' 
#' @param seq.res a matrix where each column represents a PWM and each row a result for a different sequence. 
motifRecoveryAUC = function(seq.res){
	# calculate the ranks
	r = t(apply(seq.res, 1, rank))
	
	# previous point needed to calculate AUC-ROC
	auc = rep(0, ncol(r))
	names(auc) = colnames(r)
	
	# do it for all whole steps
	rec = t(sapply(1:ncol(r), function(i) colSums(r<=i)))
	
	list("rec"=rec, "auc"=colSums(rec))
}

#' Calculate PR-AUC for motifs ranked according to some scoring scheme
#'
#' Note that this function asssumes that smaller values are better!
#' 
#' @param seq.res a matrix where each column represents a PWM and each row a result for a different sequence. 
motifPrAUC = function(seq.res){
	# calculate the ranks
	r = t(apply(seq.res, 1, rank))
	
	# previous point needed to calculate AUC-ROC
	auc = rep(0, ncol(r))
	names(auc) = colnames(r)
	
	# do it for all whole steps
	rec = t(sapply(1:ncol(r), function(i) colSums(r<=i)))
	
	recall = rec / nrow(r)
	prec = recall / 1:nrow(rec)
	
	list("prec"=prec, "recall"=recall)
}



