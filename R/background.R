# Functions to create background distributions

#' check consistency of bg.seq input parameter
#'
#' @param bg.seq a set of background sequences, either a list of DNAString object or DNAStringSet object
.normalize.bg.seq = function(bg.seq){
	# check if the sequences are in the right format
	if(class(bg.seq) != "DNAStringSet"){
		if(class(bg.seq) == "list"){
			if(class(bg.seq[[1]]) != "DNAString")
				stop("bg.seq needs to be a list of DNAString objects or a DNAStringSet object")
		} else {
			stop("bg.seq needs to be a list of DNAString objects or a DNAStringSet object")
		}
	} #else {
		#bg.seq = DNAStringSetToList(bg.seq)
	#}
	
	return(bg.seq)
}

#' Make priors from background sequences
#'
#' These priors serve both as background nucleotide frequencies and pseudo-counts
#' for PWMs. 
#'
#' @param bg.seq a set of background sequences
#' @param bg.pseudo.count the total pseudocount shared between nucleotides
#' @export
#' @examples
#' # some example sequences
#' sequences = list(DNAString("AAAGAGAGTGACCGATGAC"), DNAString("ACGATGAGGATGAC"))
#' # make priors with pseudo-count of 1 shared between them
#' makePriors(sequences, 1)
makePriors = function(bg.seq, bg.pseudo.count){
	bg.seq = .normalize.bg.seq(bg.seq)

	# start with the count of the first sequence
	acgt.count = alphabetFrequency(bg.seq[[1]])[c("A", "C", "G", "T")]
	acgt.count = acgt.count[c("A", "C", "G", "T")] + acgt.count[c("T", "G", "C", "A")]
	# add ACGT counts for the rest
	if(length(bg.seq) > 1 ){
		for(i in 2:length(bg.seq)){
			cont = alphabetFrequency(bg.seq[[i]])[c("A", "C", "G", "T")]
			# add the counts on the other strand
			cont[c("A", "C", "G", "T")] = cont[c("A", "C", "G", "T")] + cont[c("T", "G", "C", "A")]
			# sum up
			acgt.count = acgt.count + cont
		}
	}
	
	# scale to bg.pseudo.count
	prior = (acgt.count / sum(acgt.count)) * bg.pseudo.count

	return(prior)	
}


#' Make a lognormal background distribution
#'
#' Construct a lognormal background distribution for a set of sequences. 
#' Sequences concatenated are binned in 'bg.len' chunks and lognormal distribution 
#' fitted to them. 
#'
#' @param bg.seq a set of background sequences, either a list of DNAString object or DNAStringSet object
#' @param motifs a set of motifs, either a list of frequency matrices, or a list of PWM objects. If
#'               frequency matrices are given, the background distribution is fitted from bg.seq. 
#' @param bg.pseudo.count the pseudo count which is shared between nucleotides when frequency matrices are given
#' @param bg.len the length of background chunks
#' @param bg.source a free-form textual description of how the background was generated
#' @param verbose if to produce verbose output
#' @export
#' @examples
#' \dontrun{
#' if(require("PWMEnrich.Dmelanogaster.background")){
#'    data(MotifDb.Dmel.PFM)
#'
#'    # make background for MotifDb motifs using 2kb promoters of all D. melanogaster transcripts 
#' 	  if(require("BSgenome.Dmelanogaster.UCSC.dm3")) 
#'      makePWMLognBackground(Dmelanogaster$upstream2000, MotifDb.Dmel.PFM)
#' }
#' }
makePWMLognBackground = function(bg.seq, motifs, bg.pseudo.count=1, bg.len=1000, bg.source="", verbose=TRUE){
	# check if the sequences are in the right format
	bg.seq = .normalize.bg.seq(bg.seq)

	# convert to list if a single motif is given
	if(!is.list(motifs))
		motifs = list(motifs)

	# concatenate all the background sequences into a single long sequence
	bg.seq.all = concatenateSequences(bg.seq)
	# generate start and end positions
	bg.seq.start = seq(1, nchar(bg.seq.all)+1, bg.len)
	bg.seq.end = bg.seq.start - 1
	bg.seq.start = bg.seq.start[1:(length(bg.seq.start)-1)]
	bg.seq.end = bg.seq.end[2:length(bg.seq.end)]
	
	# split into bg.len sequences
	bg = DNAStringSet(bg.seq.all, start=bg.seq.start, end=bg.seq.end)
	
	if(length(bg)<1000)
		warnings(paste("The number of chunks of size", bg.len, "is smaller than 1000. This might lead to an non-robust estimate of lognormal distribution parameters."))
	
	# convert motifs to PWM format if neccessary
	if(class(motifs[[1]]) != "PWM"){
		prior = makePriors(list(DNAString(bg.seq.all)), bg.pseudo.count)
		pwms = PFMtoPWM(motifs, prior.params = prior)
	} else {
		pwms = motifs
	}
	
	# finally, do motif scanning and calculate mean and sd
	bg.res = motifScores(bg, pwms, verbose=verbose)
	bg.mean = colMeans(bg.res)
	bg.sd = apply(bg.res, 2, sd)
	
	pwm.len = sapply(pwms, length)
	bg.len.real = bg.len - pwm.len + 1
	
	new("PWMLognBackground", bg.source=bg.source, bg.len=bg.len.real, bg.mean=bg.mean, bg.sd=bg.sd, pwms=pwms)
}


#' Make a cutoff background
#'
#' Make a background based on number of motifs hits above a certain threshold. 
#'
#' @param bg.seq a set of background sequences, either a list of DNAString object or DNAStringSet object
#' @param motifs a set of motifs, either a list of frequency matrices, or a list of PWM objects. If
#'               frequency matrices are given, the background distribution is fitted from bg.seq. 
#' @param cutoff the cutoff at which the background should be made, i.e. at which a motif hit is called significant
#' @param bg.pseudo.count the pseudo count which is shared between nucleotides when frequency matrices are given
#' @param bg.source a free-form textual description of how the background was generated
#' @param verbose if to produce verbose output
#' @export
#' @examples
#' \dontrun{
#' if(require("PWMEnrich.Dmelanogaster.background")){
#'    data(MotifDb.Dmel.PFM)
#'
#'    # make background for MotifDb motifs using 2kb promoters of all D. melanogaster transcripts using cutoff of 5
#' 	  if(require("BSgenome.Dmelanogaster.UCSC.dm3")) 
#'      makePWMCutoffBackground(Dmelanogaster$upstream2000, MotifDb.Dmel.PFM, cutoff=log2(exp(5)))
#' }
#' }
makePWMCutoffBackground = function(bg.seq, motifs, cutoff=log2(exp(4)), bg.pseudo.count=1, bg.source="", verbose=TRUE){
	# check if the sequences are in the right format
	bg.seq = .normalize.bg.seq(bg.seq)
	
	# convert to list if a single motif is given
	if(!is.list(motifs))
		motifs = list(motifs)
	
	# make priors and PWMs
	if(class(motifs[[1]]) != "PWM"){
		prior = makePriors(bg.seq, bg.pseudo.count)
		pwms = PFMtoPWM(motifs, prior.params = prior)
	} else {
		pwms = motifs
	}
	
	# scan with cutoff - this gives motif *hit counts*
	bg.res = motifScores(bg.seq, pwms, verbose=verbose, cutoff=cutoff)
	
	seq.len = sapply(bg.seq, length)
	pwm.len = sapply(pwms, length)
	
	# total length of the scanned sequences for each sequence
	total.len = sum(seq.len) - length(seq.len)*(pwm.len-1)
			
	# total number of hits
	total.hits = colSums(bg.res)
	
	# density of binding sites
	bg.P = total.hits / total.len
	
	new("PWMCutoffBackground", bg.source=bg.source, bg.cutoff=cutoff, bg.P=bg.P, pwms=pwms)
}

#' Make an empirical P-value background
#'
#' Make a background appropriate for empirical P-value calculation. The provided set of background
#' sequences is contcatenated into a single long sequence which is then scanned with the motifs
#' and raw scores are saved. This object can be very large. 
#'
#' For reliable P-value calculation the size of the background set needs to be at least seq.len / min.P.value.
#' For instance, to get P-values at a resolution of 0.001 for a single sequence of 500bp, we would need 
#' a background of at least 500/0.001 = 50kb. This ensures that we can make 1000 independent 500bp samples from
#' this background to properly estimate the P-value. For a group of sequences, we would take seq.len to be the
#' total length of all sequences in a group. 
#'
#' @param bg.seq a set of background sequences, either a list of DNAString object or DNAStringSet object
#' @param motifs a set of motifs, either a list of frequency matrices, or a list of PWM objects. If
#'               frequency matrices are given, the background distribution is fitted from bg.seq. 
#' @param bg.pseudo.count the pseudo count which is shared between nucleotides when frequency matrices are given
#' @param bg.source a free-form textual description of how the background was generated
#' @param verbose if to produce verbose output
#' @param ... currently unused (this is for convenience for makeBackground function)
#' @export
#' @examples
#' \dontrun{
#' if(require("PWMEnrich.Dmelanogaster.background")){
#'    data(MotifDb.Dmel.PFM)
#'
#'    # make empirical background by saving raw scores for each bp in the sequence - this can be very large in memory!
#' 	  if(require("BSgenome.Dmelanogaster.UCSC.dm3")) 
#'      makePWMEmpiricalBackground(Dmelanogaster$upstream2000[1:100], MotifDb.Dmel.PFM)
#' }
#' }
makePWMEmpiricalBackground = function(bg.seq, motifs, bg.pseudo.count=1, bg.source="", verbose=TRUE, ...){
	# check if the sequences are in the right format
	bg.seq = .normalize.bg.seq(bg.seq)
	
	# convert to list if a single motif is given
	if(!is.list(motifs))
		motifs = list(motifs)
	
	# make priors and PWMs
	if(class(motifs[[1]]) != "PWM"){
		prior = makePriors(bg.seq, bg.pseudo.count)
		pwms = PFMtoPWM(motifs, prior.params = prior)
	} else {
		pwms = motifs
	}
	
	# scan one very long sequence that is formed by concatenating all sequences
	bg.seq.all = DNAString(concatenateSequences(bg.seq))
	
	bg.res = motifScores(bg.seq.all, pwms, raw.scores=TRUE, verbose=verbose)
	
	bg.len = length(bg.seq.all)
	bg.fwd = bg.res[[1]][1:bg.len,]
	bg.rev = bg.res[[1]][(bg.len+1):(bg.len*2),]
	
	if(length(motifs) == 1){
		bg.fwd = matrix(bg.fwd, ncol=1, dimnames=list(NULL, names(motifs)))
		bg.rev = matrix(bg.rev, ncol=1, dimnames=list(NULL, names(motifs)))
	}
	
	new("PWMEmpiricalBackground", bg.source=bg.source, bg.fwd=bg.fwd, bg.rev=bg.rev, pwms=pwms)
}

#' Construct a cutoff background from empirical background
#'
#' This function takes already calculated empirical background distribution and chooses
#' cutoff for each motif based on P-value cutoff for individual sites. 
#'
#' @param bg.p an object of class PWMEmpiricalBackground
#' @param p.value the P-value used to find cuttoffs for each of the motifs
#' @param bg.source textual description of background source
#' @return an object of type PWMCutoffBackground
#' @export
#' @examples
#' \dontrun{
#' if(require("PWMEnrich.Dmelanogaster.background")){
#'    data(MotifDb.Dmel.PFM)
#'
#'    # make empirical background - here we use only 100 sequences for illustrative purposes
#' 	  if(require("BSgenome.Dmelanogaster.UCSC.dm3")) 
#'      bg.p = makePWMEmpiricalBackground(Dmelanogaster$upstream2000[1:100], MotifDb.Dmel.PFM)
#'
#'    # use the empirical background to pick a threshold and make cutoff background
#'    makePWMPvalCutoffBackground(bg.p, 0.001)
#' }
#' }
makePWMPvalCutoffBackground = function(bg.p, p.value=1e-3, bg.source=""){
	# extract stuffs
	bg.fwd = bg.p@bg.fwd
	bg.rev = bg.p@bg.rev
	pwms = bg.p@pwms
	
	bg.len = nrow(bg.fwd)
	
	if(bg.len * p.value < 1){
		stop("The P-value is too small for background of this size. Should have at least 1/p.value nucleotides in the background.")
	}
	
	# get cutoffs
	cutoff = structure(rep(0, length(pwms)), names=names(pwms))
	P = structure(rep(0, length(pwms)), names=names(pwms))
	for(i in 1:length(pwms)){
		# join both strands
		all.data = log2(c(bg.fwd[!is.na(bg.fwd[,i]),i], bg.rev[!is.na(bg.rev[,i]),i]))
		
		# stop if not all data is finite, because then there is an error!
		stopifnot(length(is.finite(all.data)) == length(all.data))
		
		# find out the cutoff for 1-p.value
		cutoff[i] = quantile(ecdf(all.data), 1 - p.value)
		# work out the actual P value
		# NOTE: since the distribution is decrete, this won't be exactly the same as 2*P.value but will be quite close
		P[i] = 2 * (sum(all.data >= cutoff[i]) / length(all.data))
	}

	# final object
	new("PWMCutoffBackground", bg.source=bg.source, bg.cutoff=cutoff, bg.P=P, pwms=pwms)
}

#' Construct a P-value cutoff background from a set of sequences
#'
#' This function creates a P-value cutoff background for motif enrichment. 
#'
#' @param bg.seq a set of background sequences, either a list of DNAString object or DNAStringSet object
#' @param motifs a set of motifs, either a list of frequency matrices, or a list of PWM objects. If
#'               frequency matrices are given, the background distribution is fitted from bg.seq. 
#' @param p.value the P-value used to find cuttoffs for each of the motifs
#' @param bg.pseudo.count the pseudo count which is shared between nucleotides when frequency matrices are given
#' @param bg.source textual description of background source
#' @param verbose if to print verbose output
#' @return an object of type PWMCutoffBackground
#' @export
#' @examples
#' \dontrun{
#' if(require("PWMEnrich.Dmelanogaster.background")){
#'    data(MotifDb.Dmel.PFM)
#'
#'    # use the empirical background to pick a threshold and make cutoff background
#'    makePWMPvalCutoffBackground(Dmelanogaster$upstream2000, 0.001)
#' }
#' }
makePWMPvalCutoffBackgroundFromSeq = function(bg.seq, motifs, p.value=1e-3, bg.pseudo.count=1, bg.source="", verbose=TRUE){
	# check if the sequences are in the right format
	bg.seq = .normalize.bg.seq(bg.seq)
	
	# convert to list if a single motif is given
	if(!is.list(motifs))
		motifs = list(motifs)
	
	# make priors and PWMs
	if(class(motifs[[1]]) != "PWM"){
		prior = makePriors(bg.seq, bg.pseudo.count)
		pwms = PFMtoPWM(motifs, prior.params = prior)
	} else {
		pwms = motifs
	}
	
	# pre-calculate the big sequence motifScoresBigMemory() is faster
	seq.input = .inputParamSequences(bg.seq)
	seq.all = DNAString(concatenateSequences(seq.input))
	
	# get cutoffs
	cutoff = structure(rep(0, length(pwms)), names=names(pwms))
	P = structure(rep(0, length(pwms)), names=names(pwms))
	for(i in 1:length(pwms)){	
		# get the raw values for this motif only!
		bg.res = motifScoresBigMemory(bg.seq, pwms[i], verbose=verbose, raw.scores=TRUE, seq.all=seq.all)
		bg.res = na.omit(unlist(bg.res))
		
		# transform to log2 to get the final data
		all.data = log2(bg.res)
		
		# find out the cutoff for 1-p.value
		cutoff[i] = quantile(ecdf(all.data), 1 - p.value)
		# work out the actual P value
		# NOTE: since the distribution is decrete, this won't be exactly the same as 2*P.value but will be quite close
		P[i] = 2 * (sum(all.data >= cutoff[i]) / length(all.data))
	}

	# final object
	new("PWMCutoffBackground", bg.source=bg.source, bg.cutoff=cutoff, bg.P=P, pwms=pwms)
}

#' Divide total.len into fragments of length len by providing start,end positions
#'
#' @param total.len total available length to be subdivided
#' @param len size of the individual chunk
#' @return a data.frame containing paired up start,end positions
makeStartEndPos = function(total.len, len){
	# make start positions
	start = seq(1, total.len+1, len)
	# make end positions
	end = start - 1
	# pair them up
	start = start[1:(length(start)-1)]
	end = end[2:length(end)]
	
	data.frame(start, end)
}

#' Make a GEV background distribution
#'
#' Construct a lognormal background distribution for a set of sequences. 
#' Sequences concatenated are binned in 'bg.len' chunks and lognormal distribution 
#' fitted to them. 
#'
#' @param bg.seq a set of background sequences, either a list of DNAString object or DNAStringSet object
#' @param motifs a set of motifs, either a list of frequency matrices, or a list of PWM objects. If
#'               frequency matrices are given, the background distribution is fitted from bg.seq. 
#' @param bg.pseudo.count the pseudo count which is shared between nucleotides when frequency matrices are given
#' @param bg.len the length range of background chunks
#' @param bg.source a free-form textual description of how the background was generated
#' @param verbose if to produce verbose output
#' @param fit.log if to fit log odds (instead of odds)
#' @export
#' @examples
#' \dontrun{
#' if(require("PWMEnrich.Dmelanogaster.background")){
#'    data(MotifDb.Dmel.PFM)
#'
#'    # make background for MotifDb motifs using 2kb promoters of all D. melanogaster transcripts 
#' 	  if(require("BSgenome.Dmelanogaster.UCSC.dm3")) 
#'      makePWMGEVBackground(Dmelanogaster$upstream2000, MotifDb.Dmel.PFM)
#' }
#' }
makePWMGEVBackground = function(bg.seq, motifs, bg.pseudo.count=1, bg.len=seq(200,2000,200), bg.source="", verbose=TRUE, fit.log=TRUE){
	# check if the sequences are in the right format
	bg.seq = .normalize.bg.seq(bg.seq)

	# convert to list if a single motif is given
	if(!is.list(motifs))
		motifs = list(motifs)
		
	# give automatic names
	if(is.null(names(motifs))){
		names(motifs) = paste("Motif", 1:length(motifs), sep="")
	}


	# concatenate all the background sequences into a single long sequence
	bg.seq.all = concatenateSequences(bg.seq)
	
	# convert motifs to PWM format if neccessary
	if(class(motifs[[1]]) != "PWM"){
		prior = makePriors(list(DNAString(bg.seq.all)), bg.pseudo.count)
		pwms = PFMtoPWM(motifs, prior.params = prior)
	} else {
		pwms = motifs
	}
		
	# GEV parameters that need to be fitted in linear regression
	params = array(0, dim=c(length(bg.len), 3, length(pwms)), 
		dimnames=list(bg.len, c("loc", "scale", "shape"), names(pwms)))
	
	# iterate over background lengths and fit GEV distributions
	for(i in 1:length(bg.len)){
		if(verbose){
			message("Scanning background of size ", bg.len[i])
		}
		bg.pos = makeStartEndPos(nchar(bg.seq.all), bg.len[i])
		bg = DNAStringSet(bg.seq.all, start=bg.pos$start, end=bg.pos$end)
		
		# scan
		bg.res = motifScores(bg, pwms, verbose=verbose)
		
		# fit GEV for each PWM
		params[i,,] = apply(bg.res, 2, function(x) {
			if(fit.log)
				x = log(x)
			fgev(x, std.err=FALSE)$param[c("loc", "scale", "shape")]
		})
	}
	
	## now fit linear regression for each PWM
	
	bg.loc = lapply(names(pwms), function(pwm.name){
		loc = params[,"loc",pwm.name]
		log.len = log(bg.len)
		lm(loc~log.len)
	})
	names(bg.loc) = names(pwms)
	
	bg.scale = lapply(names(pwms), function(pwm.name){
		scale = params[,"scale",pwm.name]
		log.len = log(bg.len)
		lm(scale~log.len)
	})
	names(bg.scale) = names(pwms)
	
	bg.shape = lapply(names(pwms), function(pwm.name){
		shape = params[,"shape",pwm.name]
		log.len = log(bg.len)
		lm(shape~log.len)
	})
	names(bg.shape) = names(pwms)
	
	# return the object	
	new("PWMGEVBackground", bg.source=bg.source, bg.loc=bg.loc, bg.scale=bg.scale, bg.shape=bg.shape, pwms=pwms)
}

#' A helper function to pick a genome for an organism
#' 
#' @param organism either organism name (such as "dm3") or a BSgenome object
#' @return a BSgenome object
pickGenome = function(organism){
	# pick the genome object based on the 
	if(is(organism, "BSgenome")){ 
		genome = organism
	} else if(organism == "dm3"){
		if(!require("BSgenome.Dmelanogaster.UCSC.dm3"))
			stop("This functions requires the BSgenome.Dmelanogaster.UCSC.dm3 package to build background for D. melanogaster")
		genome = Dmelanogaster
	} else {
		stop("Please pick one of the valid organisms: \"dm3\" or provide a BSgenome object of the target genome.")
	}
	
	genome
}

#' A helper function to select a set of non-redundant promoters from a genome object
#'
#' @param genome a BSgenome object
#' @return a vector of inidicies of non-redundant promoters
selectPromoters = function(genome){
	promoters = genome$upstream2000
	if(organism(genome) == "Drosophila melanogaster"){
		promoter.genes = sapply(strsplit(names(promoters), "-"), function(x) x[1])
		select.one = tapply(1:length(promoter.genes), promoter.genes, function(x) x[1])
	} else{
		select.one = 1:length(promoters)
	}
	
	select.one
}

#' Make a background for a set of position frequency matrices
#'
#' This is a convenience front-end function to compile new backgrounds for a set of PFMs. 
#' Currently only supports D. melanogaster, but in the future should support other common organisms as well. 
#'
#' @param motifs a list of position frequency matrices (4xL matrices)
#' @param organism either a name of the organisms for which the background should be compiled 
#'                 (currently only supported name is "dm3" for Drosophila Melanogaster), or a \code{BSgenome} object (see \code{BSgenome} package). 
#' @param type the type of background to be compiled. Possible types are: 
#'             \itemize{
#'                 \item "logn" - estimate a lognormal background
#'                 \item "cutoff" - estimate a Z-score background with fixed log-odds cutoff (in log2)
#'                 \item "pval" - estimate a Z-score background with a fixed P-value cutoff. Note that this may require a lot of memory
#'                                since the P-value of motif hits is first estimated from the empirical distribution. 
#'                 \item "empirical" - create an empirical P-value background. Note that this may require a lot of memory (up to 10GB in
#'                                     default "slow" mode (quick=FALSE) for 126 JASPAR motifs and 1000 D. melanogaster promoters).
#'                 \item "GEV" - estimate a generalized extreme value (GEV) distribution background by fitting linear regression to distribution 
#'                               parameters in log space
#'             }
#' @param quick if to preform fitting on a reduced set of 100 promoters. This will not give as good results but is much quicker than fitting to all the promoters (~10k). 
#'              Usage of this parameter is recommended only for testing and rough estimates. 
#' @param ... other named parameters that backend function makePWM***Background functions take.
#' @export
#' @author Robert Stojnic, Diego Diez
#' @examples
#' 
#' # load in the two example de-novo motifs
#' motifs = readMotifs(system.file(package="PWMEnrich", dir="extdata", file="example.transfac"), remove.acc=TRUE)
#' 
#' \dontrun{
#'   # construct lognormal background
#'   bg.logn = makeBackground(motifs, organism="dm3", type="logn")
#'
#'   # alternatively, any BSgenome object can also be used
#'   if(require("BSgenome.Dmelanogaster.UCSC.dm3"))
#'     bg.logn = makeBackground(motifs, organism=Dmelanogaster, type="logn")
#'
#'   # construct a Z-score of hits with P-value background
#'   bg.pval = makeBackground(motifs, organism="dm3", type="pval", p.value=1e-3)
#'
#'   # now we can use them to scan for enrichment in sequences (in this case there is a consensus Tin binding site)
#'   motifEnrichment(DNAString("TGCATCAAGTGTGTAGTG"), bg.logn)
#'   motifEnrichment(DNAString("TGCATCAAGTGTGTAGTG"), bg.pval)
#' }
#' 
makeBackground = function(motifs, organism="dm3", type="logn", quick=FALSE, ...){
	# check input parameters
	valid.types = c("logn", "cutoff", "pval", "empirical", "GEV")
	if(!(type %in% valid.types)){
		stop(paste("Invalid type, please choose from:", paste(valid.types, collapse=", ")))
	}

	if(!exists("bg.source"))
		bg.source = NULL

	genome = pickGenome(organism)
	
	# take only a single promoter from each of the genes (this works only for Drosophila)
	select.one = selectPromoters(genome)
	
	if(quick){
		bg.seq = genome$upstream2000[select.one[seq(1, length(select.one), length.out=100)]]
		if(is.null(bg.source))
			bg.source = paste(organism(genome), " (", providerVersion(genome), ") 100 unique 2kb promoters", sep="")
	} else {
		if(type %in% c("pval")){
			bg.seq = genome$upstream2000[select.one[seq(1, length(select.one), length.out=500)]]
			if(is.null(bg.source))
				bg.source = paste(organism(genome), " (", providerVersion(genome), ") 500 unique 2kb promoters", sep="")
		} else {
			bg.seq = genome$upstream2000[select.one]
			if(is.null(bg.source))
				bg.source = paste(organism(genome), " (", providerVersion(genome), ") ", length(select.one), " unique 2kb promoters", sep="")
		}
	}
	
	#bg.seq = DNAStringSetToList(bg.seq)
	
	## now run the appropriate backend function
	if(type == "logn"){
		bg = makePWMLognBackground(bg.seq, motifs, bg.source=bg.source, ...)
	} else if(type == "cutoff"){
		bg = makePWMCutoffBackground(bg.seq, motifs, bg.source=bg.source, ...)
	} else if(type == "pval"){
		bg = makePWMPvalCutoffBackgroundFromSeq(bg.seq, motifs, bg.source=bg.source, ...)
	} else if(type == "empirical"){
		bg = makePWMEmpiricalBackground(bg.seq, motifs, bg.source=bg.source, ...)
	} else if(type == "GEV"){
		bg = makePWMGEVBackground(bg.seq, motifs, bg.source=bg.source, ...)
	}
	
	return(bg)
}

#' Get the four nucleotides background frequencies
#'
#' Estimate the background frequencies of A,C,G,T on a set of promoters from an organism
#'
#' @param organism either a name of the organisms for which the background should be compiled 
#'                 (currently only supported name is "dm3" for Drosophila Melanogaster), or a \code{BSgenome} object (see \code{BSgenome} package).
#' @param pseudo.count the number to which the frequencies sum up to, by default 1
#' @param quick if to preform fitting on a reduced set of 100 promoters. This will not give as good results but is much quicker than fitting to all the promoters (~10k). 
#'              Usage of this parameter is recommended only for testing and rough estimates.
#' @export
#' @author Robert Stojnic, Diego Diez
#' @examples
#' \dontrun{
#'   getBackgroundFrequencies("dm3")
#' }
getBackgroundFrequencies = function(organism="dm3", pseudo.count=1, quick=FALSE){
	# pick the set of background sequences
	genome = pickGenome(organism)

	# take only a single promoter from each of the genes
	select.one = selectPromoters(genome)
		
	if(quick){
		bg.seq = genome$upstream2000[select.one[seq(1, length(select.one), length.out=100)]]
	} else {
		bg.seq = genome$upstream2000[select.one]
	}
	
	#bg.seq = DNAStringSetToList(bg.seq)
	
	makePriors(bg.seq, pseudo.count)
}

#' Calculate the empirical distribution score distribution for a set of motifs
#'
#' @param motifs a set of motifs, either a list of frequency matrices, or a list of PWM objects. If
#'               frequency matrices are given, the background distribution is fitted from bg.seq. 
#' @param organism either a name of the organisms for which the background should be compiled 
#'                 (currently only supported name is "dm3" for Drosophila Melanogaster), or a \code{BSgenome} object (see \code{BSgenome} package).
#' @param bg.seq a set of background sequence (either this or organism needs to be specified!). Can be a DNAString or DNAStringSet object.
#' @param quick if to do the fitting only on a small subset of the data (only in combination with \code{organism}). Useful only for code testing!
#' @param pseudo.count the pseudo count which is shared between nucleotides when frequency matrices are given
#' @return a list of \code{ecdf} objects (see help page for \code{ecdf} for usage). 
#' @export
motifEcdf = function(motifs, organism=NULL, bg.seq=NULL, quick=FALSE, pseudo.count=1){
	if(is.null(organism) && is.null(bg.seq)){
		stop("Either the 'organism' or 'bg.seq' parameter need to be specified!")
	}
	
	if(!is.list(motifs))
		motifs = list(motifs)
	
	if(!is.null(bg.seq)){
		# check if the sequences are in the right format
		bg.seq = .normalize.bg.seq(bg.seq)
	} else {
		# pick the set of background sequences
		genome = pickGenome(organism)

		# take only a single promoter from each of the genes
		select.one = selectPromoters(genome)
	
		if(quick){
			bg.seq = genome$upstream2000[select.one[seq(1, length(select.one), length.out=100)]]
		} else {
			bg.seq = genome$upstream2000[select.one]
		}
	}
	
	# make priors and PWMs
	if(class(motifs[[1]]) != "PWM"){
		prior = makePriors(bg.seq, pseudo.count)
		pwms = PFMtoPWM(motifs, prior.params = prior)
	} else {
		pwms = motifs
	}	
	
	# always use the big memory implementation as it is much faster in typical usage!
	s = motifScoresBigMemory(bg.seq, pwms, raw.scores=TRUE)
	
	# group together distributions for motifs
	s = lapply(1:length(pwms), function(i) na.omit(unlist(lapply(s, function(x) x[,i]))))
	
	# create empirical CDFs and return
	e = lapply(s, ecdf)
	names(e) = names(pwms)
	
	e
}
