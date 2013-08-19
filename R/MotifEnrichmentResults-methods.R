
#' Name of different pieces of information associated with MotifEnrichmentResults
#'
#' @title Names of variables
#' @name names,MotifEnrichmentResults
#' @aliases names,MotifEnrichmentResults-method
#' @param x the MotifEnrichmentResults object
#' @return the names of the variables
#' @rdname operators-MotifEnrichmentResults
setMethod("names", signature=signature(x="MotifEnrichmentResults"), function(x) names(x@res))

#' Access a property by name
#'
#' @aliases $,MotifEnrichmentResults-method
#' @param x the MotifEnrichmentResults object
#' @param name the variable name
#' @rdname operators-MotifEnrichmentResults
setMethod("$", signature=signature(x="MotifEnrichmentResults"), function(x, name){
	x@res[[name]]
})


#' show method for MotifEnrichmentResults
#' @param object the MotifEnrichmentResults object
setMethod("show", signature=signature(object="MotifEnrichmentResults"), function(object){
	res = object@res
	is.group.only = length(grep("sequence[.]", names(res))) == 0
	
	cat("An object of class 'MotifEnrichmentResults':\n")
	cat("* created with '", res$score, "' scoring function with '", res$bg, "' background correction\n", sep="")
	cat("* on a set of ", length(res$sequences), " sequence(s) and ", length(res$pwms), " PWMs\n", sep="")
	cat("Result sets for the group:",paste("$", names(res)[grep("group[.]", names(res))], sep="", collapse=", "), "\n")
	if(!is.group.only){
		cat("Result sets for individual sequences:", paste("$", names(res)[grep("sequence[.]", names(res))], sep="", collapse=", "), "\n")
		cat("Report methods: groupReport(), sequenceReport()\n")
		#cat("Methods to extract data: motifRankingForGroup(), motifRankingForSequence()\n")
		#cat("Methods to plot data: plotTopMotifsGroup(), plotTopMotifsSequence()\n")
	} else {
		cat("Report method: sequenceReport()\n")
		#cat("Methods to extract data: motifRankingForGroup()\n")
		#cat("Methods to plot data: plotTopMotifsGroup()\n")
	}
	
})

#' A helper function for motifRankingForGroup and motifRankingForSequence with the common code
#'
#' @param res the list of results from MotifEnrichmentResults object
#' @param r the vector of raw results that needs to be processed
#' @param id if to return IDs instead of names
#' @param order if to return the ordering of motifs
#' @param rank if to return the rank of motifs
#' @param unique if to remove duplicates
#' @param decreasing specifies the sorting order
rankingProcessAndReturn = function(res, r, id, order, rank, unique, decreasing){
	if((unique & id) || (unique & order)){
		stop("Parameter 'unique' can be set to TRUE only if 'id' and 'order' are FALSE.")
	}


	# set either the names or IDs
	all.ids = sapply(res$pwms, function(x) x@id)
	all.names = sapply(res$pwms, function(x) x@name)
	
	if(id){
		names(r) = all.ids
	} else {
		if(!all(all.names == ""))
			names(r) = all.names
	}
	
	# do the required transformation on the results
	ret = NULL
	if(rank){
		x = base::rank(r)
		if(decreasing)
			x = length(x) - x + 1
		names(x) = names(r)
		ret = x
	} else if(order){
		x = base::order(r, decreasing=decreasing)
		names(x) = names(r)[x]
		ret = x
	} else {
		ret = sort(r, decreasing=decreasing)
	}
	
	# remove duplicates if applicable	
	if(unique){
		original.names = names(ret)		
		if(decreasing & !rank)
			ret = tapply(ret, names(ret), max)
		else
			ret = tapply(ret, names(ret), min)
			
		if(rank){
			# re-rank the ranks... 
			ret = rank(ret) 
		}
		
		# unique names with correct ordering
		uniq.names = c()
		for(i in 1:length(original.names)){
			if(!(original.names[i] %in% uniq.names))
				uniq.names[length(uniq.names)+1] = original.names[i]
		}
		
		# re-order in the original ordering	
		ret = ret[uniq.names]
	}
	
	ret

}

#' Get a ranking of motifs by their enrichment in the whole set of sequences
#'
#' @name motifRankingForGroup,MotifEnrichmentResults-method
#' @aliases motifRankingForGroup
#' @param obj a MotifEnrichmentResults object
#' @param bg if to use background corrected P-values to do the ranking (if available)
#' @param id if to show PWM IDs instead of target TF names
#' @param order if to output the ordering of PWMs instead of actual P-values or raw values
#' @param rank if the output should be rank of a PWM instead of actual P-values or raw values
#' @param unique if TRUE, only the best rank is taken for each TF (only when id = FALSE, order = FALSE)
#' @param ... currently unused
#'
#' @return a vector of P-values or raw enrichments sorted such that the first motif is most enriched
#' @export
#' @examples
#' if(require("PWMEnrich.Dmelanogaster.background")){
#'    ###
#'    # load the pre-compiled lognormal background
#'    data(PWMLogn.dm3.MotifDb.Dmel)
#'
#'    # scan two sequences for motif enrichment
#'    sequences = list(DNAString("GAAGTATCAAGTGACCAGTAAGTCCCAGATGA"), DNAString("AGGTAGATAGAACAGTAGGCAATGAAGCCGATG"))
#'    res = motifEnrichment(sequences, PWMLogn.dm3.MotifDb.Dmel)
#'
#'    # most enriched in both sequences (sorted by lognormal background P-value)
#'    head(motifRankingForGroup(res))
#'
#'    # Return a non-redundant set of TFs
#'    head(motifRankingForGroup(res, unique=TRUE))
#'
#'    # sorted by raw affinity instead of P-value
#'    head(motifRankingForGroup(res, bg=FALSE))
#'
#'    # show IDs instead of target TF names
#'    head(motifRankingForGroup(res, id=TRUE))
#'
#'    # output the rank instead of P-value
#'    head(motifRankingForGroup(res, rank=TRUE))
#' }
setMethod("motifRankingForGroup", signature=signature(obj="MotifEnrichmentResults"), function(obj, bg=TRUE, id=FALSE, order=FALSE, rank=FALSE, unique=FALSE, ...){
	res = obj@res
	
	# vector of scores
	if(bg && "group.bg" %in% names(res)){
		r = res$group.bg
		if(res$score == "cutoff" | res$bg == "ms")
			decreasing = TRUE
		else
			decreasing = FALSE
	} else {
		#if(bg) 
		#	warning("Parameter 'bg' is TRUE but this MotifEnrichmentResults object has no background correction, ignoring parameter.")
		r = res$group.nobg
		decreasing = TRUE
	}
	
	rankingProcessAndReturn(res, r, id, order, rank, unique, decreasing)
})

#' Get a ranking of motifs by their enrichment in one specific sequence
#'
#' @name motifRankingForSequence,MotifEnrichmentResults-method
#' @aliases motifRankingForSequence
#' @param obj a MotifEnrichmentResults object
#' @param seq.id either the sequence number or sequence name 
#' @param bg if to use background corrected P-values to do the ranking (if available)
#' @param id if to show PWM IDs instead of target TF names
#' @param order if to output the ordering of PWMs instead of actual P-values or raw values
#' @param rank if the output should be rank of a PWM instead of actual P-values or raw values
#' @param unique if TRUE, only the best rank is taken for each TF (only when id = FALSE, order = FALSE)
#' @param ... currently unused
#'
#' @export
#' @return a vector of P-values or raw enrichments sorted such that the first motif is most enriched
#'
#' @examples
#' if(require("PWMEnrich.Dmelanogaster.background")){
#'    ###
#'    # load the pre-compiled lognormal background
#'    data(PWMLogn.dm3.MotifDb.Dmel)
#'
#'    # scan two sequences for motif enrichment
#'    sequences = list(DNAString("GAAGTATCAAGTGACCAGTAAGTCCCAGATGA"), DNAString("AGGTAGATAGAACAGTAGGCAATGAAGCCGATG"))
#'    res = motifEnrichment(sequences, PWMLogn.dm3.MotifDb.Dmel)
#'
#'    # most enriched in the second sequences (sorted by lognormal background P-value)
#'    head(motifRankingForSequence(res, 2))
#'
#'    # return unique TFs enriched in sequence 2
#'    head(motifRankingForSequence(res, 2, unique=TRUE))
#'
#'    # sorted by raw affinity instead of P-value
#'    head(motifRankingForSequence(res, 2, bg=FALSE))
#'
#'    # show IDs instead of target TF names
#'    head(motifRankingForSequence(res, 2, id=TRUE))
#'
#'    # output the rank instead of P-value
#'    head(motifRankingForSequence(res, 2, rank=TRUE))
#' }
setMethod("motifRankingForSequence", signature=signature(obj="MotifEnrichmentResults"), function(obj, seq.id, bg=TRUE, id=FALSE, order=FALSE, rank=FALSE, unique=FALSE, ...){
	res = obj@res
	
	if(missing(seq.id)){
		stop("Please specify the sequence number with 'seq.id'")
	}
	
	# vector of scores
	if(bg && "sequence.bg" %in% names(res)){
		r = res$sequence.bg[seq.id,]
		if(res$score == "cutoff" | res$bg == "ms")
			decreasing = TRUE
		else
			decreasing = FALSE
	} else {
		#if(bg) 
		#	warning("Parameter 'bg' is TRUE but this MotifEnrichmentResults object has no background correction, ignoring parameter.")
		r = res$sequence.nobg[seq.id,]
		decreasing = TRUE
	}
	
	rankingProcessAndReturn(res, r, id, order, rank, unique, decreasing)

})

#' Plot the top N enrichment motifs in a group of sequences
#'
#' @name plotTopMotifsGroup,MotifEnrichmentResults-method
#' @aliases plotTopMotifsGroup
#' @param obj a MotifEnrichmentResults object
#' @param n the number of top ranked motifs to plot
#' @param bg if to use background corrected P-values to do the ranking (if available)
#' @param id if to show PWM IDs instead of target TF names
#' @param ... other parameters passed to \code{plotMultipleMotifs()}
#' @export
#' @examples
#' if(require("PWMEnrich.Dmelanogaster.background")){
#'    ###
#'    # load the pre-compiled lognormal background
#'    data(PWMLogn.dm3.MotifDb.Dmel)
#'
#'    # scan two sequences for motif enrichment
#'    sequences = list(DNAString("GAAGTATCAAGTGACCAGTAAGTCCCAGATGA"), DNAString("AGGTAGATAGAACAGTAGGCAATGAAGCCGATG"))
#'
#'    res = motifEnrichment(sequences, PWMLogn.dm3.MotifDb.Dmel)
#'
#'    # plot the top 4 motifs in a 2x2 grid
#'    plotTopMotifsGroup(res, 4)
#'
#'    # plot top 3 motifs in a single row
#'    plotTopMotifsGroup(res, 3, row=1, cols=3)
#' }
setMethod("plotTopMotifsGroup", signature=signature(obj="MotifEnrichmentResults"), function(obj, n, bg=TRUE, id=FALSE, ...){
	o = motifRankingForGroup(obj, bg, id, order=TRUE)
	
	plotMultipleMotifs(obj@res$pwms[o[1:n]], names(o)[1:n], ...)
})

#' Plot the top N enrichment motifs in a single sequence
#'
#' @name plotTopMotifsSequence,MotifEnrichmentResults-method
#' @aliases plotTopMotifsSequence
#' @param obj a MotifEnrichmentResults object
#' @param seq.id either the sequence number or sequence name 
#' @param n the number of top ranked motifs to plot
#' @param bg if to use background corrected P-values to do the ranking (if available)
#' @param id if to show PWM IDs instead of target TF names
#' @param ... other parameters passed to \code{plotMultipleMotifs()}
#' @export
#' @examples
#' if(require("PWMEnrich.Dmelanogaster.background")){
#'    ###
#'    # load the pre-compiled lognormal background
#'    data(PWMLogn.dm3.MotifDb.Dmel)
#'
#'    # scan two sequences for motif enrichment
#'    sequences = list(DNAString("GAAGTATCAAGTGACCAGTAAGTCCCAGATGA"), DNAString("AGGTAGATAGAACAGTAGGCAATGAAGCCGATG"))
#'
#'    res = motifEnrichment(sequences, PWMLogn.dm3.MotifDb.Dmel)
#'
#'    # plot the top 4 motifs in a 2x2 grid
#'    plotTopMotifsSequence(res, 1, 4)
#'
#'    # plot top 3 motifs in a single row
#'    plotTopMotifsSequence(res, 1, 3, row=1, cols=3)
#' }
setMethod("plotTopMotifsSequence", signature=signature(obj="MotifEnrichmentResults"), function(obj, seq.id, n, bg=TRUE, id=FALSE, ...){
	o = motifRankingForSequence(obj, seq.id, bg, id, order=TRUE)
	
	plotMultipleMotifs(obj@res$pwms[o[1:n]], names(o)[1:n], ...)
})

#' Generate a motif enrichment report for the whole group of sequences together
#'
#' @name groupReport,MotifEnrichmentResults-method
#' @aliases groupReport
#' @param obj a MotifEnrichmentResults object
#' @param top what proportion of top motifs should be examined in each individual sequence (by default 5%)
#' @param bg if to use background corrected P-values to do the ranking (if available)
#' @param ... unused
#' @return a MotifEnrichmentReport object containing a table with the following columns: 
#' \itemize{
#'    \item 'rank' - The rank of the PWM's enrichment in the whole group of sequences together
#'    \item 'target' - The name of the PWM's target gene, transcript or protein complex. 
#'    \item 'id' - The unique identifier of the PWM (if set during PWM creation).
#'    \item 'raw.score' - The raw score before P-value calculation
#'    \item 'p.value' - The P-value of motif enrichment (if available)
#'    \item 'top.motif.prop' - The proportion (between 0 and 1) of sequences where the motif is within \code{top} proportion of enrichment motifs. 
#' } 
#' @export
#' @examples
#' if(require("PWMEnrich.Dmelanogaster.background")){
#'    ###
#'    # load the pre-compiled lognormal background
#'    data(PWMLogn.dm3.MotifDb.Dmel)
#'
#'    # scan two sequences for motif enrichment
#'    sequences = list(DNAString("GAAGTATCAAGTGACCAGTAAGTCCCAGATGA"), DNAString("AGGTAGATAGAACAGTAGGCAATGAAGCCGATG"))
#'
#'    res = motifEnrichment(sequences, PWMLogn.dm3.MotifDb.Dmel)
#'
#'    # produce a report for all sequences taken together
#'    r.default = groupReport(res)
#'
#'    # produce a report where the last column takes top 1% motifs
#'    r = groupReport(res, top=0.01)
#'
#'    # view the results
#'    r
#'
#'    # plot the top 10 most enriched motifs
#'    plot(r[1:10])
#' 
#' }
setMethod("groupReport", signature=signature(obj="MotifEnrichmentResults"), function(obj, top=0.05, bg=TRUE, ...){
	pwms = obj$pwms	
		
	# correct ordering of motifs
	o = motifRankingForGroup(obj, rank=TRUE, bg=bg)
	
	targets = sapply(obj$pwms, function(x) x@name)
	ids = sapply(obj$pwms, function(x) x@id)
	
	# construct the initial data frame object
	df = data.frame(rank=o, target=targets, id=ids, raw.score=obj$group.nobg, stringsAsFactors=FALSE)
	
	if(bg && "group.bg" %in% names(obj)){
		p.value = obj$group.bg
	} else {
		p.value = as.numeric(NA)
	}

	df = data.frame(df, p.value=p.value, stringsAsFactors=FALSE)
	
	# calculate in how many top sequences is the motif present
	if("sequence.nobg" %in% names(obj)){
		pct.top = rep(0, nrow(df))
		pct.cutoff = nrow(df) * top
		num.seq = nrow(obj$sequence.nobg)
		for(i in 1:num.seq){
			os = motifRankingForSequence(obj, i, rank=TRUE, bg=bg)
			pct.top = pct.top + (os <= pct.cutoff)
		}
		pct.top = pct.top / num.seq
		
		top.motif = pct.top
	} else {
		top.motif = as.numeric(NA)
	}
	
	df = data.frame(df, "top.motif.prop"=top.motif, stringsAsFactors=FALSE)
	
	# sort and return
	correct.order = order(df$rank)
	
	df = df[correct.order,]
	rownames(df) = NULL
	pwms = pwms[correct.order]
	
	# return a MotifEnrichmentReport object
	new("MotifEnrichmentReport", d=df, pwms=pwms)
})

#' Generate a motif enrichment report for a single sequence
#'
#' @name sequenceReport,MotifEnrichmentResults-method
#' @aliases sequenceReport
#' @param obj a MotifEnrichmentResults object
#' @param seq.id the sequence index or name
#' @param bg if to use background corrected P-values to do the ranking (if available)
#' @param ... unused
#' @return a MotifEnrichmentReport object containing a table with the following columns: 
#' \itemize{
#'    \item 'rank' - The rank of the PWM's enrichment in the sequence
#'    \item 'target' - The name of the PWM's target gene, transcript or protein complex. 
#'    \item 'id' - The unique identifier of the PWM (if set during PWM creation).
#'    \item 'raw.score' - The raw score before P-value calculation
#'    \item 'p.value' - The P-value of motif enrichment (if available)
#' } 
#' @export
#' @examples
#' if(require("PWMEnrich.Dmelanogaster.background")){
#'    ###
#'    # load the pre-compiled lognormal background
#'    data(PWMLogn.dm3.MotifDb.Dmel)
#'
#'    # scan two sequences for motif enrichment
#'    sequences = list(DNAString("GAAGTATCAAGTGACCAGTAAGTCCCAGATGA"), DNAString("AGGTAGATAGAACAGTAGGCAATGAAGCCGATG"))
#'
#'    res = motifEnrichment(sequences, PWMLogn.dm3.MotifDb.Dmel)
#'
#'    # reports for the two sequences
#'    r1 = sequenceReport(res, 1)
#'    r2 = sequenceReport(res, 2)
#'
#'    # view the results
#'    r1
#'    r2
#'
#'    # plot the top 10 most enriched motifs in the first, and then second sequence
#'    plot(r1[1:10])
#'    plot(r2[1:10])
#' 
#' }
setMethod("sequenceReport", signature=signature(obj="MotifEnrichmentResults"), function(obj, seq.id, bg=TRUE, ...){
	if(missing(seq.id)){
		stop("Please specify the sequence number with 'seq.id'")
	}

	pwms = obj$pwms	
		
	# correct ordering of motifs
	o = motifRankingForSequence(obj, seq.id, rank=TRUE, bg=bg)
	
	targets = sapply(obj$pwms, function(x) x@name)
	ids = sapply(obj$pwms, function(x) x@id)
	
	# construct the initial data frame object
	df = data.frame(rank=o, target=targets, id=ids, raw.score=obj$sequence.nobg[seq.id,], stringsAsFactors=FALSE)
	
	if(bg && "sequence.bg" %in% names(obj)){
		p.value = obj$sequence.bg[seq.id,]
	} else {
		p.value = as.numeric(NA)
	}

	df = data.frame(df, p.value=p.value, stringsAsFactors=FALSE)
		
	# sort and return
	correct.order = order(df$rank)
	
	df = df[correct.order,]
	rownames(df) = NULL
	pwms = pwms[correct.order]
	
	# return a MotifEnrichmentReport object
	new("MotifEnrichmentReport", d=df, pwms=pwms)
})




