### Functions for differential enrichment

#' Test for differential enrichment between two groups of sequences
#'
#' This function calls \code{motifEnrichment} on two groups of sequences and calculates
#' the difference statistics when possible. 
#'
#' @title Differential motif enrichment
#' @param sequences1 First set of sequences. Can be either a single sequence 
#'        (an object of class DNAString), or a list of DNAString objects, or a DNAStringSet object.
#' @param sequences2 Second set of sequences. Can be either a single sequence 
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
#'           \item "affinity" - use threshold-free affinity scores without a background. The \code{pwms}
#'                 parameter can either be a list of frequency matrices, \code{PWM} objects, or a 
#'                 \code{PWMLognBackground} object. 
#'           \item "cutoff" - use number of motif hits above a score cutoff as a measure of enrichment. 
#'                 No background correction is performed. The \code{pwms}
#'                 parameter can either be a list of frequency matrices, \code{PWM} objects, or a 
#'                 \code{PWMCutoffBackground} object.
#'         }
#'            
#' @param bg this parameter determines which background correction to use, if any. 
#'        \itemize{
#'           \item "autodetect" - default value. Background correction is determined based on the type
#'                  of the \code{pwms} parameter. 
#'           \item "logn" - use a lognormal distribution background pre-computed for a set of PWMs.
#'                 This requires \code{pwms} to be of class \code{PWMLognBackground}.
#'           \item "z" - use a z-score for the number of significant motif hits compared to background number of hits.
#'                 This requires \code{pwms} to be of class \code{PWMCutoffBackground}.
#'           \item "none" - no background correction
#'         }
#' @param cutoff the score cutoff for a significant motif hit if scoring scheme "cutoff" is selected. 
#' @param res1 the output of \code{motifEnrichment} if already calculated for \code{sequences1}
#' @param res2 the output of \code{motifEnrichment} if already calculated for \code{sequences2}
#' @param verbose if to produce verbose output
#' @export
#' @examples
#'
#' if(require("PWMEnrich.Dmelanogaster.background")){
#'   # load the background file for drosophila and lognormal correction
#'   data(PWMLogn.dm3.MotifDb.Dmel)
#'
#'   # get the differential enrichment
#'   diff = motifDiffEnrichment(DNAString("TGCATCAAGTGTGTAGTGTGAGATTAGT"), DNAString("TGAACGAGTAGGACGATGAGAGATTGATG"), PWMLogn.dm3.MotifDb.Dmel, verbose=FALSE)
#'
#'   # motifs differentially enriched in the first sequence (with lognormal background correction)
#'   head(sort(diff$group.bg, decreasing=TRUE))
#'
#'   # motifs differentially enriched in the second sequence (with lognormal background correction)
#'   head(sort(diff$group.bg))
#' }
motifDiffEnrichment = function(sequences1, sequences2, pwms, score="autodetect", bg="autodetect", cutoff=log2(exp(4)), verbose=TRUE,
	res1=NULL, res2=NULL){
	# check the input param constrains
	if(!inherits(pwms, c("NULL", "list", "PWMLognBackground", "PWMCutoffBackground"))) {
		stop("pwms needs to be either a list of frequency matrices, list of PWM objects, or an object of class PWMLognBackground or PWMCutoffBackground.")
	}
	
	if(!(score %in% c("autodetect", "affinity", "cutoff"))){
		stop("score needs to be one of: autodetect, affinity, cutoff")
	}
	
	if(!(bg %in% c("autodetect", "logn", "z", "none"))){
		stop("score needs to be one of: autodetect, logn, z, none")
	}

	# results for both groups
	if(is.null(res1))
		res1 = motifEnrichment(sequences1, pwms, score=score, bg=bg, cutoff=cutoff, verbose=verbose)
	if(is.null(res2))
		res2 = motifEnrichment(sequences2, pwms, score=score, bg=bg, cutoff=cutoff, verbose=verbose)
		
	if((res1$score != res2$score) | (res1$bg != res2$bg))
		stop("Supplied motif enrichments for sequences1 and sequences2 use different scoring schemes and/or background corrections")
	
	# we use the same rules to decided on score and bg.. 
	if(score == "autodetect")
		score = res1$score
	if(bg == "autodetect")
		bg = res1$bg

	# calculated differential enrichment
	res = list()
	res$group.nobg = res1$group.nobg - res2$group.nobg	
	
	if(bg == "none")
		res$group.bg = NULL
	else if(bg == "logn")
		res$group.bg = res1$group.norm - res2$group.norm
	else
		res$group.bg = res1$group.bg - res2$group.bg
	
	return(res)
}


