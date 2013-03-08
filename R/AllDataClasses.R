
#' A class that represents a Position Weight Matrix (PWM)
#'
#' @slot id a systematic ID given to this PWM, could include the source, version, etc
#' @slot name the name of the transcription factor (TF) to which the PWM corresponds to
#' @slot pfm Position Frequency Matrix (PFM) from which the PWM is derived
#' @slot prior.params Defines prior frequencies of the four bases (A,C,G,T), a named vector. These will be added 
#'                    to individual values for the PFM and at the same time used as background probabilities
#' @slot pwm Final Position Weight Matrix (PWM) constructed using prior.params with logarithm base 2
#' @export
setClass("PWM", 
	representation = representation(
		id = "character",
		name = "character",
		pfm = "matrix",
		prior.params = "vector",
		pwm = "matrix"
	)
)

#' Lognormal background distribution for a set of PWMs
#'
#' @slot bg.source textual description of where the background distribution is derived from
#' @slot bg.len the length to which the background is normalized to. This is a vector of values, 
#'              can have a different value for each motif. 
#' @slot bg.mean the mean value of the lognormal distribution at bg.len
#' @slot bg.sd the standard deviation of the lognormal distribution at bg.len
#' @slot pwms the pwms for which the background has been compiled
#' @export
setClass("PWMLognBackground",
	representation = representation(
		bg.source = "character",
		bg.len = "numeric",
		bg.mean = "numeric",
		bg.sd = "numeric",
		pwms = "list"				
	)	
)

#' Hit count background distribution for a set of PWMs
#'
#' @slot bg.source textual description of where the background distribution is derived from
#' @slot bg.cutoff the cutoff score used to find significant motif hits (in log2 odds),
#'                 either a single value or a vector of values
#' @slot bg.P the density of significant motif hits per nucleotide in background
#' @slot pwms the pwms for which the background has been compiled
#' @export
setClass("PWMCutoffBackground",
	representation = representation(
		bg.source = "character",
		bg.cutoff = "numeric",
		bg.P = "numeric",
		pwms = "list"
	)
) 

#' Background for calculating empirical P-values
#'
#' This object contains raw scores for one very long sequence, thus it can be very large. 
#'
#' @slot bg.source textual description of where the background distribution is derived from
#' @slot bg.fwd affinity scores (odds) for the forward strand. PWMs as columns. 
#' @slot bg.rev affinity scores (odds) for the reverse strand. PWMs as columns.
#' @slot pwms the pwms for which the background has been compiled
#' @export
setClass("PWMEmpiricalBackground",
	representation = representation(
		bg.source = "character",
		bg.fwd = "matrix",
		bg.rev = "matrix",
		pwms = "list"
	)
) 

#' Generalized Extreme Values (GEV) background for P-values
#'
#' The three parameters of the GEV distribution are fitted by doing linear regression
#' on log of sequence length. 
#'
#' @slot bg.source textual description of where the background distribution is derived from
#' @slot bg.loc linear regression model for estimating the location parameter based on log(L), list of lm objects of PWMs
#' @slot bg.scale linear regression model for estimating the scale parameter based on log(L), list of lm objects of PWMs
#' @slot bg.shape linear regression model for estimating the shape parameter based on log(L), list of lm objects of PWMs
#' @slot pwms the pwms for which the background has been compiled
#' @export
setClass("PWMGEVBackground",
	representation = representation(
		bg.source = "character",
		bg.loc = "list",
		bg.scale = "list",
		bg.shape = "list",
		pwms = "list"				
	)	
)

#' A wrapper class for results of motifEnrichment() that should make it easier to access the results. 
#'
#' Note that this is only a wrapper around a list which is the return value in PWMEnrich 1.3 and as such it
#' provides the same interface as a list (for backward compatibility), with some additional methods. 
#'
#' @slot res a list of old results with elements such as: sequence.bg, sequence.nobg, group.bg, group.nobg
#' @export
setClass("MotifEnrichmentResults",
	representation = representation(
		res = "list"
	)
)


