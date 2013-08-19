
#' A class that represents a Position Weight Matrix (PWM)
#'
#' @aliases PWM-class
#' @section Slots: 
#'  \describe{
#'    \item{\code{id}:}{a systematic ID given to this PWM, could include the source, version, etc}
#'    \item{\code{name}:}{the name of the transcription factor (TF) to which the PWM corresponds to}
#'    \item{\code{pfm}:}{Position Frequency Matrix (PFM) from which the PWM is derived}
#'    \item{\code{prior.params}:}{Defines prior frequencies of the four bases (A,C,G,T), a named vector. These will be added to individual values for the PFM and at the same time used as background probabilities}
#'    \item{\code{pwm}:}{Final Position Weight Matrix (PWM) constructed using prior.params with logarithm base 2}
#'  }
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
#' @aliases PWMLognBackground-class
#' @section Slots: 
#'  \describe{
#'    \item{\code{bg.source}:}{textual description of where the background distribution is derived from}
#'    \item{\code{bg.len}:}{the length to which the background is normalized to. This is a vector of values, can have a different value for each motif.}
#'    \item{\code{bg.mean}:}{the mean value of the lognormal distribution at bg.len}
#'    \item{\code{bg.sd}:}{the standard deviation of the lognormal distribution at bg.len}
#'    \item{\code{pwms}:}{the pwms for which the background has been compiled}
#'  }
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
#' @aliases PWMCutoffBackground-class
#' @section Slots: 
#'  \describe{
#'    \item{\code{bg.source}:}{textual description of where the background distribution is derived from}
#'    \item{\code{bg.cutoff}:}{the cutoff score used to find significant motif hits (in log2 odds), either a single value or a vector of values}
#'    \item{\code{bg.P}:}{the density of significant motif hits per nucleotide in background}
#'    \item{\code{pwms}:}{the pwms for which the background has been compiled}
#'  }
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
#' @aliases PWMEmpiricalBackground-class
#' @section Slots: 
#'  \describe{
#'    \item{\code{bg.source}:}{textual description of where the background distribution is derived from}
#'    \item{\code{bg.fwd}:}{affinity scores (odds) for the forward strand. PWMs as columns}
#'    \item{\code{bg.rev}:}{affinity scores (odds) for the reverse strand. PWMs as columns}
#'    \item{\code{pwms}:}{the pwms for which the background has been compiled}
#'  }
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
#' @aliases PWMGEVBackground-class
#' @section Slots: 
#'  \describe{
#'    \item{\code{bg.source}:}{textual description of where the background distribution is derived from}
#'    \item{\code{bg.loc}:}{linear regression model for estimating the location parameter based on log(L), list of lm objects of PWMs}
#'    \item{\code{bg.scale}:}{linear regression model for estimating the scale parameter based on log(L), list of lm objects of PWMs}
#'    \item{\code{bg.shape}:}{linear regression model for estimating the shape parameter based on log(L), list of lm objects of PWMs}
#'    \item{\code{pwms}:}{the pwms for which the background has been compiled}
#'  }
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
#' @aliases MotifEnrichmentResults-class
#' @section Slots: 
#'  \describe{
#'    \item{\code{res}:}{a list of old results with elements such as: sequence.bg, sequence.nobg, group.bg, group.nobg}
#'  }
#' @export
setClass("MotifEnrichmentResults",
	representation = representation(
		res = "list"
	)
)

#' A report class with formatted results of motif enrichment
#'
#' The columns stored in this object will depend on the type of the report (either for group of sequences, or individual sequences).
#'
#' @aliases MotifEnrichmentReport-class
#' @section Slots: 
#'  \describe{
#'    \item{\code{d}:}{a DataFrame object that contains the main tabular report data}
#'    \item{\code{pwms}:}{a list of \code{PWM} objects corresponding to rows of \code{d}}
#'  }
#' @export
setClass("MotifEnrichmentReport",
	representation = representation(
		d = "data.frame",
		pwms = "list"
	)
)



