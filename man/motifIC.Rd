\name{motifIC}
\alias{motifIC}
\title{Information content for a PWM or PFM...}
\usage{motifIC(motif, prior.params=c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
    bycol=FALSE)
}
\description{Information content for a PWM or PFM}
\value{information content in bits (i.e. log2)}
\arguments{\item{motif}{a matrix of frequencies, or a PWM object}
\item{prior.params}{the prior parameters to use when a matrix is given (ignored if motif is already a PWM)}
\item{bycol}{if to return values separately for each column}
}
\examples{if(require("PWMEnrich.Dmelanogaster.background")){
data(jaspar.insects)
data(jaspar.insects.PFM)

motifIC(jaspar.insects$ttk) # the nucleotide distribution is taken from the PWM (in this case genomic background)
motifIC(jaspar.insects.PFM$ttk) # information content with default uniform background because the input is a matrix, not PWM object
}}

