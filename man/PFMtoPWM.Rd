\name{PFMtoPWM}
\alias{PFMtoPWM}
\title{Convert frequencies into motifs using PWMUnscaled...}
\usage{PFMtoPWM(motifs, ...)
}
\description{Convert frequencies into motifs using PWMUnscaled}
\arguments{\item{motifs}{a list of motifs represented as matrices of frequencies (PFM)}
\item{...}{other parameters to PWMUnscaled}
}
\examples{if(require("PWMEnrich.Dmelanogaster.background")){
data(jaspar.insects.PFM)

PFMtoPWM(jaspar.insects.PFM) # convert to PWM with uniform background

prior = getBackgroundFrequencies("dm3", quick=TRUE) # get background for drosophila (quick mode on a reduced dataset)
PFMtoPWM(jaspar.insects.PFM, prior=prior) # convert with genomic background 
}}

