\name{PFMtoPWM}
\alias{PFMtoPWM}
\title{Convert frequencies into motifs using PWMUnscaled...}
\usage{PFMtoPWM(motifs, id=names(motifs), name=names(motifs), ...)
}
\description{Convert frequencies into motifs using PWMUnscaled}
\arguments{\item{motifs}{a list of motifs represented as matrices of frequencies (PFM)}
\item{id}{the set of IDs for the motifs (defaults to names of the 'motifs' list)}
\item{name}{the set of names for the motifs (defaults to names of the 'motifs' list)}
\item{...}{other parameters to PWMUnscaled}
}
\examples{if(require("PWMEnrich.Dmelanogaster.background")){
data(MotifDb.Dmel.PFM)

PFMtoPWM(MotifDb.Dmel.PFM) # convert to PWM with uniform background

prior = getBackgroundFrequencies("dm3", quick=TRUE) # get background for drosophila (quick mode on a reduced dataset)
PFMtoPWM(MotifDb.Dmel.PFM, prior.params=prior) # convert with genomic background 
}}

