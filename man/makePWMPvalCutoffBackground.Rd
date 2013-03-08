\name{makePWMPvalCutoffBackground}
\alias{makePWMPvalCutoffBackground}
\title{Construct a cutoff background from empirical background...}
\usage{makePWMPvalCutoffBackground(bg.p, p.value=0.001, bg.source="")
}
\description{Construct a cutoff background from empirical background}
\details{This function takes already calculated empirical background distribution and chooses
cutoff for each motif based on P-value cutoff for individual sites.}
\value{an object of type PWMCutoffBackground}
\arguments{\item{bg.p}{an object of class PWMEmpiricalBackground}
\item{p.value}{the P-value used to find cuttoffs for each of the motifs}
\item{bg.source}{textual description of background source}
}
\examples{\dontrun{
if(require("PWMEnrich.Dmelanogaster.background")){
data(MotifDb.Dmel.PFM)

# make empirical background - here we use only 100 sequences for illustrative purposes
if(require("BSgenome.Dmelanogaster.UCSC.dm3")) 
bg.p = makePWMEmpiricalBackground(Dmelanogaster$upstream2000[1:100], MotifDb.Dmel.PFM)

# use the empirical background to pick a threshold and make cutoff background
makePWMPvalCutoffBackground(bg.p, 0.001)
}
}}

