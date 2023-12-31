% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/background.R
\name{makePWMPvalCutoffBackground}
\alias{makePWMPvalCutoffBackground}
\title{Construct a cutoff background from empirical background}
\usage{
makePWMPvalCutoffBackground(bg.p, p.value = 0.001, bg.source = "")
}
\arguments{
\item{bg.p}{an object of class PWMEmpiricalBackground}

\item{p.value}{the P-value used to find cuttoffs for each of the motifs}

\item{bg.source}{textual description of background source}
}
\value{
an object of type PWMCutoffBackground
}
\description{
This function takes already calculated empirical background distribution and chooses
cutoff for each motif based on P-value cutoff for individual sites.
}
\examples{
\dontrun{
if(requireNamespace("PWMEnrich.Dmelanogaster.background")){
   data(MotifDb.Dmel.PFM, package = "PWMEnrich.Dmelanogaster.background")

   # make empirical background - here we use only 100 sequences for illustrative purposes
	  if(requireNamespace("BSgenome.Dmelanogaster.UCSC.dm3")) 
     bg.p = makePWMEmpiricalBackground(Dmelanogaster$upstream2000[1:100], MotifDb.Dmel.PFM)

   # use the empirical background to pick a threshold and make cutoff background
   makePWMPvalCutoffBackground(bg.p, 0.001)
}
}
}
