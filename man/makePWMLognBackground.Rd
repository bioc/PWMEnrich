% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/background.R
\name{makePWMLognBackground}
\alias{makePWMLognBackground}
\title{Make a lognormal background distribution}
\usage{
makePWMLognBackground(
  bg.seq,
  motifs,
  bg.pseudo.count = 1,
  bg.len = 250,
  bg.len.sizes = 2^(0:4),
  bg.source = "",
  verbose = TRUE,
  algorithm = "default"
)
}
\arguments{
\item{bg.seq}{a set of background sequences, either a list of DNAString object or DNAStringSet object}

\item{motifs}{a set of motifs, either a list of frequency matrices, or a list of PWM objects. If
frequency matrices are given, the background distribution is fitted from bg.seq.}

\item{bg.pseudo.count}{the pseudo count which is shared between nucleotides when frequency matrices are given}

\item{bg.len}{background sequences will be split into tiles of this length (default: 250bp)}

\item{bg.len.sizes}{background tiles will be joined into bigger tiles containing this much smaller tiles. 
The default is \code{2^(0:4)}, which with \code{bg.len} translates into 
250bp, 500bp, 1000bp, 1500bp, 2000bp, 4000bp. Note this is only used in the "human" algorithm.}

\item{bg.source}{a free-form textual description of how the background was generated}

\item{verbose}{if to produce verbose output}

\item{algorithm}{type of algorithm to use, valid values are: "default" and "human".}
}
\description{
Construct a lognormal background distribution for a set of sequences. 
Sequences concatenated are binned in 'bg.len' chunks and lognormal distribution 
fitted to them.
}
\examples{
\dontrun{
if(requireNamespace("PWMEnrich.Dmelanogaster.background")){
   data(MotifDb.Dmel.PFM, package = "PWMEnrich.Dmelanogaster.background")

   # make background for MotifDb motifs using 2kb promoters of all D. melanogaster transcripts 
	  if(requireNamespace("BSgenome.Dmelanogaster.UCSC.dm3")) 
     makePWMLognBackground(Dmelanogaster$upstream2000, MotifDb.Dmel.PFM)
}
}
}
