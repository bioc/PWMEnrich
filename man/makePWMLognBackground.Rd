\name{makePWMLognBackground}
\alias{makePWMLognBackground}
\title{Make a lognormal background distribution}
\usage{
  makePWMLognBackground(bg.seq, motifs,
    bg.pseudo.count = 1, bg.len = 1000, bg.source = "",
    verbose = TRUE)
}
\arguments{
  \item{bg.seq}{a set of background sequences, either a
  list of DNAString object or DNAStringSet object}

  \item{motifs}{a set of motifs, either a list of frequency
  matrices, or a list of PWM objects. If frequency matrices
  are given, the background distribution is fitted from
  bg.seq.}

  \item{bg.pseudo.count}{the pseudo count which is shared
  between nucleotides when frequency matrices are given}

  \item{bg.len}{the length of background chunks}

  \item{bg.source}{a free-form textual description of how
  the background was generated}

  \item{verbose}{if to produce verbose output}
}
\description{
  Construct a lognormal background distribution for a set
  of sequences. Sequences concatenated are binned in
  'bg.len' chunks and lognormal distribution fitted to
  them.
}
\examples{
\dontrun{
if(require("PWMEnrich.Dmelanogaster.background")){
   data(MotifDb.Dmel.PFM)

   # make background for MotifDb motifs using 2kb promoters of all D. melanogaster transcripts
	  if(require("BSgenome.Dmelanogaster.UCSC.dm3"))
     makePWMLognBackground(Dmelanogaster$upstream2000, MotifDb.Dmel.PFM)
}
}
}

