\name{motifSimilarity}
\alias{motifSimilarity}
\title{Calculates similarity between two PFMs.}
\usage{
  motifSimilarity(m1, m2, trim = 0.4, self.sim = FALSE)
}
\arguments{
  \item{m1}{matrix with four rows representing the
  frequency matrix of first motif}

  \item{m2}{matrix with four rows representing the
  frequency matrix of second motif}

  \item{trim}{bases with information content smaller than
  this value will be trimmed off both motif ends}

  \item{self.sim}{if to calculate self similarity (i.e.
  without including offset=0 in alignment)}
}
\description{
  This function calculates the normalized motif correlation
  as a measure of motif frequency matrix similarity.
}
\details{
  This score is essentially a normalized version of the sum
  of column correlations as proposed by Pietrokovski
  (1996). The sum is normalized by the average motif length
  of m1 and m2, i.e. (ncol(m1)+ncol(m2))/2. Thus, for two
  idential motifs this score is going to be 1. For
  unrelated motifs the score is going to be typically
  around 0.

  Motifs need to aligned for this score to be calculated.
  The current implementation tries all possible ungapped
  alignment with a minimal of two basepair matching, and
  the maximal score over all alignments is returned.

  Motif 1 is aligned both to Motif 2 and its reverse
  complement. Thus, the motif similarities are the same if
  the reverse complement of any of the two motifs is given.
}
\examples{
if(require("PWMEnrich.Dmelanogaster.background")){
   data(MotifDb.Dmel.PFM)

    # calculate the similarity of tin and vnd motifs (which are almost identical)
    motifSimilarity(MotifDb.Dmel.PFM[["Dmelanogaster-JASPAR_CORE-tin-MA0247.1"]], MotifDb.Dmel.PFM[["Dmelanogaster-JASPAR_CORE-vnd-MA0253.1"]])

    # similarity of two unrelated motifs
    motifSimilarity(MotifDb.Dmel.PFM[["Dmelanogaster-JASPAR_CORE-tin-MA0247.1"]], MotifDb.Dmel.PFM[["Dmelanogaster-JASPAR_CORE-ttk-MA0460.1"]])
}
}
\references{
  Pietrokovski S. Searching databases of conserved sequence
  regions by aligning protein multiple-alignments. Nucleic
  Acids Res 1996;24:3836-3845.
}

