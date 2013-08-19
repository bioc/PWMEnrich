\name{motifScoresBigMemory}
\alias{motifScoresBigMemory}
\title{This is a memory intensive version of motifScore() which is about 2 times faster}
\usage{
  motifScoresBigMemory(sequences, motifs,
    raw.scores = FALSE, verbose = TRUE, cutoff = NULL,
    seq.all = NULL)
}
\arguments{
  \item{sequences}{set of input sequences}

  \item{motifs}{set of input PWMs or PFMs}

  \item{raw.scores}{if to return scores for each base-pair}

  \item{verbose}{if to produce verbose output}

  \item{cutoff}{the cutoff for calling binding sites (in
  base 2 log).}

  \item{seq.all}{already concatenated sequences if already
  available (used to internally speed up things)}
}
\description{
  The parameters and functionality are the same as
  \code{\link{motifScores}}. Please refer to documentation
  of this function for detailed explanation of
  functionality.
}
\details{
  This function is not meant to be called directly, but is
  indirectly called by motifScores() once a global
  parameters useBigMemory is set.
}
\seealso{
  \code{\link{motifScores}}
}

