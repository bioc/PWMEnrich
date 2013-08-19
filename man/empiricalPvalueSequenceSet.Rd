\name{empiricalPvalueSequenceSet}
\alias{empiricalPvalueSequenceSet}
\title{Empirical P-value for a set of sequences}
\usage{
  empiricalPvalueSequenceSet(scores, seq.len, pwm.len,
    bg.fwd, bg.rev, cutoff = NULL, B = 10000,
    verbose = FALSE)
}
\arguments{
  \item{scores}{a matrix of scores, rows for sequences,
  columns for PWMs}

  \item{seq.len}{the lengths of sequences}

  \item{pwm.len}{the lengths of PWMs}

  \item{bg.fwd}{raw odds scores for the forward strand of
  background}

  \item{bg.rev}{raw odds scores for the reverse strand of
  background}

  \item{cutoff}{if not NULL, will use hit count above this
  cutoff. The cutoff should be specified in log2.}

  \item{B}{the number of random replicates}

  \item{verbose}{if to give verbose progress reports}
}
\description{
  Calculate empirical P-value for a set of sequences, using
  either affinity or cutoff. When cutoff is used, the score
  is a number of motif hits above a certain log-odds
  cutoff.
}

