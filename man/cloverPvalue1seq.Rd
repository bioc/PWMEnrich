\name{cloverPvalue1seq}
\alias{cloverPvalue1seq}
\title{Calculate the Clover P-value as described in the Clover paper}
\usage{
  cloverPvalue1seq(scores, seq.len, pwm.len, bg.fwd,
    bg.rev, B = 1000, verbose = TRUE, clover = NULL)
}
\arguments{
  \item{scores}{the affinity scores for individual
  sequences}

  \item{seq.len}{lengths of sequences}

  \item{pwm.len}{lengths of PWMs}

  \item{bg.fwd}{the raw score of forward strand}

  \item{bg.rev}{the raw scores of reverse strand}

  \item{B}{the number of random replicates}

  \item{verbose}{if to give verbose progress reports}

  \item{clover}{the clover scores if already calculated}
}
\value{
  P-value
}
\description{
  This function only take one background sequence as input,
  it also just calculates the P-value so it is more
  efficient.
}

