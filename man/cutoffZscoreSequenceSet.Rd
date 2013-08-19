\name{cutoffZscoreSequenceSet}
\alias{cutoffZscoreSequenceSet}
\title{Z-score calculation for cutoff hits for group of sequences}
\usage{
  cutoffZscoreSequenceSet(scores, seq.len, pwm.len, bg.P)
}
\arguments{
  \item{scores}{the hit counts for the sequences}

  \item{seq.len}{the length distribution of sequences}

  \item{pwm.len}{the length distribution of the PWMs}

  \item{bg.P}{background probabilities of observing a motif
  hit at nucleotide resolution}
}
\value{
  Z-score
}
\description{
  The Z-score is calculated as if the sequence came for one
  very long sequence
}

