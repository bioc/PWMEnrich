\name{maxAligned}
\alias{maxAligned}
\title{Returned the aligned motif parts}
\usage{
  maxAligned(m1, m2, offset)
}
\arguments{
  \item{m1}{frequency matrix of first motif}

  \item{m2}{frequency matrix of second motif}

  \item{offset}{a number of nucleotides by which the first
  motif is offsetted compared to the second}
}
\value{
  a list of column-trimmed motifs m1, m2
}
\description{
  This function takes the offset of first motif relative to
  second and chops off the end of both motifs that are not
  aligned. It returns a list containing only the columns
  that align.
}

