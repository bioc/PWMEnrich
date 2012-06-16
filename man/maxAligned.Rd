\name{maxAligned}
\alias{maxAligned}
\title{Returned the aligned motif parts...}
\usage{maxAligned(m1, m2, offset)
}
\description{Returned the aligned motif parts}
\details{This function takes the offset of first motif relative to second and 
chops off the end of both motifs that are not aligned. It returns a list
containing only the columns that align.}
\value{a list of column-trimmed motifs m1, m2}
\arguments{\item{m1}{frequency matrix of first motif}
\item{m2}{frequency matrix of second motif}
\item{offset}{a number of nucleotides by which the first motif is offsetted compared to the second}
}

