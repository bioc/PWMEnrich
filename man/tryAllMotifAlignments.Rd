\name{tryAllMotifAlignments}
\alias{tryAllMotifAlignments}
\title{Try all motif alignments and return max score}
\usage{
  tryAllMotifAlignments(m1, m2, min.align = 2,
    exclude.zero = FALSE)
}
\arguments{
  \item{m1}{frequency matrix of motif 1}

  \item{m2}{frequency matrix of motif 2}

  \item{min.align}{minimal number of basepairs that need to
  align}

  \item{exclude.zero}{if to exclude offset=0, useful for
  calculating self-similarity}
}
\value{
  single maximal score
}
\description{
  This function tries all offsets of motif1 compared to
  motif2 and returns the maximal (unnormalized) correlation
  score.
}
\details{
  The correlation score is essentially the sum of
  correlations of individual aligned columns as described
  in Pietrokovski (1996).
}
\references{
  Pietrokovski S. Searching databases of conserved sequence
  regions by aligning protein multiple-alignments. Nucleic
  Acids Res 1996;24:3836-3845.
}

