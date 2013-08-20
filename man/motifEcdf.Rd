\name{motifEcdf}
\alias{motifEcdf}
\title{Calculate the empirical distribution score distribution for a set of motifs}
\usage{
  motifEcdf(motifs, organism = NULL, bg.seq = NULL,
    quick = FALSE, pseudo.count = 1)
}
\arguments{
  \item{motifs}{a set of motifs, either a list of frequency
  matrices, or a list of PWM objects. If frequency matrices
  are given, the background distribution is fitted from
  bg.seq.}

  \item{organism}{either a name of the organisms for which
  the background should be compiled (currently only
  supported name is "dm3" for Drosophila Melanogaster), or
  a \code{BSgenome} object (see \code{BSgenome} package).}

  \item{bg.seq}{a set of background sequence (either this
  or organism needs to be specified!). Can be a DNAString
  or DNAStringSet object.}

  \item{quick}{if to do the fitting only on a small subset
  of the data (only in combination with \code{organism}).
  Useful only for code testing!}

  \item{pseudo.count}{the pseudo count which is shared
  between nucleotides when frequency matrices are given}
}
\value{
  a list of \code{ecdf} objects (see help page for
  \code{ecdf} for usage).
}
\description{
  Calculate the empirical distribution score distribution
  for a set of motifs
}

