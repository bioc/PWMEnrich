\name{readTRANSFAC}
\alias{readTRANSFAC}
\title{Read in motifs in TRANSFAC format}
\usage{
  readTRANSFAC(file, remove.acc = TRUE)
}
\arguments{
  \item{file}{the filename}

  \item{remove.acc}{if to ignore transfac accession
  numbers}
}
\value{
  a list of matrices representing motifs (with four
  nucleotides as rows)
}
\description{
  Read in motifs in TRANSFAC format
}

