\name{readJASPAR}
\alias{readJASPAR}
\title{Read motifs in JASPAR format}
\usage{
  readJASPAR(file, remove.ids = FALSE)
}
\arguments{
  \item{file}{the filename}

  \item{remove.ids}{if to strip JASPAR ID's from motif
  names, e.g. "MA0211.1 bap" would become just "bap"}
}
\value{
  a list of matrices representing motifs (with four
  nucleotides as rows)
}
\description{
  Read motifs in JASPAR format
}

