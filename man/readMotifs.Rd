\name{readMotifs}
\alias{readMotifs}
\title{Read in motifs in JASPAR or TRANSFAC format...}
\usage{readMotifs(file, remove.acc=FALSE)
}
\description{Read in motifs in JASPAR or TRANSFAC format}
\details{The format is autodetected based on file format. If the autodetection fail
then the file cannot be read.}
\value{a list of 4xL matrices representing motifs (four nucleotides as rows)}
\arguments{\item{file}{the filename}
\item{remove.acc}{if to remove accession numbers. If TRUE, the AC entry in TRANSFAC files is ignored,
and the accession is stripped from JASPAR, e.g. motif with name "MA0211.1 bap" would 
become just "bap". If FALSE, botht he AC and ID are used to generate the TRANSFAC name
and the original motif names are preserved in JASPAR files.}
}
\examples{# read in example TRANSFAC motifs without accession codes (just IDs)
readMotifs(system.file(package="PWMEnrich", dir="extdata", file="example.transfac"), remove.acc=TRUE)

# read in the JASPAR insects motifs provided as example
readMotifs(system.file(package="PWMEnrich", dir="extdata", file="jaspar-insecta.jaspar"), remove.acc=TRUE)}

