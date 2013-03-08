\name{getBackgroundFrequencies}
\alias{getBackgroundFrequencies}
\title{Get the four nucleotides background frequencies...}
\usage{getBackgroundFrequencies(organism="dm3", pseudo.count=1, quick=FALSE)
}
\description{Get the four nucleotides background frequencies}
\details{Estimate the background frequencies of A,C,G,T on a set of promoters from an organism}
\author{Robert Stojnic, Diego Diez}
\arguments{\item{organism}{either a name of the organisms for which the background should be compiled 
(currently only supported name is "dm3" for Drosophila Melanogaster), or a \code{BSgenome} object (see \code{BSgenome} package).}
\item{pseudo.count}{the number to which the frequencies sum up to, by default 1}
\item{quick}{if to preform fitting on a reduced set of 100 promoters. This will not give as good results but is much quicker than fitting to all the promoters (~10k). 
Usage of this parameter is recommended only for testing and rough estimates.}
}
\examples{\dontrun{
getBackgroundFrequencies("dm3")
}}

