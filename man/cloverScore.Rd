\name{cloverScore}
\alias{cloverScore}
\title{Calculate the Clover score using the recursive formula from Frith et al...}
\usage{cloverScore(scores, lr3=FALSE, verbose=FALSE)
}
\description{Calculate the Clover score using the recursive formula from Frith et al}
\value{the LR4 score, which is the mean of LR3 scores over subset sizes}
\arguments{\item{scores}{a matrix of average odds scores, where columns are motifs, and rows sequences}
\item{lr3}{if to return a matrix of LR3 scores, where columns correpond to motifs, and rows to subset sizes}
\item{verbose}{if to produce verbose output of progress}
}

