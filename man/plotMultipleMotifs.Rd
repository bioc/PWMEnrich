\name{plotMultipleMotifs}
\alias{plotMultipleMotifs}
\title{Plot mulitple motifs in a single plot...}
\usage{plotMultipleMotifs(pwms, titles=names(pwms), rows=ceiling(sqrt(length(pwms))),
    cols=ceiling(sqrt(length(pwms))), xmargin.scale=1, ymargin.scale=1,
    ...)
}
\description{Plot mulitple motifs in a single plot}
\details{Individual motif logos are plotted on a rows x cols grid. This function is a convenience
interface for the \code{seqLogoGrid} function that deals with viewpoint placement in a 
matrix-like grid layout. 

By default will try to make a square grid plot that would fit all the motifs and use
list names as captions.}
\arguments{\item{pwms}{a list of PWM objects or frequency matrices}
\item{titles}{a characater vector of titles for each of the plots}
\item{rows}{number of rows in the grid}
\item{cols}{number or cols in the grid}
\item{xmargin.scale}{the scaling parameter for the X-axis margin. Useful when plotting more than one logo on a page}
\item{ymargin.scale}{the scaling parameter for the Y-axis margin. Useful when plotting more than one logo on a page}
\item{...}{other parameters passed to seqLogoGrid()}
}

