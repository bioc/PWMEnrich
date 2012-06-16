\name{seqLogoGrid}
\alias{seqLogoGrid}
\title{Draw a motif logo on an existing viewport...}
\usage{seqLogoGrid(pwm, ic.scale=TRUE, xaxis=TRUE, yaxis=TRUE, xfontsize=15, yfontsize=15,
    xmargin.scale=1, ymargin.scale=1, title="", titlefontsize=15)
}
\description{Draw a motif logo on an existing viewport}
\details{This function comes from the seqLogo package. It has been modified to remove 
some unneccessary code as suggested by W Huber (https://stat.ethz.ch/pipermail/bioconductor/2010-September/035267.html).

Use this function for more advanced plotting where the viewports are directly set up and maintained (see package \code{grid}).}
\arguments{\item{pwm}{numeric The 4xW position weight matrix.}
\item{ic.scale}{logical If TRUE, the height of each column is proportional to its information
content. Otherwise, all columns have the same height.}
\item{xaxis}{logical If TRUE, an X-axis will be plotted.}
\item{yaxis}{logical If TRUE, a Y-axis will be plotted.}
\item{xfontsize}{numeric Font size to be used for the X-axis.}
\item{yfontsize}{numeric Font size to be used for the Y-axis.}
\item{xmargin.scale}{the scaling parameter for the X-axis margin. Useful when plotting more than one log on a page}
\item{ymargin.scale}{the scaling parameter for the Y-axis margin. Useful when plotting more than one log on a page}
\item{title}{to be shown on the top}
\item{titlefontsize}{the fontsize of the title}
}

