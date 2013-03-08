\name{plot,PWM,missing-method}
\alias{plot,PWM,missing-method}
\alias{plot}
\title{Plotting for the PWM class...}
\usage{\S4method{plot}{PWM,missing}(x, y, ...)
}
\description{Plotting for the PWM class}
\details{This function produces a sequence logo (via package seqLogo).}
\arguments{\item{x}{the PWM object}
\item{y}{unused}
\item{...}{other parameters to pass to seqLogo's \code{plot} function}
}
\examples{if(require("PWMEnrich.Dmelanogaster.background")){
data(MotifDb.Dmel)

# plot the tinman motif from MotifDb
plot(MotifDb.Dmel$tin)
}}

