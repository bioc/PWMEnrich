\name{reverseComplement,PWM-method}
\alias{reverseComplement,PWM-method}
\alias{reverseComplement}
\title{Reverse complement for the PWM object...}
\usage{\S4method{reverseComplement}{PWM}(x, ...)
}
\description{Reverse complement for the PWM object}
\details{Finds the reverse complement of the PWM}
\value{an object of type PWM that is reverse complement of x}
\arguments{\item{x}{an object of type PWM}
\item{...}{unused}
}
\examples{if(require("PWMEnrich.Dmelanogaster.background")){
data(jaspar.insects.PFM)

reverseComplement(jaspar.insects.PFM$ttk) # reverse complement of the ttk PWM
}}

