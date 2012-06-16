\name{PWMGEVBackground-class}
\alias{PWMGEVBackground-class}
\alias{PWMGEVBackground}
\title{Generalized Extreme Values (GEV) background for P-values...}
\description{Generalized Extreme Values (GEV) background for P-values}
\details{The three parameters of the GEV distribution are fitted by doing linear regression
on log of sequence length.}
\section{Slots}{\describe{\item{\code{bg.source}:}{(\code{\link{character}}) textual description of where the background distribution is derived from}\item{\code{bg.loc}:}{(\code{\link{list}}) linear regression model for estimating the location parameter based on log(L), list of lm objects of PWMs}\item{\code{bg.scale}:}{(\code{\link{list}}) linear regression model for estimating the scale parameter based on log(L), list of lm objects of PWMs}\item{\code{bg.shape}:}{(\code{\link{list}}) linear regression model for estimating the shape parameter based on log(L), list of lm objects of PWMs}\item{\code{pwms}:}{(\code{\link{list}}) the pwms for which the background has been compiled}}}
\section{Methods}{\describe{\item{\code{\link[=show,PWMGEVBackground-method]{show}}}{\code{signature(object = "PWMGEVBackground")}: show method for PWMGEVBackground}\item{\code{\link[=names,PWMGEVBackground-method]{names}}}{\code{signature(x = "PWMGEVBackground")}: Name of different pieces of information associated with PWMGEVBackground}\item{\code{\link[=$,PWMGEVBackground-method]{$}}}{\code{signature(x = "PWMGEVBackground")}: Access a property by name}}}

