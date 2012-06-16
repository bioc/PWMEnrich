\name{PWMLognBackground-class}
\alias{PWMLognBackground-class}
\alias{PWMLognBackground}
\title{Lognormal background distribution for a set of PWMs...}
\description{Lognormal background distribution for a set of PWMs}
\section{Slots}{\describe{\item{\code{bg.source}:}{(\code{\link{character}}) textual description of where the background distribution is derived from}\item{\code{bg.len}:}{(\code{\link{numeric}}) the length to which the background is normalized to. This is a vector of values, 
can have a different value for each motif.}\item{\code{bg.mean}:}{(\code{\link{numeric}}) the mean value of the lognormal distribution at bg.len}\item{\code{bg.sd}:}{(\code{\link{numeric}}) the standard deviation of the lognormal distribution at bg.len}\item{\code{pwms}:}{(\code{\link{list}}) the pwms for which the background has been compiled}}}
\section{Methods}{\describe{\item{\code{\link[=show,PWMLognBackground-method]{show}}}{\code{signature(object = "PWMLognBackground")}: show method for PWMLognBackground}\item{\code{\link[=names,PWMLognBackground-method]{names}}}{\code{signature(x = "PWMLognBackground")}: Name of different pieces of information associated with PWMLognBackground}\item{\code{\link[=$,PWMLognBackground-method]{$}}}{\code{signature(x = "PWMLognBackground")}: Access a property by name}}}

