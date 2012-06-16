\name{PWMCutoffBackground-class}
\alias{PWMCutoffBackground-class}
\alias{PWMCutoffBackground}
\title{Hit count background distribution for a set of PWMs...}
\description{Hit count background distribution for a set of PWMs}
\section{Slots}{\describe{\item{\code{bg.source}:}{(\code{\link{character}}) textual description of where the background distribution is derived from}\item{\code{bg.cutoff}:}{(\code{\link{numeric}}) the cutoff score used to find significant motif hits (in log2 odds),
either a single value or a vector of values}\item{\code{bg.P}:}{(\code{\link{numeric}}) the density of significant motif hits per nucleotide in background}\item{\code{pwms}:}{(\code{\link{list}}) the pwms for which the background has been compiled}}}
\section{Methods}{\describe{\item{\code{\link[=show,PWMCutoffBackground-method]{show}}}{\code{signature(object = "PWMCutoffBackground")}: show method for PWMCutoffBackground}\item{\code{\link[=names,PWMCutoffBackground-method]{names}}}{\code{signature(x = "PWMCutoffBackground")}: Name of different pieces of information associated with PWMCutoffBackground}\item{\code{\link[=$,PWMCutoffBackground-method]{$}}}{\code{signature(x = "PWMCutoffBackground")}: Access a property by name}}}

