\name{PWMEmpiricalBackground-class}
\alias{PWMEmpiricalBackground-class}
\alias{PWMEmpiricalBackground}
\title{Background for calculating empirical P-values...}
\description{Background for calculating empirical P-values}
\details{This object contains raw scores for one very long sequence, thus it can be very large.}
\section{Slots}{\describe{\item{\code{bg.source}:}{(\code{\link{character}}) textual description of where the background distribution is derived from}\item{\code{bg.fwd}:}{(\code{\link{matrix}}) affinity scores (odds) for the forward strand. PWMs as columns.}\item{\code{bg.rev}:}{(\code{\link{matrix}}) affinity scores (odds) for the reverse strand. PWMs as columns.}\item{\code{pwms}:}{(\code{\link{list}}) the pwms for which the background has been compiled}}}
\section{Methods}{\describe{\item{\code{\link[=show,PWMEmpiricalBackground-method]{show}}}{\code{signature(object = "PWMEmpiricalBackground")}: show method for PWMEmpiricalBackground}\item{\code{\link[=names,PWMEmpiricalBackground-method]{names}}}{\code{signature(x = "PWMEmpiricalBackground")}: Name of different pieces of information associated with PWMEmpiricalBackground}\item{\code{\link[=$,PWMEmpiricalBackground-method]{$}}}{\code{signature(x = "PWMEmpiricalBackground")}: Access a property by name}}}

