\name{PWM-class}
\alias{PWM-class}
\alias{PWM}
\title{A class that represents a Position Weight Matrix (PWM)...}
\description{A class that represents a Position Weight Matrix (PWM)}
\section{Slots}{\describe{\item{\code{id}:}{(\code{\link{character}}) a systematic ID given to this PWM, could include the source, version, etc}\item{\code{name}:}{(\code{\link{character}}) the name of the transcription factor (TF) to which the PWM corresponds to}\item{\code{pfm}:}{(\code{\link{matrix}}) Position Frequency Matrix (PFM) from which the PWM is derived}\item{\code{prior.params}:}{(\code{\link{vector}}) Defines prior frequencies of the four bases (A,C,G,T), a named vector. These will be added 
to individual values for the PFM and at the same time used as background probabilities}\item{\code{pwm}:}{(\code{\link{matrix}}) Final Position Weight Matrix (PWM) constructed using prior.params with logarithm base 2}}}
\section{Methods}{\describe{\item{\code{\link[=plot,PWM,missing-method]{plot}}}{\code{signature(x = "PWM", y = "missing")}: Plotting for the PWM class}\item{\code{\link[=names,PWM-method]{names}}}{\code{signature(x = "PWM")}: Name of different pieces of information associated with PWM}\item{\code{\link[=$,PWM-method]{$}}}{\code{signature(x = "PWM")}: Access a property by name}\item{\code{\link[=length,PWM-method]{length}}}{\code{signature(x = "PWM")}: Length of the motif}\item{\code{\link[=reverseComplement,PWM-method]{reverseComplement}}}{\code{signature(x = "PWM")}: Reverse complement for the PWM object}\item{\code{\link[=show,PWM-method]{show}}}{\code{signature(object = "PWM")}: show method for PWM}}}

