\name{PWM}
\alias{PWM}
\alias{PWM-class}
\title{A class that represents a Position Weight Matrix (PWM)}
\description{
  A class that represents a Position Weight Matrix (PWM)
}
\section{Slots}{
  \describe{ \item{\code{id}:}{a systematic ID given to
  this PWM, could include the source, version, etc}
  \item{\code{name}:}{the name of the transcription factor
  (TF) to which the PWM corresponds to}
  \item{\code{pfm}:}{Position Frequency Matrix (PFM) from
  which the PWM is derived}
  \item{\code{prior.params}:}{Defines prior frequencies of
  the four bases (A,C,G,T), a named vector. These will be
  added to individual values for the PFM and at the same
  time used as background probabilities}
  \item{\code{pwm}:}{Final Position Weight Matrix (PWM)
  constructed using prior.params with logarithm base 2} }
}

