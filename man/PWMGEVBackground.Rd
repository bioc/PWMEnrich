\name{PWMGEVBackground}
\alias{PWMGEVBackground}
\alias{PWMGEVBackground-class}
\title{Generalized Extreme Values (GEV) background for P-values}
\description{
  The three parameters of the GEV distribution are fitted
  by doing linear regression on log of sequence length.
}
\section{Slots}{
  \describe{ \item{\code{bg.source}:}{textual description
  of where the background distribution is derived from}
  \item{\code{bg.loc}:}{linear regression model for
  estimating the location parameter based on log(L), list
  of lm objects of PWMs} \item{\code{bg.scale}:}{linear
  regression model for estimating the scale parameter based
  on log(L), list of lm objects of PWMs}
  \item{\code{bg.shape}:}{linear regression model for
  estimating the shape parameter based on log(L), list of
  lm objects of PWMs} \item{\code{pwms}:}{the pwms for
  which the background has been compiled} }
}

