\name{PWMEmpiricalBackground}
\alias{PWMEmpiricalBackground}
\alias{PWMEmpiricalBackground-class}
\title{Background for calculating empirical P-values}
\description{
  This object contains raw scores for one very long
  sequence, thus it can be very large.
}
\section{Slots}{
  \describe{ \item{\code{bg.source}:}{textual description
  of where the background distribution is derived from}
  \item{\code{bg.fwd}:}{affinity scores (odds) for the
  forward strand. PWMs as columns}
  \item{\code{bg.rev}:}{affinity scores (odds) for the
  reverse strand. PWMs as columns} \item{\code{pwms}:}{the
  pwms for which the background has been compiled} }
}

