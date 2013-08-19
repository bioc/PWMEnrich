\name{logNormPval}
\alias{logNormPval}
\title{Calculate the P-value from lognormal distribution with background of equal length}
\usage{
  logNormPval(scores, seq.len, pwm.len, bg.mean, bg.sd,
    bg.len)
}
\arguments{
  \item{scores}{affinity scores for the PWMs, can contain
  scores for more than one sequence (as rows), P-values are
  extracted separately}

  \item{seq.len}{the length distribution of the sequences}

  \item{pwm.len}{the leggths of PWMs}

  \item{bg.mean}{the mean values from the background for
  PWMs}

  \item{bg.sd}{the sd values from the background}

  \item{bg.len}{the length distribution of the background
  (we currently support only constant length)}
}
\description{
  Calculate the P-value from lognormal distribution with
  background of equal length
}

