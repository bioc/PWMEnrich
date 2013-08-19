\name{gevPerSequence}
\alias{gevPerSequence}
\title{Apply GEV background normalization per every sequence}
\usage{
  gevPerSequence(scores, seq.len, pwm.len, bg.loc,
    bg.scale, bg.shape)
}
\arguments{
  \item{scores}{affinity scores for the PWMs, can contain
  scores for more than one sequence (as rows), P-values are
  extracted separately}

  \item{seq.len}{the length distribution of the sequences}

  \item{pwm.len}{the lengths of PWMs}

  \item{bg.loc}{list of linear regression for location
  parameter}

  \item{bg.scale}{list of linear regression for scale
  parameter}

  \item{bg.shape}{list of linear regression for shape
  parameter}
}
\description{
  Apply GEV background normalization per every sequence
}

