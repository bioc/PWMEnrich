\name{matrixShuffleZscorePerSequence}
\alias{matrixShuffleZscorePerSequence}
\title{Obtain z-score for motif column shuffling...}
\usage{matrixShuffleZscorePerSequence(scores, sequences, pwms, cutoff, B=30)
}
\description{Obtain z-score for motif column shuffling}
\details{All PWMs are shuffled at the same time. This function would be too slow to produce
empirical P-values, thus we return a z-score from a small number of shuffles. 

The z-scores are calculated for each sequence individually.}
\arguments{\item{scores}{a set of already calculated scores}
\item{sequences}{either one sequence or a list/set of sequences (objects of type DNAString or DNAStringSet)}
\item{pwms}{a list of PWMs}
\item{cutoff}{if NULL, will use affinity, otherwise will use number of hits over this log2 odds cutoff}
\item{B}{number of replicates, i.e. PWM column shuffles}
}

