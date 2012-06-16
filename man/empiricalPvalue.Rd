\name{empiricalPvalue}
\alias{empiricalPvalue}
\title{Calculate the empirical P-value by affinity of cutoff.}
\usage{empiricalPvalue(scores, seq.len, pwm.len, bg.fwd, bg.rev, cutoff, B=10000,
    verbose=FALSE, exact.length=FALSE)
}
\description{Calculate the empirical P-value by affinity of cutoff.}
\details{This is the new backend function for empirical P-values for either affinity or cutoff. 
The function only works on single sequences.}
\arguments{\item{scores}{the scores obtained for the sequence}
\item{seq.len}{the length of the sequence, if a single value will take a single sequence
of given length. If a vector of values, will take sequences of given lengths
and joint them together}
\item{pwm.len}{the lengths of PWMs}
\item{bg.fwd}{raw odds scores for the forward strand of background}
\item{bg.rev}{raw odds scores for the reverse strand of background}
\item{cutoff}{if not NULL, will use hit count above this cutoff. The cutoff should be specified in log2.}
\item{B}{the number of random replicates}
\item{verbose}{if to give verbose progress reports}
\item{exact.length}{if to take into consideration that the actual sequence lengths differ for different PWMs.
For very long sequences (i.e. seq.len >> pwm.len) this make very little difference, however
the run time with exact.length is much longer.}
}

