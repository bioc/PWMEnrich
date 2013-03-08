\name{scanWithPWM}
\alias{scanWithPWM}
\title{Scan the whole sequence on both strands...}
\usage{scanWithPWM(pwm, dna, pwm.rev, odds.score=FALSE, both.strands=FALSE,
    strand.fun="mean")
}
\description{Scan the whole sequence on both strands}
\details{The whole sequence is scanned with a PWM and scores returned beginning at each position. Partial motif
matches are not done, thus the last #[length of motif]-1 scores are NA. 

The function returns either an odds average (*not* log-odds average), maximal score on each strand,
or scores on both strands. 

The function by default returns the score in log2 following the package \code{Biostrings}.}
\value{a vector representing scores starting at each position, or a matrix with score in the two strands}
\arguments{\item{pwm}{PWM object}
\item{dna}{a DNAString or other sequence from Biostrings}
\item{pwm.rev}{the reverse complement for a pwm (if it is already pre-computed)}
\item{odds.score}{if to return raw scores in odds (not logodds) space}
\item{both.strands}{if to return results on both strands}
\item{strand.fun}{which function to use to summarise values over two strands (default is "mean")}
}
\examples{if(require("PWMEnrich.Dmelanogaster.background")){
data(MotifDb.Dmel)

scanWithPWM(MotifDb.Dmel$ttk, DNAString("CGTAGGATAAAGTAACT")) # odds average over the two strands expressed as log2-odds
scanWithPWM(MotifDb.Dmel$ttk, DNAString("CGTAGGATAAAGTAACT"), both.strands=TRUE) # log2-odds scores on both strands
}}

