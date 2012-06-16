\name{makePWMEmpiricalBackground}
\alias{makePWMEmpiricalBackground}
\title{Make an empirical P-value background...}
\usage{makePWMEmpiricalBackground(bg.seq, motifs, bg.pseudo.count=1, bg.source="", verbose=TRUE, ...)
}
\description{Make an empirical P-value background}
\details{Make a background appropriate for empirical P-value calculation. The provided set of background
sequences is contcatenated into a single long sequence which is then scanned with the motifs
and raw scores are saved. This object can be very large. 

For reliable P-value calculation the size of the background set needs to be at least seq.len / min.P.value.
For instance, to get P-values at a resolution of 0.001 for a single sequence of 500bp, we would need 
a background of at least 500/0.001 = 50kb. This ensures that we can make 1000 independent 500bp samples from
this background to properly estimate the P-value. For a group of sequences, we would take seq.len to be the
total length of all sequences in a group.}
\arguments{\item{bg.seq}{a set of background sequences, either a list of DNAString object or DNAStringSet object}
\item{motifs}{a set of motifs, either a list of frequency matrices, or a list of PWM objects. If
frequency matrices are given, the background distribution fitted from bg.seq. 
Same ratios are used for pseudo counts that sum up to bg.pseudo.count for the 4 nucleotides.}
\item{bg.pseudo.count}{the pseudo count which is shared between nucleotides when frequency matrices are given}
\item{bg.source}{a free-form textual description of how the background was generated}
\item{verbose}{if to produce verbose output}
\item{...}{currently unused (this is for convenience for makeBackground function)}
}
\examples{\dontrun{
if(require("PWMEnrich.Dmelanogaster.background")){
data(jaspar.insects.PFM)

# make empirical background by saving raw scores for each bp in the sequence - this can be very large in memory!
if(require("BSgenome.Dmelanogaster.UCSC.dm3")) 
makePWMEmpiricalBackground(Dmelanogaster$upstream2000[1:100], jaspar.insects.PFM)
}
}}

