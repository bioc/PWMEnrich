\name{makePWMCutoffBackground}
\alias{makePWMCutoffBackground}
\title{Make a cutoff background...}
\usage{makePWMCutoffBackground(bg.seq, motifs, cutoff=log2(exp(4)), bg.pseudo.count=1, bg.source="",
    verbose=TRUE)
}
\description{Make a cutoff background}
\details{Make a background based on number of motifs hits above a certain threshold.}
\arguments{\item{bg.seq}{a set of background sequences, either a list of DNAString object or DNAStringSet object}
\item{motifs}{a set of motifs, either a list of frequency matrices, or a list of PWM objects. If
frequency matrices are given, the background distribution fitted from bg.seq. 
Same ratios are used for pseudo counts that sum up to bg.pseudo.count for the 4 nucleotides.}
\item{cutoff}{the cutoff at which the background should be made, i.e. at which a motif hit is called significant}
\item{bg.pseudo.count}{the pseudo count which is shared between nucleotides when frequency matrices are given}
\item{bg.source}{a free-form textual description of how the background was generated}
\item{verbose}{if to produce verbose output}
}
\examples{\dontrun{
if(require("PWMEnrich.Dmelanogaster.background")){
data(jaspar.insects.PFM)

# make background for JASPAR motifs using 2kb promoters of all D. melanogaster transcripts using cutoff of 5
if(require("BSgenome.Dmelanogaster.UCSC.dm3")) 
makePWMCutoffBackground(Dmelanogaster$upstream2000, jaspar.insects.PFM, cutoff=log2(exp(5)))
}
}}

