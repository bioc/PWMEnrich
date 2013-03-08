\name{plotTopMotifsSequence,MotifEnrichmentResults-method}
\alias{plotTopMotifsSequence,MotifEnrichmentResults-method}
\alias{plotTopMotifsSequence}
\title{Plot the top N enrichment motifs in a single sequence...}
\usage{\S4method{plotTopMotifsSequence}{MotifEnrichmentResults}(obj, seq.id, n, bg=TRUE, id=FALSE, ...)
}
\description{Plot the top N enrichment motifs in a single sequence}
\arguments{\item{obj}{a MotifEnrichmentResults object}
\item{seq.id}{either the sequence number or sequence name}
\item{n}{the number of top ranked motifs to plot}
\item{bg}{if to use background P-values to do the ranking (if available)}
\item{id}{if to show PWM IDs instead of target TF names}
\item{...}{other parameters passed to \code{plotMultipleMotifs()}}
}
\examples{\dontrun{
if(require("PWMEnrich.Dmelanogaster.background")){
###
# load the pre-compiled lognormal background
data(PWMLogn.dm3.MotifDb.Dmel)

# scan two sequences for motif enrichment
sequences = list(DNAString("GAAGTATCAAGTGACCAGTAAGTCCCAGATGA"), DNAString("AGGTAGATAGAACAGTAGGCAATGAAGCCGATG"))

res = motifEnrichment(sequences, PWMLogn.dm3.MotifDb.Dmel)

# plot the top 4 motifs in a 2x2 grid
plotTopMotifsSequence(res, 1, 4)

# plot top 3 motifs in a single row
plotTopMotifsSequence(res, 1, 3, row=1, cols=3)
}
}}

