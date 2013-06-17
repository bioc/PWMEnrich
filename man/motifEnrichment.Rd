\name{motifEnrichment}
\alias{motifEnrichment}
\title{Motif enrichment}
\usage{motifEnrichment(sequences, pwms, score="autodetect", bg="autodetect", cutoff,
    verbose=TRUE, motif.shuffles=30, B=1000, group.only=FALSE)
}
\description{Calculate motif enrichment using one of available scoring algorithms and background corrections.}
\details{This function provides and interface to all algorithms available in PWMEnrich
to find motif enrichment in a single or a group of sequences with/without
background correction. 

Since for all algorithms the first step involves calculating raw scores without background correction, the output
always contains the scores without background correction together with (optional) background-corrected
scores. 

Unless otherwise specified the scores are returned both separately for each sequence (without/with background) and
for the whole group of sequences (without/with background). 

To use a background correction you need to supply a set of PWMs with precompiled background distribution parameters
(see function \code{\link{makeBackground}}). When such an object is supplied as the \code{pwm} parameter, the scoring 
scheme and background correction are automatically determined. 

There are additional packages with already pre-computed background (e.g. see package \code{PWMEnrich.Dmelanogaster.background}).

Please refer to (Stojnic & Adryan, 2012) for more details on the algorithms.}
\value{a MotifEnrichmentResults object containing a subset following elements:
\itemize{
\item "score" - scoring scheme used
\item "bg" - background correction used
\item "params" - any additional parameters
\item "sequences" - the set of sequences used
\item "pwms" - the set of pwms used
\item "sequence.nobg" - per-sequence scores without any background correction. 
For "affinity" and "clover" a matrix of mean affinity scores; for
"cutoff" number of significant hits above a cutoff 
\item "sequence.bg" - per-sequence scores after background correction. For "logn" and "pval" the P-value (smaller is better);
for "z" and "ms" background corrections the z-scores (bigger is better). 
\item "group.nobg" - aggregate scores for the whole group of sequences without background correction. For "affinity"
and "clover" the mean affinity over all sequences in the set; for "cutoff" the total number of hits in all
sequences.
\item "group.bg" - aggregate scores for the whole group of sequences with background correction. For "logn" and "pval",
the P-value for the whole group (smaller is better); for "z" and "ms" the z-score for the whole set (bigger is better).
\item "sequence.norm" - (only for "logn") the length-normalized scores for each of the sequences. Currently only implemented
for "logn", where it returns the values normalized from LogN(0,1) distribution
\item "group.norm" - (only for "logn") similar to sequence.norm, but for the whole group of sequences
}}
\references{\itemize{ 
\item R. Stojnic & B. Adryan: Identification of functional DNA motifs using a binding affinity lognormal background distribution, submitted.
\item MC Frith et al: Detection of functional DNA motifs via statistical over-representation, Nucleid Acid Research (2004).
}}
\arguments{\item{sequences}{the sequences to be scanned for enrichment. Can be either a single sequence 
(an object of class DNAString), or a list of DNAString objects, or a DNAStringSet object.}
\item{pwms}{this parameter can take multiple values depending on the scoring scheme and background correction used.
When the \code{method} parameter is set to "autodetect", the following default algorithms
are going to be used: 
\itemize{
\item if \code{pwms} is a list containing either frequency matrices or a list of PWM objects then
the "affinity" algorithm is selected. If frequency matrices are given, they are converted 
to PWMs using uniform background. For best performance, convert frequency matrices to PWMs
before calling this function using realistic genomic background. 
\item Otherwise, appropriate scoring scheme and background correction are selected based on the
class of the object (see below). 
}}
\item{score}{this parameter determines which scoring scheme to use. Following scheme as available:
\itemize{
\item "autodetect" - default value. Scoring method is determined based
on the type of \code{pwms} parameter. 
\item "affinity" - use threshold-free affinity scores without a background. The \code{pwms}
parameter can either be a list of frequency matrices, \code{PWM} objects, or a 
\code{PWMLognBackground} object. 
\item "cutoff" - use number of motif hits above a score cutoff as a measure of enrichment. 
No background correction is performed. The \code{pwms}
parameter can either be a list of frequency matrices, \code{PWM} objects, or a 
\code{PWMCutoffBackground} object.
\item "clover" - use the Clover algorithm (Frith et al, 2004). The Clover score of a single 
sequence is identical to the affinity score, while for a group of sequences is an
average of products of affinities over all sequence subsets. 
}}
\item{bg}{this parameter determines which background correction to use, if any. 
\itemize{
\item "autodetect" - default value. Background correction is determined based on the type
of the \code{pwms} parameter. 
\item "logn" - use a lognormal distribution background pre-computed for a set of PWMs.
This requires \code{pwms} to be of class \code{PWMLognBackground}.
\item "z" - use a z-score for the number of significant motif hits compared to background number of hits.
This requires \code{pwms} to be of class \code{PWMCutoffBackground}.
\item "pval" - use empirical P-value based on a set of background sequences. This requires
\code{pwms} to be of class \code{PWMEmpiricalBackground}. Note that PWMEmpiricalBackground
objects tend to be very large so that the empirical P-value can be calculated in reasonable time.
\item "ms" - shuffle columns of motif matrices and use that as basis for P-value calculation. Note that
since the sequences need to rescanned with all of the new shuffled motifs this can be very slow. 
Also, this also works only no *individual* sequences, not groups. 
\item "none" - no background correction
}}
\item{cutoff}{the score cutoff for a significant motif hit if scoring scheme "cutoff" is selected.}
\item{verbose}{if to print verbose output}
\item{motif.shuffles}{number of times to shuffle motifs if using "ms" background correction}
\item{B}{number of replicates when calculating empirical P-value}
\item{group.only}{if to return statistics only for the group of sequences, not individual sequences. In the case of
empirical background the P-values for individual sequences are not calculated (thus saving time), for other 
backgrounds they are calculated but not returned.}
}
\examples{if(require("PWMEnrich.Dmelanogaster.background")){
###
# load the pre-compiled lognormal background
data(PWMLogn.dm3.MotifDb.Dmel)

# scan two sequences for motif enrichment
sequences = list(DNAString("GAAGTATCAAGTGACCAGTAGATTGAAGTAGACCAGTC"), DNAString("AGGTAGATAGAACAGTAGGCAATGGGGGAAATTGAGAGTC"))
res = motifEnrichment(sequences, PWMLogn.dm3.MotifDb.Dmel)

# most enriched in both sequences (lognormal background P-value)
head(motifRankingForGroup(res))

# most enriched in both sequences (raw affinity, no background)
head(motifRankingForGroup(res, bg=FALSE))

# most enriched in the first sequence (lognormal background P-value)
head(motifRankingForSequence(res, 1))

# most enriched in the first sequence (raw affinity, no background)
head(motifRankingForSequence(res, 1, bg=FALSE))

###
# Load the pre-compiled background for hit-based motif counts with cutoff of P-value = 0.001 
data(PWMPvalueCutoff1e3.dm3.MotifDb.Dmel)

res.count = motifEnrichment(sequences, PWMPvalueCutoff1e3.dm3.MotifDb.Dmel)

# Enrichment in the whole group, z-score for the number of motif hits
head(motifRankingForGroup(res))

# First sequence, sorted by number of motif hits with P-value < 0.001
head(motifRankingForSequence(res, 1, bg=FALSE))

}}

