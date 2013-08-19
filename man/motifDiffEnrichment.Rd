\name{motifDiffEnrichment}
\alias{motifDiffEnrichment}
\title{Differential motif enrichment}
\usage{
  motifDiffEnrichment(sequences1, sequences2, pwms,
    score = "autodetect", bg = "autodetect",
    cutoff = log2(exp(4)), verbose = TRUE, res1 = NULL,
    res2 = NULL)
}
\arguments{
  \item{sequences1}{First set of sequences. Can be either a
  single sequence (an object of class DNAString), or a list
  of DNAString objects, or a DNAStringSet object.}

  \item{sequences2}{Second set of sequences. Can be either
  a single sequence (an object of class DNAString), or a
  list of DNAString objects, or a DNAStringSet object.}

  \item{pwms}{this parameter can take multiple values
  depending on the scoring scheme and background correction
  used.  When the \code{method} parameter is set to
  "autodetect", the following default algorithms are going
  to be used: \itemize{ \item if \code{pwms} is a list
  containing either frequency matrices or a list of PWM
  objects then the "affinity" algorithm is selected. If
  frequency matrices are given, they are converted to PWMs
  using uniform background. For best performance, convert
  frequency matrices to PWMs before calling this function
  using realistic genomic background.  \item Otherwise,
  appropriate scoring scheme and background correction are
  selected based on the class of the object (see below).
  }}

  \item{score}{this parameter determines which scoring
  scheme to use. Following scheme as available: \itemize{
  \item "autodetect" - default value. Scoring method is
  determined based on the type of \code{pwms} parameter.
  \item "affinity" - use threshold-free affinity scores
  without a background. The \code{pwms} parameter can
  either be a list of frequency matrices, \code{PWM}
  objects, or a \code{PWMLognBackground} object.  \item
  "cutoff" - use number of motif hits above a score cutoff
  as a measure of enrichment.  No background correction is
  performed. The \code{pwms} parameter can either be a list
  of frequency matrices, \code{PWM} objects, or a
  \code{PWMCutoffBackground} object.  }}

  \item{bg}{this parameter determines which background
  correction to use, if any.  \itemize{ \item "autodetect"
  - default value. Background correction is determined
  based on the type of the \code{pwms} parameter.  \item
  "logn" - use a lognormal distribution background
  pre-computed for a set of PWMs.  This requires
  \code{pwms} to be of class \code{PWMLognBackground}.
  \item "z" - use a z-score for the number of significant
  motif hits compared to background number of hits.  This
  requires \code{pwms} to be of class
  \code{PWMCutoffBackground}.  \item "none" - no background
  correction }}

  \item{cutoff}{the score cutoff for a significant motif
  hit if scoring scheme "cutoff" is selected.}

  \item{res1}{the output of \code{motifEnrichment} if
  already calculated for \code{sequences1}}

  \item{res2}{the output of \code{motifEnrichment} if
  already calculated for \code{sequences2}}

  \item{verbose}{if to produce verbose output}
}
\description{
  Test for differential enrichment between two groups of
  sequences
}
\details{
  This function calls \code{motifEnrichment} on two groups
  of sequences and calculates the difference statistics
  when possible.
}
\examples{
if(require("PWMEnrich.Dmelanogaster.background")){
  # load the background file for drosophila and lognormal correction
  data(PWMLogn.dm3.MotifDb.Dmel)

  # get the differential enrichment
  diff = motifDiffEnrichment(DNAString("TGCATCAAGTGTGTAGTGTGAGATTAGT"), DNAString("TGAACGAGTAGGACGATGAGAGATTGATG"), PWMLogn.dm3.MotifDb.Dmel, verbose=FALSE)

  # motifs differentially enriched in the first sequence (with lognormal background correction)
  head(sort(diff$group.bg, decreasing=TRUE))

  # motifs differentially enriched in the second sequence (with lognormal background correction)
  head(sort(diff$group.bg))
}
}

