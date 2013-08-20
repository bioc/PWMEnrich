\name{groupReport,MotifEnrichmentResults-method}
\alias{groupReport}
\alias{groupReport,MotifEnrichmentResults-method}
\title{Generate a motif enrichment report for the whole group of sequences together}
\arguments{
  \item{obj}{a MotifEnrichmentResults object}

  \item{top}{what proportion of top motifs should be
  examined in each individual sequence (by default 5\%)}

  \item{bg}{if to use background corrected P-values to do
  the ranking (if available)}

  \item{...}{unused}
}
\value{
  a MotifEnrichmentReport object containing a table with
  the following columns: \itemize{ \item 'rank' - The rank
  of the PWM's enrichment in the whole group of sequences
  together \item 'target' - The name of the PWM's target
  gene, transcript or protein complex.  \item 'id' - The
  unique identifier of the PWM (if set during PWM
  creation).  \item 'raw.score' - The raw score before
  P-value calculation \item 'p.value' - The P-value of
  motif enrichment (if available) \item 'top.motif.prop' -
  The proportion (between 0 and 1) of sequences where the
  motif is within \code{top} proportion of enrichment
  motifs. }
}
\description{
  Generate a motif enrichment report for the whole group of
  sequences together
}
\examples{
if(require("PWMEnrich.Dmelanogaster.background")){
   ###
   # load the pre-compiled lognormal background
   data(PWMLogn.dm3.MotifDb.Dmel)

   # scan two sequences for motif enrichment
   sequences = list(DNAString("GAAGTATCAAGTGACCAGTAAGTCCCAGATGA"), DNAString("AGGTAGATAGAACAGTAGGCAATGAAGCCGATG"))

   res = motifEnrichment(sequences, PWMLogn.dm3.MotifDb.Dmel)

   # produce a report for all sequences taken together
   r.default = groupReport(res)

   # produce a report where the last column takes top 1\\% motifs
   r = groupReport(res, top=0.01)

   # view the results
   r

   # plot the top 10 most enriched motifs
   plot(r[1:10])

}
}

