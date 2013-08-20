\name{sequenceReport,MotifEnrichmentResults-method}
\alias{sequenceReport}
\alias{sequenceReport,MotifEnrichmentResults-method}
\title{Generate a motif enrichment report for a single sequence}
\arguments{
  \item{obj}{a MotifEnrichmentResults object}

  \item{seq.id}{the sequence index or name}

  \item{bg}{if to use background corrected P-values to do
  the ranking (if available)}

  \item{...}{unused}
}
\value{
  a MotifEnrichmentReport object containing a table with
  the following columns: \itemize{ \item 'rank' - The rank
  of the PWM's enrichment in the sequence \item 'target' -
  The name of the PWM's target gene, transcript or protein
  complex.  \item 'id' - The unique identifier of the PWM
  (if set during PWM creation).  \item 'raw.score' - The
  raw score before P-value calculation \item 'p.value' -
  The P-value of motif enrichment (if available) }
}
\description{
  Generate a motif enrichment report for a single sequence
}
\examples{
if(require("PWMEnrich.Dmelanogaster.background")){
   ###
   # load the pre-compiled lognormal background
   data(PWMLogn.dm3.MotifDb.Dmel)

   # scan two sequences for motif enrichment
   sequences = list(DNAString("GAAGTATCAAGTGACCAGTAAGTCCCAGATGA"), DNAString("AGGTAGATAGAACAGTAGGCAATGAAGCCGATG"))

   res = motifEnrichment(sequences, PWMLogn.dm3.MotifDb.Dmel)

   # reports for the two sequences
   r1 = sequenceReport(res, 1)
   r2 = sequenceReport(res, 2)

   # view the results
   r1
   r2

   # plot the top 10 most enriched motifs in the first, and then second sequence
   plot(r1[1:10])
   plot(r2[1:10])

}
}

