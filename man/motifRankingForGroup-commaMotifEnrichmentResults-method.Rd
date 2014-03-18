\name{motifRankingForGroup,MotifEnrichmentResults-method}
\alias{motifRankingForGroup}
\alias{motifRankingForGroup,MotifEnrichmentResults-method}
\title{Get a ranking of motifs by their enrichment in the whole set of sequences}
\arguments{
  \item{obj}{a MotifEnrichmentResults object}

  \item{bg}{if to use background corrected P-values to do
  the ranking (if available)}

  \item{id}{if to show PWM IDs instead of target TF names}

  \item{order}{if to output the ordering of PWMs instead of
  actual P-values or raw values}

  \item{rank}{if the output should be rank of a PWM instead
  of actual P-values or raw values}

  \item{unique}{if TRUE, only the best rank is taken for
  each TF (only when id = FALSE, order = FALSE)}

  \item{...}{currently unused}
}
\value{
  a vector of log(P-values), P-values or raw enrichments 
  sorted such that the first motif is most enriched
}
\description{
  Get a ranking of motifs by their enrichment in the whole
  set of sequences
}
\examples{
if(require("PWMEnrich.Dmelanogaster.background")){
   ###
   # load the pre-compiled lognormal background
   data(PWMLogn.dm3.MotifDb.Dmel)

   # scan two sequences for motif enrichment
   sequences = list(DNAString("GAAGTATCAAGTGACCAGTAAGTCCCAGATGA"), DNAString("AGGTAGATAGAACAGTAGGCAATGAAGCCGATG"))
   res = motifEnrichment(sequences, PWMLogn.dm3.MotifDb.Dmel)

   # most enriched in both sequences (sorted by lognormal background P-value)
   head(motifRankingForGroup(res))

   # Return a non-redundant set of TFs
   head(motifRankingForGroup(res, unique=TRUE))

   # sorted by raw affinity instead of P-value
   head(motifRankingForGroup(res, bg=FALSE))

   # show IDs instead of target TF names
   head(motifRankingForGroup(res, id=TRUE))

   # output the rank instead of P-value
   head(motifRankingForGroup(res, rank=TRUE))
}
}

