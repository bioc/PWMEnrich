\name{rankingProcessAndReturn}
\alias{rankingProcessAndReturn}
\title{A helper function for motifRankingForGroup and motifRankingForSequence with the common code}
\usage{
  rankingProcessAndReturn(res, r, id, order, rank, unique,
    decreasing)
}
\arguments{
  \item{res}{the list of results from
  MotifEnrichmentResults object}

  \item{r}{the vector of raw results that needs to be
  processed}

  \item{id}{if to return IDs instead of names}

  \item{order}{if to return the ordering of motifs}

  \item{rank}{if to return the rank of motifs}

  \item{unique}{if to remove duplicates}

  \item{decreasing}{specifies the sorting order}
}
\description{
  A helper function for motifRankingForGroup and
  motifRankingForSequence with the common code
}

