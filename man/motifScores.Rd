\name{motifScores}
\alias{motifScores}
\title{Motif affinity of number of hits over a threshold}
\usage{
  motifScores(sequences, motifs, raw.scores = FALSE,
    verbose = TRUE, cutoff = NULL)
}
\arguments{
  \item{sequences}{a set of sequences to be scanned, a list
  of DNAString or other scannable objects}

  \item{motifs}{a list of motifs either as frequency
  matrices (PFM) or as PWM objects. If PFMs are specified
  they are converted to PWMs using uniform background.}

  \item{raw.scores}{if to return raw scores (odds) for each
  position in the sequence. Note that scores for forward
  and reverse strand are concatenated into a single long
  vector of scores (twice the length of the sequence)}

  \item{verbose}{if to print verbose output}

  \item{cutoff}{if not NULL, will count number of matches
  with score above value specified (instead of returning
  the average affinity).  Can either be one value, or a
  vector of values for each of the motifs.}
}
\value{
  if raw.scores=FALSE, returns a matrix of mean scores
  (after cutoff if any), where columns are motifs.  The
  returned values are either mean odd scores (not log-odd),
  or number of hits above a threshold; otherwise if
  raw.scores=TRUE, returns a list of raw score values
  (before cutoff)
}
\description{
  Scan a number of sequences either to find overall
  affinity, or a number of hits over a score threshold.
}
\examples{
if(require("PWMEnrich.Dmelanogaster.background")){
   data(MotifDb.Dmel)

   affinity = motifScores(DNAString("CGTAGGATAAAGTAACTAGTTGATGATGAAAG"), MotifDb.Dmel) # affinity scores
   counts = motifScores(DNAString("CGTAGGATAAAGTAACTAGTTGATGATGAAAG"), MotifDb.Dmel, cutoff=log2(exp(4))) # motif hit count with Patser score of 4
   print(affinity)
   print(counts)

   # scanning multiple sequences
   sequences = list(DNAString("CGTAGGATAAAGTAACTAGTTGATGATGAAAG"), DNAString("TGAGACGAAGGGGATGAGATGCGGAAGAGTGAAA"))
   affinity2 = motifScores(sequences, MotifDb.Dmel)
   print(affinity2)
}
}

