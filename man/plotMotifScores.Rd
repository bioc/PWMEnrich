\name{plotMotifScores}
\alias{plotMotifScores}
\title{Plot the raw motifs scores as returned by motifScores()}
\usage{
  plotMotifScores(scores, sel.motifs = NULL,
    seq.names = NULL, cols = NULL, cutoff = NULL,
    log.fun = log2, main = "", legend.space = 0.3,
    max.score = NULL, trans = 0.5, text.cex = 0.9,
    legend.cex = 0.9, motif.names = NULL,
    seq.len.spacing = 8, shape="rectangle")
}
\arguments{
  \item{scores}{the list of motifs scores. Each element of
  the list is a matrix of scores for one sequences. The
  columns in the matrix correspond to different motifs.
  Each column contains the odds (not log-odds!) scores over
  both strands. For example, for a sequence of length 5,
  scores for a 3 bp motifs could be: \code{c(0.1, 1, 4, NA,
  NA, 1, 0.3, 2, NA, NA)}. The first 3 numbers are odds
  scores starting at first three bases, and the second lot
  of 3 numbers is the scores starting at the same positions
  but with the reverse complement of the motif. The last
  two values are NA on both strands because we do not
  support partial motif hits.}

  \item{sel.motifs}{a vector of motif names. Use this
  parameter to show the motif hits to only a subset of
  motifs for which the scores are available.}

  \item{seq.names}{a vector of sequence names to show in
  the graph. If none specified, the sequences will be named
  Sequence 1, Sequence 2, ...}

  \item{cols}{a vector of colours to use to colour code
  motif hits. If none are specified, the current palette
  will be used.}

  \item{cutoff}{either a single value, or a vector of
  values. The values are PWM cutoffs after \code{log.fun}
  (see below). Only motif scores above these cutoffs will
  be shown.  If a single values is specified, it will be
  used for all PWMs, otherwise the vector needs to specify
  one cutoff per PWM.}

  \item{log.fun}{the logarithm function to use to calculate
  log-odds. By default log2 is used for consistency with
  Biostrings.}

  \item{main}{the main title}

  \item{legend.space}{the proportion of horizontal space to
  reserve for the legend. The default is 30\%.}

  \item{max.score}{the maximal log-odds score used to scale
  all other scores. By default this values is automatically
  determined, but it can also be set manually to make
  multiple plots comparable.}

  \item{trans}{the level of transparency. By default 50\%
  transparency to be able to see overlapping binding sites}

  \item{text.cex}{the scaling factor for sequence names}

  \item{legend.cex}{the scaling factor for the legend}

  \item{motif.names}{optional vector of motif names to show
  instead of those present as column names in
  \code{scores}}

  \item{seq.len.spacing}{the spacing (in bp units) between
  the end of the sequence line and the text showing the
  length in bp}
  
  \item{shape}{the shape to use to draw motif occurances, 
  valid values are "rectangle" (default), "line" and "triangle"}
}
\description{
  This function visualises the motif scores for one or more
  sequences. Sequences are drawn as lines, and scores are
  plotted as triangles at both sides of the line
  (corresponding to the two strands). The width of the base
  of the triangle corresponds to motif width and the height
  to the motif \code{log(score)} that is positive and
  greater than the \code{cutoff} parameter (if specified).
  All scores have the same y-axis, so the heights of bars
  are comparable between sequences and motifs.
}
\examples{
if(require("PWMEnrich.Dmelanogaster.background")){
   ###
   # Load Drosophila PWMs
   data(MotifDb.Dmel)

   # two sequences of interest
   sequences = list(DNAString("GAAGTATCAAGTGACCAGGTGAAGTCCCAGATGA"), DNAString("AGGTAGATAGAACAGTAGGCAATGAAGCCGATG"))

   # select the tinman and snail motifs
   pwms = MotifDb.Dmel[c("tin", "sna")]

   # get the raw score that will be plotted
   scores = motifScores(sequences, pwms, raw.scores=TRUE)

   # plot the scores in both sequences, green for tin and blue for sna
   plotMotifScores(scores, cols=c("green", "blue"))

}
}

