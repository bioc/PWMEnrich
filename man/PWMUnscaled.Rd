\name{PWMUnscaled}
\alias{PWMUnscaled}
\title{Create a PWM from PFM}
\usage{
  PWMUnscaled(x, id = "", name = "",
    type = c("log2probratio", "prob"),
    prior.params = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
    pseudo.count = prior.params, unit.scale = FALSE,
    seq.count = NULL)
}
\arguments{
  \item{x}{the integer count matrix representing the motif,
  rows as nucleotides}

  \item{id}{a systematic ID given to this PWM, could
  include the source, version, etc}

  \item{name}{the name of the transcription factor (TF) to
  which the PWM corresponds to}

  \item{type}{the type of PWM calculation, either as
  log2-odds, or posterior probability (frequency matrix)}

  \item{prior.params}{the pseudocounts for each of the
  nucleotides}

  \item{pseudo.count}{the pseudo-count values if different
  from priors}

  \item{unit.scale}{if to unit.scale the pwm (default is no
  unit scaling)}

  \item{seq.count}{if x is a normalised PFM (i.e. with
  probabilities instead of sequence counts), then this
  sequence count will be used to convert \code{x} into a
  count matrix}
}
\value{
  a new PWM object representing the PWM
}
\description{
  The PWM function from Biostrings without unit scaling
}
\details{
  By default the Biostrings package scales the log-odds
  score so it is within 0 and 1. In this function we take a
  more traditional approach with no unit scaling and offer
  unit scaling as an additional parameter.

  See ?PWM from Biostrings for more information on input
  arguments.
}
\examples{
if(require("PWMEnrich.Dmelanogaster.background")){
    data(MotifDb.Dmel.PFM)

    ttk = MotifDb.Dmel.PFM[["ttk"]]

    PWMUnscaled(ttk, id="ttk-JASPAR", name="ttk") # make a PWM with uniform background
    PWMUnscaled(ttk, id="ttk-JASPAR", name="ttk", prior.params=c("A"=0.2, "C"=0.3, "G"=0.3, "T"=0.2)) # custom background

    prior = getBackgroundFrequencies("dm3", quick=TRUE) # get background for drosophila (quick mode on a reduced dataset)
    PWMUnscaled(ttk, id="ttk-JASPAR", name="ttk", prior.params=prior) # convert using genomic background
}
}

