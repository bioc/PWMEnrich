\name{plot}
\alias{plot}
\alias{plot,MotifEnrichmentReport,missing-method}
\alias{plot,PWM,missing-method}
\title{Plotting for the PWM class}
\arguments{
  \item{x}{the PWM object}

  \item{y}{unused}

  \item{...}{other parameters to pass to seqLogo's
  \code{plot} function}

  \item{x}{a MotifEnrichmentReport object}

  \item{y}{unused}

  \item{fontsize}{font size to use in the plot}

  \item{header.fontsize}{font size of the header}

  \item{widths}{the relative widths of columns}

  \item{...}{unused
  if(require("PWMEnrich.Dmelanogaster.background")){ ### #
  load the pre-compiled lognormal background
  data(PWMLogn.dm3.MotifDb.Dmel)

  # scan two sequences for motif enrichment sequences =
  list(DNAString("GAAGTATCAAGTGACCAGTAAGTCCCAGATGA"),
  DNAString("AGGTAGATAGAACAGTAGGCAATGAAGCCGATG"))

  res = motifEnrichment(sequences,
  PWMLogn.dm3.MotifDb.Dmel)

  # produce a report for all sequences taken together r =
  groupReport(res)

  # plot the top 10 most enriched motifs plot(r[1:10])

  }}
}
\description{
  This function produces a sequence logo (via package
  seqLogo).

  Plots a graphical version of the motif enrichment report.
  Note that all values are plotted, if you want to plot
  only a subset of a report, first select this subset (see
  examples).
}
\examples{
if(require("PWMEnrich.Dmelanogaster.background")){
  data(MotifDb.Dmel)

  # plot the tinman motif from MotifDb
  plot(MotifDb.Dmel$tin)
}
}

