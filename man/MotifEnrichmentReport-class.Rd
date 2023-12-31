% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllDataClasses.R
\docType{class}
\name{MotifEnrichmentReport-class}
\alias{MotifEnrichmentReport-class}
\title{A report class with formatted results of motif enrichment}
\description{
The columns stored in this object will depend on the type of the report (either for group of sequences, or individual sequences).
}
\section{Slots}{
 
 \describe{
   \item{\code{d}:}{a DataFrame object that contains the main tabular report data}
   \item{\code{pwms}:}{a list of \code{PWM} objects corresponding to rows of \code{d}}
 }
}

