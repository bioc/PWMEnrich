\name{reverseComplement}
\alias{reverseComplement}
\alias{reverseComplement,PWM-method}
\title{Reverse complement for the PWM object}
\arguments{
  \item{x}{an object of type PWM}

  \item{...}{unused}
}
\value{
  an object of type PWM that is reverse complement of x
}
\description{
  Finds the reverse complement of the PWM
}
\examples{
if(require("PWMEnrich.Dmelanogaster.background")){
   data(MotifDb.Dmel.PFM)

   reverseComplement(MotifDb.Dmel.PFM[["Dmelanogaster-JASPAR_CORE-ttk-MA0460.1"]]) # reverse complement of the ttk PWM
}
}

