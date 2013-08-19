\name{useBigMemoryPWMEnrich}
\alias{useBigMemoryPWMEnrich}
\title{If to use a faster implementation of motif scanning that requires abount 5 to 10 times more memory}
\usage{
  useBigMemoryPWMEnrich(useBigMemory = FALSE)
}
\arguments{
  \item{useBigMemory}{a boolean value denoting if to use
  big memory implementation}
}
\description{
  If to use a faster implementation of motif scanning that
  requires abount 5 to 10 times more memory
}
\examples{
\dontrun{
useBigMemoryPWMEnrich(TRUE) # switch to big memory implementation globally
useBigMemoryPWMEnrich(FALSE) # switch back to default implementation
}
}

