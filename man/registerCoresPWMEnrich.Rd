\name{registerCoresPWMEnrich}
\alias{registerCoresPWMEnrich}
\title{Register than PWMEnrich can use parallel CPU cores...}
\usage{registerCoresPWMEnrich(numCores=NA)
}
\description{Register than PWMEnrich can use parallel CPU cores}
\details{Certain functions (like motif scanning) can be parallelized in PWMEnrich. This function
registers a number of parallel cores (via core package parallel) to be used in
code that can be parallelized. After this function is called, all further PWMEnrich
function calls will run in parallel if possible. 

By default parallel execution is turned off. To turn it off after using it, call this
function by passing NULL.}
\arguments{\item{numCores}{number of cores to use (default to take all cores), or NULL if no parallel execution is to be used}
}
\examples{\dontrun{
registerCoresPWMEnrich(4) # use 4 CPU cores in PWMEnrich
registerCoresPWMEnrich() # use maximal number of CPUs
registerCoresPWMEnrich(NULL) # do not use parallel execution
}}

