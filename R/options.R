# managing parallel execution of code

# this global variable records options, e.g. for parallel execution etc
.PWMEnrich.Options = new.env(parent=emptyenv())

#' Register than PWMEnrich can use parallel CPU cores
#'
#' Certain functions (like motif scanning) can be parallelized in PWMEnrich. This function
#' registers a number of parallel cores (via core package parallel) to be used in
#' code that can be parallelized. After this function is called, all further PWMEnrich
#' function calls will run in parallel if possible. 
#'
#' By default parallel execution is turned off. To turn it off after using it, call this
#' function by passing NULL. 
#'
#' @param numCores number of cores to use (default to take all cores), or NULL if no parallel execution is to be used
#' @export
#' @examples
#' \dontrun{
#' registerCoresPWMEnrich(4) # use 4 CPU cores in PWMEnrich
#' registerCoresPWMEnrich() # use maximal number of CPUs
#' registerCoresPWMEnrich(NULL) # do not use parallel execution
#' }
registerCoresPWMEnrich = function(numCores=NA){
	if (!require("parallel"))
	    stop("Parallel execution requires package parallel")
	
	if(!is.null(numCores) && is.na(numCores))
		numCores = detectCores()
	
	assign("numCores", numCores, pos=.PWMEnrich.Options)
}

#' If to use a faster implementation of motif scanning that requires abount 5 to 10 times more memory
#'
#' @param useBigMemory a boolean value denoting if to use big memory implementation
#'
#' @export
#' @examples
#' \dontrun{
#' useBigMemoryPWMEnrich(TRUE) # switch to big memory implementation globally
#' useBigMemoryPWMEnrich(FALSE) # switch back to default implementation
#' }
useBigMemoryPWMEnrich = function(useBigMemory=FALSE){
	assign("useBigMemory", useBigMemory, pos=.PWMEnrich.Options)
}
