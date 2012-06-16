# managing parallel execution of code

# this global variable records options for parallel execution of code
.PWMEnrich.Parallel = new.env(parent=emptyenv())

#' Register than PWMEnrich can use parallel CPU cores
#'
#' Certain functions (like motif scanning) can be parallelized in PWMEnrich. This function
#' registers a number of parallel cores (via packages doMC and foreach) to be used in
#' code that can be parallelized. After this function is called, all further PWMEnrich
#' function calls will run in parallel if possible. 
#'
#' By default parallel execution is turned off. To turn it off after using it, call this
#' function without any argument or by passing NULL. 
#'
#' @param numCores number of cores to use, or NULL if no parallel execution is to be used
#' @export
#' @examples
#' \dontrun{
#' registerCoresPWMEnrich(4) # use 4 CPU cores in PWMEnrich
#' }
registerCoresPWMEnrich = function(numCores=NULL){
	if (!require("foreach") | !require("doMC"))
	    stop("Parallel execution requires package foreach and doMC")
	
	registerDoMC(cores=numCores)
	
	assign("numCores", numCores, pos=.PWMEnrich.Parallel)
}
