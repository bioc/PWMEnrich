#' Name of different pieces of information associated with PWM
#'
#' @title Names of variables
#' @name names,PWM
#' @aliases names,PWM-method
#' @param x the PWM object
#' @return the names of the variables
#' @rdname operators-PWM
setMethod("names", signature=signature(x="PWM"), function(x) slotNames(x))

#' Access a property by name
#'
#' @aliases $,PWM-method
#' @param x the PWM object
#' @param name the variable name
#' @rdname operators-PWM
setMethod("$", signature=signature(x="PWM"), function(x, name){
	slot(x, name)
})

#' Length of the motif
#'
#' Returns the motif length, i.e. the number of columns in the PWM.
#'
#' @param x the PWM object
#' @aliases length,PWM-method
#' @rdname operators-PWM
setMethod("length", signature=signature(x="PWM"), function(x){
	ncol(x@pwm)
})

#' Reverse complement for the PWM object
#'
#' Finds the reverse complement of the PWM
#'
#' @aliases reverseComplement,PWM-method
#' @param x an object of type PWM
#' @param ... unused
#' @return an object of type PWM that is reverse complement of x
#' @export
#' @examples
#'
#' if(require("PWMEnrich.Dmelanogaster.background")){
#'    data(MotifDb.Dmel.PFM)
#'
#'    reverseComplement(MotifDb.Dmel.PFM$ttk) # reverse complement of the ttk PWM
#' }
#' 
setMethod("reverseComplement", signature=signature(x="PWM"), function (x, ...) {
	pfm = reverseComplement(x@pfm)
	pwm = reverseComplement(x@pwm)
	
	new("PWM", id=paste(x@id, "-- reverse complement"), name=x@name, pfm=pfm, prior.params=x@prior.params, pwm=pwm)
})

#' show method for PWM
#' @param object the PWM object
setMethod("show", signature=signature(object="PWM"), function(object){
	cat("An object of class 'PWM'\n")
	cat("ID:", object$id, "\n")
	cat("Target name:", object$name, "\n")
	cat("Frequency matrix:\n")
	cat("$pfm\n")
	print(object$pfm)
	cat("Position weight matrix (PWM):\n")
	cat("$pwm\n")
	print(object$pwm)
	cat("With background nucleotide frequencies which also serve as pseudo-count:\n")
	cat("$prior.params\n")
	print(object$prior.params)
})

