# Methods for the different PWMBackground objects

#' show method for PWMLognBackground
#' @param object the PWMLognBackground object
setMethod("show", signature=signature(object="PWMLognBackground"), function(object){
	cat("An object of class '", class(object), "'\n", sep="")
	cat("Background source:", object$bg.source, "\n")
	cat("Fitted on a mean sequence length of", round(mean(object$bg.len)), 
		"for a set of", length(object$pwms), "PWMs\n")
	cat("Lognormal parameters: $bg.mean, $bg.sd\n")
	cat("PWMS: $pwms\n")
})

#' show method for PWMCutoffBackground
#' @param object the PWMCutoffBackground object
setMethod("show", signature=signature(object="PWMCutoffBackground"), function(object){
	cat("An object of class '", class(object), "'\n", sep="")
	cat("Background source:", object$bg.source, "\n")
	cat("Fitted for a set of", length(object$pwms), "PWMs\n")
	cat("Z-score parameters (cutoff is in log2): $bg.cutoff, $bg.P\n")
	cat("PWMS: $pwms\n")
})

#' show method for PWMEmpiricalBackground
#' @param object the PWMEmpiricalBackground object
setMethod("show", signature=signature(object="PWMEmpiricalBackground"), function(object){
	cat("An object of class '", class(object), "'\n", sep="")
	cat("Background source:", object$bg.source, "\n")
	cat("Raw scores for a ", nrow(object$bg.fwd), " bp sequence for ", length(object$pwms), " PWMs\n", sep="")
	cat("Raw scores (odds, not log-odds): $bg.fwd, $bg.rev\n")
	cat("PWMS: $pwms\n")
})

#' show method for PWMGEVBackground
#' @param object the PWMGEVBackground object
setMethod("show", signature=signature(object="PWMGEVBackground"), function(object){
	cat("An object of class '", class(object), "'\n", sep="")
	cat("Background source:", object$bg.source, "\n")
	cat("Generalized extreme value (GEV) distribution fitted for", length(object$pwms), "PWMs\n")
	cat("GEV parameters fitted with linear regressions: $bg.loc, $bg.scale, $bg.shape\n")
	cat("PWMS: $pwms\n")
})


#' Name of different pieces of information associated with PWMLognBackground
#'
#' @title Names of variables
#' @name names,PWMLognBackground
#' @aliases names,PWMLognBackground-method
#' @param x the PWMLognBackground object
#' @return the names of the variables
#' @rdname operators-PWMLognBackground
setMethod("names", signature=signature(x="PWMLognBackground"), function(x) slotNames(x))

#' Access a property by name
#'
#' @aliases $,PWMLognBackground-method
#' @param x the PWMLognBackground object
#' @param name the variable name
#' @rdname operators-PWMLognBackground
setMethod("$", signature=signature(x="PWMLognBackground"), function(x, name){
	slot(x, name)
})

#' Name of different pieces of information associated with PWMCutoffBackground
#'
#' @title Names of variables
#' @name names,PWMCutoffBackground
#' @aliases names,PWMCutoffBackground-method
#' @param x the PWMCutoffBackground object
#' @return the names of the variables
#' @rdname operators-PWMCutoffBackground
setMethod("names", signature=signature(x="PWMCutoffBackground"), function(x) slotNames(x))

#' Access a property by name
#'
#' @aliases $,PWMCutoffBackground-method
#' @param x the PWMCutoffBackground object
#' @param name the variable name
#' @rdname operators-PWMCutoffBackground
setMethod("$", signature=signature(x="PWMCutoffBackground"), function(x, name){
	slot(x, name)
})

#' Name of different pieces of information associated with PWMEmpiricalBackground
#'
#' @title Names of variables
#' @name names,PWMEmpiricalBackground
#' @aliases names,PWMEmpiricalBackground-method
#' @param x the PWMEmpiricalBackground object
#' @return the names of the variables
#' @rdname operators-PWMEmpiricalBackground
setMethod("names", signature=signature(x="PWMEmpiricalBackground"), function(x) slotNames(x))

#' Access a property by name
#'
#' @aliases $,PWMEmpiricalBackground-method
#' @param x the PWMEmpiricalBackground object
#' @param name the variable name
#' @rdname operators-PWMEmpiricalBackground
setMethod("$", signature=signature(x="PWMEmpiricalBackground"), function(x, name){
	slot(x, name)
})

#' Name of different pieces of information associated with PWMGEVBackground
#'
#' @title Names of variables
#' @name names,PWMGEVBackground
#' @aliases names,PWMGEVBackground-method
#' @param x the PWMGEVBackground object
#' @return the names of the variables
#' @rdname operators-PWMGEVBackground
setMethod("names", signature=signature(x="PWMGEVBackground"), function(x) slotNames(x))

#' Access a property by name
#'
#' @aliases $,PWMGEVBackground-method
#' @param x the PWMGEVBackground object
#' @param name the variable name
#' @rdname operators-PWMGEVBackground
setMethod("$", signature=signature(x="PWMGEVBackground"), function(x, name){
	slot(x, name)
})

#' Get the background for a subset of PWMs
#'
#' @aliases [,PWMLognBackground-method
#' @name [,PWMLognBackground-method
#' @param x the PWMLognBackground object
#' @param i the indicies of PWMs
#' @param j unused
#' @param ... unused
#' @param drop unused
#' @rdname subsetting-PWMLognBackground
setMethod("[", "PWMLognBackground", 
function(x, i, j, ..., drop = TRUE){		
	new("PWMLognBackground",
		pwms=x@pwms[i],
		bg.mean=x@bg.mean[i],
		bg.len=x@bg.len[i],
		bg.sd=x@bg.sd[i],
		bg.source=paste(x@bg.source, "--subset")
	)
})

#' Get the background for a subset of PWMs
#'
#' @aliases [,PWMCutoffBackground-method
#' @name [,PWMCutoffBackground-method
#' @param x the PWMCutoffBackground object
#' @param i the indicies of PWMs
#' @param j unused
#' @param ... unused
#' @param drop unused
#' @rdname subsetting-PWMCutoffBackground
setMethod("[", "PWMCutoffBackground", 
function(x, i, j, ..., drop = TRUE){		
	new("PWMCutoffBackground",
		pwms=x@pwms[i],
		bg.cutoff=x@bg.cutoff[i],
		bg.P=x@bg.P[i],
		bg.source=paste(x@bg.source, "--subset")
	)
})

#' Get the background for a subset of PWMs
#'
#' @aliases [,PWMEmpiricalBackground-method
#' @name [,PWMEmpiricalBackground-method
#' @param x the PWMEmpiricalBackground object
#' @param i the indicies of PWMs
#' @param j unused
#' @param ... unused
#' @param drop unused
#' @rdname subsetting-PWMEmpiricalBackground
setMethod("[", "PWMEmpiricalBackground",
function(x, i, j, ..., drop = TRUE)
{		
	new("PWMEmpiricalBackground",
		pwms=x@pwms[i],
		bg.fwd=x@bg.fwd[,i,drop=FALSE],
		bg.rev=x@bg.rev[,i,drop=FALSE],
		bg.source=paste(x@bg.source, "--subset")
	)
})

#' Get the background for a subset of PWMs
#'
#' @aliases [,PWMGEVBackground-method
#' @name [,PWMGEVBackground-method
#' @param x the PWMGEVBackground object
#' @param i the indicies of PWMs
#' @param j unused
#' @param ... unused
#' @param drop unused
#' @rdname subsetting-PWMGEVBackground
setMethod("[", "PWMGEVBackground",
function(x, i, j, ..., drop = TRUE)
{		
	new("PWMGEVBackground",
		pwms=x@pwms[i],
		bg.loc=x@bg.loc[i],
		bg.scale=x@bg.scale[i],
		bg.shape=x@bg.shape[i],
		bg.source=paste(x@bg.source, "--subset")
	)
})



