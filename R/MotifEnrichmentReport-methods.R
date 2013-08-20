# methods to extract data from, and subset MotifEnrichmentReport

#' Columns stored in the motif enrichment report
#'
#' @title Names of variables
#' @name names,MotifEnrichmentReport
#' @aliases names,MotifEnrichmentReport-method
#' @param x the MotifEnrichmentReport object
#' @return the names of the variables
#' @rdname operators-MotifEnrichmentReport
setMethod("names", signature=signature(x="MotifEnrichmentReport"), function(x){
	c(names(x@d), "pwms")
})

#' Access a column by name
#'
#' @aliases $,MotifEnrichmentReport-method
#' @param x the MotifEnrichmentReport object
#' @param name the variable name
#' @rdname operators-MotifEnrichmentReport
setMethod("$", signature=signature(x="MotifEnrichmentReport"), function(x, name){
	if(name == "pwms")
		return(x@pwms)
	else
		return(x@d[[name]])
})

#' Subset the report
#'
#' @aliases [,MotifEnrichmentReport-method
#' @param x the MotifEnrichmentReport object
#' @param i the row selector
#' @param j unused
#' @param ... unused
#' @param drop unused (always FALSE)
#' @rdname operators-MotifEnrichmentReport
setMethod("[", signature=signature(x="MotifEnrichmentReport"), function (x, i, j, ..., drop=TRUE){
	d = x@d[i, , drop=FALSE]
	pwms = x@pwms[i]
	
	new("MotifEnrichmentReport", d=d, pwms=pwms)
})

#' show method for MotifEnrichmentReport
#' @param object the MotifEnrichmentReport object
setMethod("show", signature=signature(object="MotifEnrichmentReport"), function(object){
	d = object@d
	pwms = object@pwms

	cat("An object of class 'MotifEnrichmentReport':\n")	
	if(nrow(d) > 20){
		# show only the first 10 and the last one
		dd = rbind(d[1:10, ], "..."=rep("...", ncol(d)), d[nrow(d),])
		print(dd)
	} else {
		print(d)
	}
})

#' Convert a MotifEnrichmentReport into a data.frame object
#'
#' @name as.data.frame,MotifEnrichmentReport-method
#' @aliases as.data.frame
#' @param x the MotifEnrichmentReport object
#' @param row.names unused
#' @param optional unused
#' @param ... unused
#' @export
setMethod("as.data.frame", signature=signature(x="MotifEnrichmentReport"), function (x, row.names = NULL, optional = FALSE, ...){
	x@d
})



