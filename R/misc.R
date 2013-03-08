#' Calculate medians of columns
#' @param x a matrix
colMedians = function(x){
	apply(x, 2, median, na.rm=T)
}

#' Calculate standard deviations of columns
#' @param x a matrix
colSds = function(x){
	apply(x, 2, sd, na.rm=T)
}

#' Divide each row of a matrix with a vector
#'
#' @param m matrix to be divided
#' @param v the vector to use for division
divideRows = function(m, v){
	t(apply(m, 1, function(x) x/v))
}


#' Convert DNAStringSet to list of DNAString objects
#'
#' as.list doesn't seem to always work for DNAStringSets, so
#' implementing this ourselves. 
#'
#' @param x an object of class DNAStringSet
DNAStringSetToList = function(x){
	res = list()
	for(i in 1:length(x)){
		res[[i]] = x[[i]]
	}
	
	return(res)
}

#' Concatenata DNA sequences into a single character object
#'
#' @param sequences either a list of DNAString objects, or a DNAStringSet
#' @return a single character string
concatenateSequences = function(sequences){
	if(is.list(sequences)){
		paste(unlist(sapply(sequences, toString)), collapse="")
	} else{
		paste(sequences, collapse="")
	}
}
