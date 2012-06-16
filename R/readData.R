## load in motifs 

#' Read in motifs in JASPAR or TRANSFAC format
#'
#' The format is autodetected based on file format. If the autodetection fail
#' then the file cannot be read. 
#'
#' @param file the filename
#' @param remove.acc if to remove accession numbers. If TRUE, the AC entry in TRANSFAC files is ignored,
#'                   and the accession is stripped from JASPAR, e.g. motif with name "MA0211.1 bap" would 
#'                   become just "bap". If FALSE, botht he AC and ID are used to generate the TRANSFAC name
#'                   and the original motif names are preserved in JASPAR files. 
#' @return a list of 4xL matrices representing motifs (four nucleotides as rows)
#' @export
#' @examples
#'
#' # read in example TRANSFAC motifs without accession codes (just IDs)
#' readMotifs(system.file(package="PWMEnrich", dir="extdata", file="example.transfac"), remove.acc=TRUE)
#'
#' # read in the JASPAR insects motifs provided as example
#' readMotifs(system.file(package="PWMEnrich", dir="extdata", file="jaspar-insecta.jaspar"), remove.acc=TRUE)               
readMotifs = function(file, remove.acc=FALSE){
	t = readLines(file)
	
	if(length(t) == 0){
		stop(paste("Non-existing or empty file", file))
	}
	
	is.transfac = length(grep("^XX", t)) > 0
	is.jaspar = length(grep("^>", t)) > 0
	
	if((is.transfac & is.jaspar) | (!is.transfac & !is.jaspar)){
		stop("Cannot detect the file format of supplied file. Please look at the JASPAR and TRANSFAC websites for examples of well-formed motif files.")
	}
	
	if(is.jaspar){
		return(readJASPAR(file, remove.acc))
	} else {
		return(readTRANSFAC(file, remove.acc)) 
	}
}


#' Read motifs in JASPAR format
#'
#' @param file the filename
#' @param remove.ids if to strip JASPAR ID's from motif names, e.g. "MA0211.1 bap" would become just "bap"
#' @return a list of matrices representing motifs (with four nucleotides as rows)
readJASPAR = function(file, remove.ids=FALSE){
	t = readLines(file)
	
	h = grep("^>", t)
	
	# ordering of nucleotides
	norder = c("A", "C", "G", "T")

	motifs = list()

	for(i in 1:length(h)){
		# header and motif from jaspar format
		header = t[h[i]]
		h.inx = (h[i]+1)
		motif = t[h.inx:(h.inx+3)]
	
		# motif in clover format
		motif.matrix = t(sapply(strsplit(motif, "[ \t\\[]+"), function(x) as.integer(x[2:(length(x)-1)])))
		nucleotides = sapply(strsplit(motif, "[ \t\\[]+"), function(x) x[1])
		# if the ordering is not the same, re-order into ACGT order
		motif.matrix = motif.matrix[match(norder, nucleotides), ]
		rownames(motif.matrix) = norder
		
		motifs[[i]] = motif.matrix
		motif.name = trim(gsub(">", "", header))
		if(remove.ids){
			motif.parts = unlist(strsplit(motif.name, " +"))
			if(length(motif.parts)>1){
				motif.name = paste(motif.parts[2:length(motif.parts)], collapse="_")
			}
		}
		names(motifs)[i] = motif.name
	
	}
	
	motifs
}

#' Read in motifs in TRANSFAC format
#'
#' @param file the filename
#' @param remove.acc if to ignore transfac accession numbers
#' @return a list of matrices representing motifs (with four nucleotides as rows)
readTRANSFAC = function(file, remove.acc=TRUE){
	t = readLines(file)
	
	# ordering of nucleotides
	norder = c("A", "C", "G", "T")
	
	# get indicies of different tags
	ac.inx = grep("^AC", t)
	id.inx = grep("^ID", t)
	po.inx = grep("^P[O0]", t)
	xx.inx = grep("^XX", t)
	
	if(length(ac.inx) != length(id.inx) || length(ac.inx) != length(po.inx) || !all(ac.inx < id.inx) || !all(id.inx < po.inx)){
		stop("Inconsistent TRANSFAC format. Every motif needs to have exactly one AC, ID and PO/P0 entry (in that order).")
	}

	# extract individual motifs	
	motifs = list()
	for(i in 1:length(ac.inx)){
		ac = trim(gsub("^AC", "", t[ac.inx[i]]))
		# find the id entry that is just after the ac entry
		id.cur = min( id.inx[id.inx > ac.inx[i]] )
		id = trim(gsub("^ID", "", t[id.cur]))
		# find boundaries of the motif
		po.cur = min( po.inx[po.inx > ac.inx[i]] )
		xx.cur = min( xx.inx[xx.inx > po.cur] )
		
		# check the ordering of nucleotides
		po = trim(gsub("^P[O0]", "", t[po.cur]))
		po = unlist(strsplit(po, "[ \t]+"))
		if(!identical(po, norder)){
			stop(paste("In all motifs the nucleotides need to be ordered: A C G T. This is not the case for motif", id))
		}
		
		# extract motif
		m = t[(po.cur+1) : (xx.cur-1)]
		m = strsplit(m, "[ \t]+")
		m = sapply(m, as.integer)[-1,]		
		rownames(m) = norder
		
		# generate the name and save
		if(remove.acc)
			name = id
		else
			name = paste(ac, id, sep="_")
		
		motifs[[name]] = m
	}
	
	motifs

}
