
#' Plot a PFM (not PWM) using seqLogo
#'
#' @param pfm a matrix where rows are the four nucleotides
#' @param ... additional parameters for plot()
plotPFM = function(pfm, ...){
	plot(makePWM(divideRows(pfm,colSums(pfm))), ...)
}

#' Plotting for the PWM class
#'
#' This function produces a sequence logo (via package seqLogo). 
#'
#' @aliases plot,PWM,missing-method
#' @param x the PWM object
#' @param y unused
#' @param ... other parameters to pass to seqLogo's \code{plot} function
#' @export
#' @examples
#' if(require("PWMEnrich.Dmelanogaster.background")){
#'   data(MotifDb.Dmel)
#'
#'   # plot the tinman motif from MotifDb
#'   plot(MotifDb.Dmel[["tin"]])
#' }
setMethod("plot", signature=signature(x="PWM", y="missing"), function(x, y, ...){
	pfm = x$pfm
	
	plot(makePWM(divideRows(pfm,colSums(pfm))), ...)	
})

#' Plot mulitple motifs in a single plot
#'
#' Individual motif logos are plotted on a rows x cols grid. This function is a convenience
#' interface for the \code{seqLogoGrid} function that deals with viewpoint placement in a 
#' matrix-like grid layout. 
#'
#' By default will try to make a square grid plot that would fit all the motifs and use
#' list names as captions. 
#'
#' @param pwms a list of PWM objects or frequency matrices
#' @param titles a characater vector of titles for each of the plots
#' @param rows number of rows in the grid
#' @param cols number or cols in the grid
#' @param xmargin.scale the scaling parameter for the X-axis margin. Useful when plotting more than one logo on a page
#' @param ymargin.scale the scaling parameter for the Y-axis margin. Useful when plotting more than one logo on a page
#' @param ... other parameters passed to seqLogoGrid()
#' @export
plotMultipleMotifs = function(pwms, titles=names(pwms), rows=ceiling(sqrt(length(pwms))), 
	cols=ceiling(sqrt(length(pwms))), xmargin.scale=0.4, ymargin.scale=0.4, ...){
	if(!is.list(pwms))
		pwms = list(pwms)
		
	if(length(pwms) != length(titles))
		stop("Number of titles in the 'titles' parameter does not match the number of input motifs")

	# start a new viewport page
	grid.newpage()
	pushViewport(viewport(layout = grid.layout(rows,cols)))
	
	# use the grid layout
	for(i in 1:rows){
		for(j in 1:cols){
			pushViewport(viewport(layout.pos.row = i, layout.pos.col = j))
			 
			# work out which PWM to plot
			inx = (i-1)*cols + j
			# only plot if available
			if( inx <= length(pwms) ){
				pwm = pwms[[inx]]
				if(inherits(pwm, "PWM"))
					pwm = pwm$pfm
				# use the backend functions to plot
				seqLogoGrid(divideRows(pwm,colSums(pwm)), 
					xmargin.scale=xmargin.scale, ymargin.scale=ymargin.scale, title=titles[inx], ...)
			}
			popViewport()
		}
	}
	popViewport()
}

#' Plot the motif enrichment report
#'
#' Plots a graphical version of the motif enrichment report. Note that all values are plotted, if you want to plot only a subset of
#' a report, first select this subset (see examples). 
#'
#' @aliases plot,MotifEnrichmentReport,missing-method
#' @param x a MotifEnrichmentReport object
#' @param y unused
#' @param fontsize font size to use in the plot
#' @param id.fontsize font size to use for the motif IDs
#' @param header.fontsize font size of the header
#' @param widths the relative widths of columns
#' @param ... unused
#' if(require("PWMEnrich.Dmelanogaster.background")){
#'    ###
#'    # load the pre-compiled lognormal background
#'    data(PWMLogn.dm3.MotifDb.Dmel)
#'
#'    # scan two sequences for motif enrichment
#'    sequences = list(DNAString("GAAGTATCAAGTGACCAGTAAGTCCCAGATGA"), DNAString("AGGTAGATAGAACAGTAGGCAATGAAGCCGATG"))
#'
#'    res = motifEnrichment(sequences, PWMLogn.dm3.MotifDb.Dmel)
#'
#'    # produce a report for all sequences taken together
#'    r = groupReport(res)
#'
#'    # plot the top 10 most enriched motifs
#'    plot(r[1:10])
#' 
#' }
setMethod("plot", signature=signature(x="MotifEnrichmentReport", y="missing"), function(x, y, fontsize=14, id.fontsize=fontsize, header.fontsize=fontsize, widths=NULL, ...){
	d = x@d
	pwms = x@pwms
	
	rows = nrow(d)+1
	cols = ncol(d)+1 
	
	# some default widths
	if(is.null(widths))
		widths = c(0.05, 0.1, 0.2, 0.36, 0.08, 0.08, 0.08)[1:(ncol(d)+1)]
		
	widths = widths / sum(widths)
	
	if("z.score" %in% colnames(d)){
		names(d) = c("Rank", "Target", "Motif ID", "Raw score", "Z score", "In top\nmotifs")[1:ncol(d)]
	} else {
		names(d) = c("Rank", "Target", "Motif ID", "Raw score", "P-value", "In top\nmotifs")[1:ncol(d)]
	}
	
	# start a new viewport page
	grid.newpage()
	pushViewport(viewport(layout = grid.layout(rows, cols, widths=widths)))
	
	# use the grid layout
	for(ii in 1:rows){
		for(j in 1:cols){
			pushViewport(viewport(layout.pos.row = ii, layout.pos.col = j))
						
			# figure out which column to plot
			if(j <= 2){ 
				inx = j
			} else {
				inx = j - 1
			}
			
			if(ii == 1){
				# header
				if(j == 3)
					grid.text("PWM", gp=gpar(fontsize=header.fontsize, fontface="bold"))
				else
					grid.text(names(d)[inx], gp=gpar(fontsize=header.fontsize, fontface="bold"))
				
			} else {
				# rest of the table
				i = ii - 1

				# 3rd column is PWM			
				if(j == 3){			
					pwm = pwms[[i]]$pfm
					# use the backend functions to plot
					seqLogoGrid(divideRows(pwm,colSums(pwm)), 
						xmargin.scale=0.01, ymargin.scale=0.01, xaxis=FALSE, yaxis=FALSE, ...)
				} else {
					if(inx == 6){
						grid.text(paste(round(d[i,inx]*100), "%"), gp=gpar(fontsize=fontsize))
					} else if(inx == 3){
						grid.text(d[i,inx], gp=gpar(fontsize=id.fontsize))
					} else if(inx > 3){
						grid.text(signif(d[i,inx],3), gp=gpar(fontsize=fontsize))
					} else {
						grid.text(d[i,inx], gp=gpar(fontsize=fontsize))
					}
				}
			}
				
			popViewport()
		}
	}
	popViewport()
})

#' Draw a motif logo on an existing viewport
#'
#' This function comes from the seqLogo package. It has been modified to remove 
#' some unneccessary code as suggested by W Huber (https://stat.ethz.ch/pipermail/bioconductor/2010-September/035267.html).
#' 
#' Use this function for more advanced plotting where the viewports are directly set up and maintained (see package \code{grid}). 
#'
#' @param pwm numeric The 4xW position weight matrix.
#' @param ic.scale logical If TRUE, the height of each column is proportional to its information
#'                 content. Otherwise, all columns have the same height.
#' @param xaxis logical If TRUE, an X-axis will be plotted.
#' @param yaxis logical If TRUE, a Y-axis will be plotted.
#' @param xfontsize numeric Font size to be used for the X-axis.
#' @param yfontsize numeric Font size to be used for the Y-axis.
#' @param xmargin.scale the scaling parameter for the X-axis margin. Useful when plotting more than one logo on a page
#' @param ymargin.scale the scaling parameter for the Y-axis margin. Useful when plotting more than one logo on a page
#' @param title to be shown on the top
#' @param titlefontsize the fontsize of the title
#' @export
seqLogoGrid <- function(pwm, ic.scale=TRUE, xaxis=TRUE, yaxis=TRUE, xfontsize=10, yfontsize=10, xmargin.scale=1, ymargin.scale=1, title="",
	titlefontsize=15){

  if (inherits(pwm, "pwm")){
    pwm <- pwm@pwm    
  }else if (is.data.frame(pwm)){
    pwm <- as.matrix(pwm)
  }else if (!is.matrix(pwm)){
    stop("pwm must be of class matrix or data.frame")
  }

  if (any(abs(1 - apply(pwm,2,sum)) > 0.01))
    stop("Columns of PWM must add up to 1.0")

  
  chars <- c("A","C","G","T")
  letters <- list(x=NULL,y=NULL,id=NULL,fill=NULL)
  npos <- ncol(pwm)

  
  if (ic.scale){
    ylim <- 2
    ylab <- "Information content"
    facs <- pwm2ic(pwm)
  }else{
    ylim <- 1
    ylab <- "Probability"
    facs <- rep(1, npos)
  }
  
  wt <- 1  
  x.pos <- 0  
  for (j in 1:npos){
    
    column <- pwm[,j]
    hts <- 0.95*column*facs[j]
    letterOrder <- order(hts)
        
    y.pos <- 0    
    for (i in 1:4){
      letter <- chars[letterOrder[i]]
      ht <- hts[letterOrder[i]]
      if (ht>0) letters <- addLetter(letters,letter,x.pos,y.pos,ht,wt)
      y.pos <- y.pos + ht + 0.01
    }
    x.pos <- x.pos + wt
  }

  bottomMargin = ifelse(xaxis, 2 + xfontsize/3.5, 2) * ymargin.scale
  leftMargin = ifelse(yaxis, 2 + yfontsize/3.5, 2) * xmargin.scale
  rightMargin = 2 * xmargin.scale
  topMargin = 2 * ymargin.scale
  
  pushViewport(plotViewport(c(bottomMargin,leftMargin,topMargin,rightMargin)))
  pushViewport(dataViewport(0:ncol(pwm),0:ylim,name="vp1"))
  grid.polygon(x=unit(letters$x,"native"), y=unit(letters$y,"native"),
               id=letters$id,
               gp=gpar(fill=letters$fill,col="transparent"))

  # put in the title
  grid.text(title, y=1.1, vjust=0, gp=gpar(fontsize=titlefontsize))
               
  if (xaxis){
    grid.xaxis(at=seq(0.5,ncol(pwm)-0.5),label=1:ncol(pwm), gp=gpar(fontsize=xfontsize))
    grid.text("Position",y=unit(-3,"lines"), gp=gpar(fontsize=xfontsize))
  }
  if (yaxis){
    grid.yaxis(gp=gpar(fontsize=yfontsize))
    grid.text(ylab,x=unit(-3,"lines"),rot=90, gp=gpar(fontsize=yfontsize))
  }
  popViewport()
  popViewport()
}

#' Plot the raw motifs scores as returned by motifScores()
#' 
#' This function visualises the motif scores for one or more sequences. Sequences are drawn as lines, and scores are plotted
#' as triangles at both sides of the line (corresponding to the two strands). The width of the base of the triangle corresponds to motif width and 
#' the height to the motif \code{log(score)} that is positive and greater than the \code{cutoff} parameter (if specified). All scores
#' have the same y-axis, so the heights of bars are comparable between sequences and motifs.
#' 
#' @param scores the list of motifs scores. Each element of the list is a matrix of scores for one sequences. The columns in the matrix
#'               correspond to different motifs. Each column contains the odds (not log-odds!) scores over both strands. For example, 
#'               for a sequence of length 5, scores for a 3 bp motifs could be: \code{c(0.1, 1, 4, NA, NA, 1, 0.3, 2, NA, NA)}. The first
#'               3 numbers are odds scores starting at first three bases, and the second lot of 3 numbers is the scores starting at the
#'               same positions but with the reverse complement of the motif. The last two values are NA on both strands because we do not
#'               support partial motif hits.
#' @param sel.motifs a vector of motif names. Use this parameter to show the motif hits to only a subset of motifs for which the scores are available.
#' @param seq.names a vector of sequence names to show in the graph. If none specified, the sequences will be named Sequence 1, Sequence 2, ...
#' @param cols a vector of colours to use to colour code motif hits. If none are specified, the current palette will be used. 
#' @param cutoff either a single value, or a vector of values. The values are PWM cutoffs after \code{log.fun} (see below). Only motif scores above these cutoffs will be shown. 
#'               If a single values is specified, it will be used for all PWMs, otherwise the vector needs to specify one cutoff per PWM. 
#' @param log.fun the logarithm function to use to calculate log-odds. By default log2 is used for consistency with Biostrings.               
#' @param main the main title
#' @param legend.space the proportion of horizontal space to reserve for the legend. The default is 30%.
#' @param max.score the maximal log-odds score used to scale all other scores. By default this values is automatically determined, but it can
#'                  also be set manually to make multiple plots comparable. 
#' @param trans the level of transparency. By default 50% transparency to be able to see overlapping binding sites
#' @param text.cex the scaling factor for sequence names
#' @param legend.cex the scaling factor for the legend
#' @param motif.names optional vector of motif names to show instead of those present as column names in \code{scores}
#' @param seq.len.spacing the spacing (in bp units) between the end of the sequence line and the text showing the length in bp
#' @param shape the shape to use to draw motif occurances, valid values are "rectangle" (default), "line" and "triangle"
#' @export
#' @examples
#' if(require("PWMEnrich.Dmelanogaster.background")){
#'    ###
#'    # Load Drosophila PWMs
#'    data(MotifDb.Dmel)
#'
#'    # two sequences of interest
#'    sequences = list(DNAString("GAAGTATCAAGTGACCAGGTGAAGTCCCAGATGA"), DNAString("AGGTAGATAGAACAGTAGGCAATGAAGCCGATG"))
#'
#'    # select the tinman and snail motifs
#'    pwms = MotifDb.Dmel[c("tin", "sna")]
#'
#'    # get the raw score that will be plotted
#'    scores = motifScores(sequences, pwms, raw.scores=TRUE)
#'
#'    # plot the scores in both sequences, green for tin and blue for sna
#'    plotMotifScores(scores, cols=c("green", "blue"))
#'     
#' }
plotMotifScores = function(scores, sel.motifs=NULL, seq.names=NULL, cols=NULL, cutoff=NULL, log.fun=log2, main="", legend.space=0.30, max.score=NULL,
	trans=0.5, text.cex=0.9, legend.cex=0.9, motif.names=NULL, seq.len.spacing=8, shape="rectangle"){
	# subset motifs
	if(!is.null(sel.motifs)){
		scores = lapply(scores, function(x) x[, sel.motifs, drop=FALSE])
	}	
	
	if(length(unique(sapply(scores, ncol)))!=1){
		stop("All elements of the 'scores' list need to have matrices with the same number of columns.")
	}
	
	if(!(shape %in% c("line", "triangle", "rectangle"))){
		stop("'shape' parameter needs to be 'rectangle', 'line' or 'triangle'")
	}
	
	# reverse the scores for plotting order!
	scores = rev(scores)

	# number of sequences
	num.seq = length(scores)
	# maximal sequence length
	seq.len = sapply(scores, nrow)/2
	max.seq.len = max(seq.len)
	
	# find out the length of each motif by the number of NAs
	motif.len = apply(scores[[1]], 2, function(x) sum(is.na(x))/2 + 1)
	num.motifs = length(motif.len)
	
	# threshold the signal
	if(is.null(cutoff)){
		cutoff = rep(0, num.motifs)
	} else if(length(cutoff) == 1) {
		cutoff = rep(cutoff, num.motifs)
	} else if(length(cutoff) != num.motifs){
		stop("The length of 'cutoff' does not match the number of shown motifs")
	}
	
	# do the log of scores
	scores = lapply(scores, log.fun)
	
	# apply the cutoff to the scores
	scores = lapply(scores, function(s){
		for(i in 1:ncol(s)){
			sel = which(s[,i] <= cutoff[i])
			if(length(sel) > 0)
				s[sel,i] = 0
		} 
		s
	})
	
	# the largest score to use to scale all scores
	if(is.null(max.score))
		max.score = max(sapply(scores, function(s) max(s, na.rm=TRUE)))
	
	############# PLOTING #############
	
	if(is.null(cols)){
		pal = palette()
		cols = pal[ (1:num.motifs) %% length(pal) ]
	}
		
	# add transparency
	cols.rgb = col2rgb(cols)
	cols.rgb = rbind(cols.rgb, "alpha"=(1-trans)*255) / 255
	
	for(i in 1:length(cols)){
		cols[i] = rgb(cols.rgb[1, i], cols.rgb[2, i], cols.rgb[3, i], cols.rgb[4, i])
	}
	
	
	# set up the empty plotting area
	ylim = c(0, 2*num.seq)
	xlim = c(0, max.seq.len/(1-legend.space)) # allow for extra space for the legend
	
	par(mar=c(0,0,2,1))
	plot(NULL, xlim=xlim, ylim=ylim, xaxt="n", yaxt="n", ylab="", xlab="", bty="n", main=main)
	
	# plot the lines corresponding to the sequences
	for(i in 1:num.seq){
		y = 1+(i-1)*2
		lines(c(1, seq.len[i]+1), c(y, y))
		# ticks at the dn
		lines(c(1,1), c(y-0.03, y+0.03))
		lines(c(seq.len[i]+1, seq.len[i]+1), c(y-0.03, y+0.03))
	}
	
	
	# now plot the signal
	for(i in 1:num.seq){
		s = scores[[i]]
		
		# the basic y axis
		y = 1+(i-1)*2
		
		# iterate over motifs
		for(j in 1:ncol(s)){
			by.strand = matrix(s[,j], ncol=2)
			for(k in 1:2){
				x.start = which(by.strand[,k] > 0)
				x.end = x.start + motif.len[j]
				
				y.start = rep(y, length(x.start))
				if(k == 1){
					# multiply by 0.8 to leave some space between sequences!
					y.end = y.start + by.strand[x.start,k] / max.score * 0.8
				} else {
					y.end = y.start - by.strand[x.start,k] / max.score * 0.8
				}
				
				for(kk in 1:length(x.start)){
					if(shape == "triangle"){
						polygon(c(x.start[kk], x.end[kk], (x.start[kk]+x.end[kk])/2), c(y.start[kk], y.start[kk], y.end[kk]), col=cols[j], border=cols[j])
					} else if(shape == "rectangle"){
						rect(x.start[kk], y.start[kk], x.end[kk], y.end[kk], col=cols[j], border=cols[j])
					} else {
						xmid = (x.start[kk]+x.end[kk])/2
						rect(xmid-0.5, y.start[kk], xmid+0.5, y.end[kk], col=cols[j], border=cols[j])
					}
				}
			}
		}
	}
	
	# find out max scores for each motif separately
	max.score.motif = sapply(1:num.motifs, function(i){
		max(sapply(scores, function(s) max(s[,i], na.rm=TRUE)))
	})
	
	# set up the legend
	if(is.null(motif.names))
		motif.names = colnames(scores[[1]])
		
	legend = paste(motif.names, " (", round(max.score.motif, 2), " max)", sep="")
	legend("topright", pch=rep(15, num.motifs), col=cols, legend=legend, cex=legend.cex)
	
	# plot sequence names
	if(is.null(seq.names)){
		if(is.null(names(scores)))
			seq.names = paste("Sequence", 1:length(scores))
		else
			seq.names = names(scores)
	}
	
	seq.names = rev(seq.names)
		
	for(i in 1:length(seq.names)){
		y = 1+(i-1)*2 + 1
		text(1, y, seq.names[i], adj=0, cex=text.cex)
	}
	
	# plot the lengths of sequences
	for(i in 1:length(seq.len)){
		y = 1+(i-1)*2
		x = seq.len[i] + seq.len.spacing
		
		text(x, y, paste(seq.len[i], "bp"), adj=0, cex=text.cex*0.8)
	}
}

