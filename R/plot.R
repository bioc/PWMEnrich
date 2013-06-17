
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
#' @param x the PWM object
#' @param y unused
#' @param ... other parameters to pass to seqLogo's \code{plot} function
#' @export
#' @examples
#' if(require("PWMEnrich.Dmelanogaster.background")){
#'   data(MotifDb.Dmel)
#'
#'   # plot the tinman motif from MotifDb
#'   plot(MotifDb.Dmel$tin)
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
		stop("Provide titles for each of the input motifs")

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
				if(class(pwm) == "PWM")
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

  if (class(pwm) == "pwm"){
    pwm <- pwm@pwm    
  }else if (class(pwm) == "data.frame"){
    pwm <- as.matrix(pwm)
  }else if (class(pwm) != "matrix"){
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

