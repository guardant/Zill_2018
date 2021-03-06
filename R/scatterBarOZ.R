#' Draw a scatter-bar plot
#' 
#' Takes in a 3-column data frame (sample number, clonal fraction, color) and plots it.
#' @param x 3-column data frame to be plotted (sample number, clonal fraction, color) 
#' @param lhist length passed to length.out
#' @return Generates a scatter bar plot of x
#' @export
scatterBarOZ <- function(x, lhist=20, ...){
    ## check input
    stopifnot(ncol(x)==3)
    ## set up layout and graphical parameters
    layMat <- matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
    layout(layMat, widths=c(5/7, 2/7), heights=c(2/7, 5/7))
    ospc <- 0.5 # outer space
    pext <- 4 # par extension down and to the left
    bspc <- 1 # space between scatter plot and bar plots
    par. <- par(mar=c(pext, pext, bspc, bspc),
                oma=rep(ospc, 4)) # plot parameters
    ## scatter plot
    plot(x[,1:2], xlim=range(x[,1]), ylim=range(x[,2]), col=x[,3], ...)
	plot.new()  #skips the top histogram
    ## 3) determine barplot and height parameter
    ## histogram (for barplot-ting the density)
    #xhist <- hist(x[,1], plot=FALSE, breaks=seq(from=min(x[,1]), to=max(x[,1]),
    #                                 length.out=lhist))
    yhist <- hist(x[,2], plot=FALSE, breaks=seq(from=min(x[,2]), to=max(x[,2]),
                                     length.out=lhist)) # note: this uses probability=TRUE
    ## determine the plot range and all the things needed for the barplots and lines
    #xx <- seq(min(x[,1]), max(x[,1]), length.out=num.dnorm) # evaluation points for the overlaid density
    #xy <- dnorm(xx, mean=mean(x[,1]), sd=sd(x[,1])) # density points
    #yx <- seq(min(x[,2]), max(x[,2]), length.out=num.dnorm)
    #yy <- dnorm(yx, mean=mean(x[,2]), sd=sd(x[,2]))
    ## barplot and line for x (top)
    #par(mar=c(0, pext, 0, 0))
    #barplot(xhist$density, axes=FALSE, ylim=c(0, max(xhist$density, xy)),
    #        space=0) # barplot
    #lines(seq(from=0, to=lhist-1, length.out=num.dnorm), xy, col=dcol) # line
    ## barplot and line for y (right)
    par(mar=c(pext, 0, 0, 0))
    barplot(yhist$density, axes=FALSE, xlim=c(0, max(yhist$density)),
            space=0, horiz=TRUE) # barplot
    #lines(yy, seq(from=0, to=lhist-1, length.out=num.dnorm), col=dcol) # line
    ## restore parameters
    par(par.)
}

#' example:
#require(mvtnorm)
#X <- rmvnorm(1000, c(0,0), matrix(c(1, 0.8, 0.8, 1), 2, 2))
#scatterBarNorm(X, xlab=expression(italic(X[1])), ylab=expression(italic(X[2])))
