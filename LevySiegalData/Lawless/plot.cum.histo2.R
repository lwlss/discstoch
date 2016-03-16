
######################################################################
plot.cum.histo2 <- function(list, save.file=NULL, ylim = c(0, 110), xlim = c(0, 0.6), text.x = 0.15, condition.names=NULL, colors=NULL, LWD=2, tcex=1, print.legend = TRUE,...){
  # Inputs:
  #    well.list -- a list of the colony IDs in each well (output from load.GR $well.list)
  #    colparam -- matrix of colony parameters (output from get.rates $colparam)
  #    save.file -- location of where to save a pdf
  #              -- if NULL, output is to quartz
  #    text.x -- the x position of where to center the text
  #    condition.names -- the names to assign to each well
  #    colors -- a vector of the plotting colors
  
  if (is.null(colors)){
    cl <- rainbow(length(list))
  }else{
    cl <- colors
  }
  if (! is.null(save.file)){
    pdf(file = save.file, height = 6, width = 6)
  }
  i <- 1
  x <- list[[1]]
  x <- x[which(! is.na(x))]
  if (length(x) > 0){
    x <- x[order(x)]
    lx <- length(x)
    y <- 1:lx/lx*100
    plot(x,y, type = "s", lwd =LWD, xlim = xlim, ylim = ylim, ylab = "Cumulative percent", xlab = expression("Relative growth rate"), col=cl[1], ...)
    if(print.legend){
      if (is.null(condition.names)){
        text(text.x, 100, paste(names(list)[1], "  ", signif(mean(x, na.rm=T), 3), "+/-", signif(sd(x, na.rm=T), 2)), col =cl[1], cex= tcex)
      }else{
        text(text.x, 100, paste(condition.names[1], "  ", signif(mean(x, na.rm=T), 3), "+/-", signif(sd(x, na.rm=T), 2)), col =cl[1], cex= tcex)
      }
    }
  }else{
      plot(0,0, type = "s", lwd =LWD, xlim = xlim, ylim = ylim, ylab = "Cumulative percent", xlab =   expression("Relative growth rate"), col=cl[1], ...)
    }
  if (length(list) > 1){
    for (i in 2:length(list)){
      x <- list[[i]]
      if (length(x) > 1){
        x <- x[which(! is.na(x))]
        if (length(x) > 0){
          x <- x[order(x)]
          lx <- length(x)
          y <- 1:lx/lx*100
          points(x,y, type = "s", lwd =LWD, col = cl[i])
          if(print.legend){
            if (is.null(condition.names)){
              text(text.x, 100 - (i-1)*6, paste(names(list)[i], "  ", round(mean(x, na.rm=T), digits=3), "+/-", round(sd(x, na.rm=T), digits=3)), col = cl[i], cex= tcex)
            }else{
              text(text.x, 100 - (i-1)*6, paste(condition.names[i], "  ", round(mean(x, na.rm=T), digits=3), "+/-", round(sd(x, na.rm = T), digits=3)), col = cl[i], cex = tcex)
                 }
          }
        }
      }
    }
  }
    if (! is.null(save.file)){
      dev.off( )
  }
}
######################################################################
