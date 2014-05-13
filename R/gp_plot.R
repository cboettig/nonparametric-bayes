
#' Plot a Gaussian Process
#' 
#' Plot the mean (line) and 2 standard deviations (ribbon)
#' @param x a gpfit object from \code{\link{gp_fit}}
#' @param true an optional data.frame with columns "x" and "y" containing target curve.  
#' @param ... additional arguments (not implemented) 
#' @return a ggplot object of the plot
#' @export
plot.gpfit <- function(x, true, ...){
  dat <- data.frame(x=x$X, y=(x$Ef), ymin=(x$Ef-2*sqrt(diag(x$Cf))), ymax=(x$Ef+2*sqrt(diag(x$Cf))))
  p <- ggplot(dat) +
    geom_ribbon(aes(x=x,y=y, ymin=ymin, ymax=ymax), fill="grey80") + # Var
    geom_line(aes(x=x,y=y), size=1) + #MEAN
    geom_point(data=gp$obs,aes(x=x,y=y)) + 
    labs(title=paste("llik =", prettyNum(gp$llik)))
  if(!is.null(true))
    p <- p + geom_line(data = true, aes(x, y), col = "red", lty = 2) 
  p
}
