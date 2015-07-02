#' Plot results of sensitivity analysis
#'
#' Plots the results of a sensitivity analysis carried out using \code{sens}.
#'
#' @param sout A list of data frames returned from a call to \code{\link{sens}}.
#' @param by Character string, either 'fun' or 'par'. If 'fun', the method
#'   creates one plot for each component of the objective function. If
#'   \code{sens} was called with a scalar-valued objective function (and thus
#'   \code{sout} is of length 1) only a single plot is created.
#'   If \code{by} is set to 'par', one plot is created for every parameter whose
#'   influence was tested. This makes sense for vector-valued objective
#'   functions only, i.e. if \code{sout} is of length > 1.
#' @param nc The number of columns used to arrange multiple plots. This info is
#'   passed to the \code{layout} function.   
#'
#' @return \code{NULL}.
#'
#' @note If \code{by='fun'}, the plots can be interpreted as follows:
#'   \itemize{
#'     \item{} The value of the objective function called with default values
#'       for all parameters is plotted as a dot.
#'     \item{} The value of the objective function called with a particular
#'       parameter at the \emph{lower} limit is plotted as a \emph{blue} whisker line having
#'       its origin in the dot.
#'     \item{} The value of the objective function called with a particular
#'       parameter at the \emph{upper} limit is plotted as a \emph{red} whisker line having
#'       its origin in the dot.
#'   }
#'
#'   If \code{by='par'}, the plots can be interpreted as follows:
#'   \itemize{
#'     \item{} The \emph{normalized} value of the objective function called
#'       with default values for all parameters is 1.
#'     \item{} The \emph{normalized} value of the objective function called
#'       with a particular parameter at the \emph{lower} limit is plotted as a \emph{blue}
#'       whisker line (having its origin at x=1).
#'     \item{} The \emph{normalized} value of the objective function called
#'       with a particular parameter at the \emph{upper} limit is plotted as a \emph{red}
#'       whisker line (having its origin at x=1).
#'   }
#'
#' @author David Kneis \email{david.kneis@@tu-dresden.de}
#'
#' @export
#'
#' @examples
#' # Sensitivity of parameters of a linear model
#' obs= data.frame(x=c(1,2), y=c(1,2))
#' model= function(p, x) { p["slope"] * x + p["intercept"] }
#' objfun= function(p, obs) { c(sse= sum((obs$y - model(p, obs$x))^2),
#'   mae= sum(abs(obs$y - model(p, obs$x)))) }
#' p= data.frame(
#'   name=c("slope","intercept"),
#'   default= c(1, 0.1),
#'   min= c(0.5, -1),
#'   max= c(2, 1)
#' )
#' s= sens(fn=objfun, p=p, obs=obs)
#' sensPlot(s, by="fun")
#' sensPlot(s, by="par")

sensPlot= function (sout, by=c("fun","par"), nc=1) {

  # separate plots for each component of the objective function
  if (by[1] == "fun") {
    tempfun= function(x, xnames) {
      rng= range(c(x$fmin, x$fmax), na.rm=TRUE)
      if (!all(is.finite(rng)))
        stop(paste0("only non-finite 'fmin' and/or 'fmax' for criterion '",xnames,"'"))
      plot(x=rng, y=c(1,nrow(x)), type="n", yaxt="n", xlab=xnames, ylab="")
      axis(side=2, at=1:nrow(x), labels=x$name, las=2)
      for (k in 1:nrow(x)) {
        lines(x=c(x$fdef[k], x$fmin[k]), y=c(k,k), col="steelblue2", lwd=3)
        lines(x=c(x$fdef[k], x$fmax[k]), y=c(k,k), col="red")
        points(x$fdef[k],k)
      }
      return(invisible(NULL))
    }
    layout(matrix(1:(ceiling(length(sout)/nc)*nc), ncol=nc, byrow=TRUE))
    omar= par("mar")
    par(mar=c(4,8,2,0.1))
    unused= mapply(FUN=tempfun, x=sout, xnames=names(sout))
    layout(matrix(1))
    par(mar=omar)

  # separate plots for each parameter
  } else if (by[1] == "par") {
    pars= sout[[1]]$name
    layout(matrix(1:(ceiling(length(pars)/nc)*nc), ncol=nc, byrow=TRUE))
    omar= par("mar")
    par(mar=c(4,8,2,0.1))
    for (i in 1:length(pars)) {
      scaledMin= unlist(lapply(sout, function(x, i){x$fmin[i] / x$fdef[i]}, i=i))
      scaledMax= unlist(lapply(sout, function(x, i){x$fmax[i] / x$fdef[i]}, i=i))
      if ((any(!is.finite(scaledMin))) || (any(!is.finite(scaledMax))))
        stop("normalization failed due to zero 'fdef'")
      plot(range(c(scaledMin, scaledMax), na.rm=TRUE), c(1,length(sout)),
        type="n", yaxt="n", xlab="normalized value", ylab="")
      axis(side=2, at=1:length(sout), labels=names(sout), las=2)
      abline(v=1, lty=3)
      for (k in 1:length(scaledMin)) {
        lines(x=c(1, scaledMin[k]), y=c(k,k), col="steelblue2", lwd=3)
        lines(x=c(1, scaledMax[k]), y=c(k,k), col="red")
      }
      mtext(side=3, pars[i], cex=par("cex"))
    }
    layout(matrix(1))
    par(mar=omar)
  } else {
    stop("supplied value of 'how' not supported")
  }
  return(invisible(NULL))
}

