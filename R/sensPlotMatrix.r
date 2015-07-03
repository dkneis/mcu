# Helper function
# low: value of obj. function for lower limit of parameter
# mid: value of obj. function for default of parameter
# upp: value of obj. function for upper limit of parameter
# thr: threshold value for relative sensitivity
symb_scalar= function(low, mid, upp, thr) {
  insensitive= list(pch=1, col="black", bg="black")
  unknown= list(pch=4, col="black", bg="white")
  middle_def= list(pch=23, col="grey", bg="grey")
  middle_unc= list(pch=5, col="grey", bg="white")
  lower_def= list(pch=25, col="steelblue", bg="steelblue")
  lower_unc= list(pch=6, col="steelblue", bg="white")
  upper_def= list(pch=24, col="salmon", bg="salmon")
  upper_unc= list(pch=2, col="salmon", bg="white")
  v= c(low, mid, upp)
  if (all(is.finite(v))) {
    if (length(unique(v)) == 1) {
      res= insensitive
    } else if (min(v) == mid) {
      maxRelSens= min(mid/upp, upp/mid, mid/low, low/mid, na.rm=TRUE)
      if (!is.finite(maxRelSens))
        stop("can't compute maximum relative sensitivity")
      if (maxRelSens <= thr) { res= middle_def } else { res= middle_unc }
    } else if (min(v) == low) {
      maxRelSens= min(low/mid, mid/low, low/upp, upp/low, na.rm=TRUE)
      if (!is.finite(maxRelSens))
        stop("can't compute maximum relative sensitivity")
      if (maxRelSens <= thr) { res= lower_def } else { res= lower_unc }
    } else if (min(v) == upp) {
      maxRelSens= min(upp/mid, mid/upp, upp/low, low/upp, na.rm=TRUE)
      if (!is.finite(maxRelSens))
        stop("can't compute maximum relative sensitivity")
      if (maxRelSens <= thr) { res= upper_def } else { res= upper_unc }
    }
  } else {
    res= unknown
  }
  return(res)
}

# Helper function; as above but for vectors 
symb_vector= function(low, mid, upp, thr) {
  if (!identical(length(low), length(mid), length(upp)))
    stop("input lengths differ")
  return(mapply(FUN=symb_scalar, low=low, mid=mid, upp=upp, MoreArgs=list(thr=thr)))
}

#' Plot sensitivity info as a matrix
#'
#' Plots the results of a sensitivity analysis carried out using \code{sens}
#' as a matrix of symbols. The output can facilitate manual model calibration
#' (see details below).
#'
#' @param sout A list of data frames returned from a call to \code{\link{sens}}.
#' @param thr A number in range \eqn{0 < thr < 1}. This is used to pragmatically
#'   distinguish between highly sensitive and less sensitive parameters.
#'   \code{thr} specifies the minimum relative decrease in the value of the
#'   objective function to consider a particular varied parameter as
#'   \emph{highly} sensitive. A value of, say, 0.9 means that a parameter is
#'   consired as highly sensitive if \eqn{ftst/fdef <= 0.9}, where \eqn{ftst} 
#'   and \eqn{fdef} denote the output of the objective function for a test value
#'   and the default value of a parameter, respectively. Reasonable values are
#'   probably between 0.8 and 0.95.
#' @param xpars Logical. Controls the plot's layout. If \code{TRUE}, the
#'   parameter names appear as column headers and the objective function(s) as
#'   row headers(s). If \code{FALSE}, the result matrix is transposed.
#'
#' @return \code{NULL}.
#'
#' @note Symbols in the created graphics have the following meaning:
#'   \itemize{
#'     \item{} Triangles pointing up: Parameter should be increased to/beyond
#'       the tested upper limit in order to reduce the value of the objective
#'       function. Filled/non-filled: Relative sensitivity is high/low.
#'     \item{} Triangles pointing down: Parameter should be decreased to/beyond
#'       the tested lower limit in order to reduce the value of the objective
#'       function. Filled/non-filled: Relative sensitivity is high/low.
#'     \item{} Diamond: Lowest value of the objective function does not
#'       (exclusively) occur at a boundary of the tested parameter range. Hence,
#'       the optimum parameter values may be inside that range. Filled/non-filled:
#'       Relative sensitivity is high/low.
#'     \item{} Circle: Objective function is not sensitive to the parameter.
#'     \item{} Cross: No information due to non-finite return values of the
#'       objective function.
#'   }
#'
#'   Note that the analysis does not account for possible parameter interactions
#'   such as compensation effects. For the example (see below), the created plot
#'   suggests that the intercept should be decreased (although the true optimum
#'   value is 0). This is due to the fact that all test values for the slope
#'   are actually too high (true optimum at 1).
#'
#'   In cases with long names of parameters/functions, it will be necessary to
#'   adjust the plot margings accordingly (see example).
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
#'   default= c(1.5, 0.1),
#'   min= c(1.1, -1),
#'   max= c(2, 1)
#' )
#' s= sens(fn=objfun, p=p, obs=obs)
#' omar= par("mar")
#' par(mar=c(0.5,6,10,0.5))
#' sensPlotMatrix(sout=s, thr=0.75, xpars=TRUE)
#' par(mar=omar)

sensPlotMatrix= function (sout, thr=0.9, xpars=TRUE) {
  nr= ifelse(xpars, length(sout), nrow(sout[[1]]))
  nc= ifelse(xpars, nrow(sout[[1]]), length(sout))
  plot(c(.5, nc+0.5), c(.5, nr+0.5), type="n", xaxt="n", yaxt="n", xlab="", ylab="")
  fnames= names(sout)
  pnames= sout[[1]]$name
  axis(side=2, at=1:nr, labels={if (xpars) fnames else pnames}, las=2, tick=FALSE)
  axis(side=3, at=1:nc, labels={if (xpars) pnames else fnames}, las=2, tick=FALSE)
  for (i in 1:length(fnames)) {
    e= sout[[fnames[i]]]
    n= match(pnames, e$name)
    sym= symb_vector(low=e$fmin[n], mid=e$fdef[n], upp=e$fmax[n], thr=thr)
    x= { if (xpars) 1:length(n) else rep(i, length(n)) }
    y= { if (xpars) rep(i, length(n)) else 1:length(n) }
    points(x=x, y=y, type="p", pch=unlist(sym["pch",]),
      col=unlist(sym["col",]), bg=unlist(sym["bg",]))
  }
  return(invisible(NULL))
}


