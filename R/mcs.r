#' Monte-Carlo simulation
#'
#' Performs a Monte-Carlo simulation with parameters sampled from a uniform
#' distribution using a latin hypercube method.
#'
#' @param fn Function of interest. It must accept as its first argument
#'   a named numeric vector of parameters. It must return a named numeric vector
#'   (which can be of length 1).
#' @param p Data frame specifying ranges and defaults for the varied parameters.
#'   There must be three columns: 'name', 'default', 'min', and 'max'.
#' @param nRuns Desired number of parameter samples. The total number of
#'   evaluations of \code{fn} is \code{nRuns} + 1.
#' @param silent If \code{TRUE}, diagnostic messages are suppressed.
#' @param ... Additional arguments passed to function \code{fn}.
#'
#' @return A list of two elements \code{p} and \code{out}, both being data
#'   frames with \code{nRuns} + 1 rows. \code{p} holds the tested parameter values
#'   with column names taken from the 'name' field of \code{ranges}.
#'   \code{out} holds the return values of \code{fn}. Each row in \code{out}
#'   corresponds to the same row of \code{p}. In the common case where \code{fn}
#'   returns a scalar result, \code{out} contains just a single column.
#'   The first row in both data frames corresponds to the default parameter set.
#'
#' @author David Kneis \email{david.kneis@@tu-dresden.de}
#'
#' @export
#'
#' @examples
#' # Analysis of the residuals' sum of squares for a linear model
#' obs= data.frame(x=c(1,2), y=c(1,2))
#' model= function(p, x) { p["slope"] * x + p["intercept"] }
#' objfun= function(p, obs) { c(sse= sum((obs$y - model(p, obs$x))^2)) }
#' p= data.frame(
#'   name=c("slope","intercept"),
#'   default= c(1, 0),
#'   min= c(0.5, -1),
#'   max= c(2, 1)
#' )
#' x= mcs(fn=objfun, p=p, obs=obs)
#' layout(matrix(1:2, ncol=2))
#' plot(x$p[,"slope"], x$out$sse, xlab="slope", ylab="SSE")
#' plot(x$p[,"intercept"], x$out$sse, xlab="intercept", ylab="SSE")
#' layout(matrix(1))

mcs= function(fn, p, nRuns=10, silent=TRUE, ...) {

  # Check inputs
  if (!is.function(fn))
    stop(paste0("'fn' must be a function"))
  if (!is.data.frame(p))
    stop("'p' must be a data frame")
  required= c("name","default","min","max")
  if (is.null(names(p)) || (!all(required %in% colnames(p))))
    stop(paste0("missing column names in data frame 'p',",
      " expecting '",paste(required,collapse="', '"),"'"))
  if (any(p$max < p$min))
    stop("parameter ranges not reasonable")

  # Sample parameters
  if (!silent)
    print("creating sample")
  tmp= improvedLHS(n=nRuns, k=nrow(p))
  colnames(tmp)= p$name
  prand= as.data.frame(tmp)
  for (i in 1:nrow(p)) {
    prand[,i]= p[i,"min"] + prand[,i] * (p[i,"max"] - p[i,"min"])
  }

  # Simulation with defaults to initialize result table
  if (!silent)
    print(paste0("initial run with defaults"))
  tmp= fn(setNames(p$default,p$name),...)
  if (is.null(names(tmp)) || (any(names(tmp) == "")))
    stop("'fn' does not return a named vector")
  out= data.frame(matrix(tmp, nrow=1))
  names(out)= names(tmp)

  # Simulations
  for (i in 1:nRuns) {
    if (!silent)
      print(paste0("run ",i," of ",nRuns," (",(i-1)/nRuns*100,"% completed)"))
    tryCatch({
      tmp= fn(setNames(unlist(prand[i,]),names(prand)),...)
      out= rbind(out,tmp)
    }, error= function(e) {
      print(e)
      out= rbind(out,rep(NA, ncol(out)))
    }, warning= function(w) {
      print(w)
      out= rbind(out,rep(NA, ncol(out)))
    })
  }

  # Add default parameters to table of parameters
  prand= rbind(as.list(setNames(p$default, p$name)), prand)

  # Return tested parameter values and function results
  return(list(p=prand, out=out))
}

