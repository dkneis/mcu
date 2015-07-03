#' Simple local sensitivity analysis
#'
#' Performs a simplistic local sensitivity analysis. Allows for a vector-valued
#' objective function, that is, the sensitivity of multiple criteria with
#' respect to parameters can be tested. See details below.
#'
#' @param fn Objective function. This function must accept as its first argument
#'   a named numeric vector of parameters. It must return a named numeric vector
#'   (which can be of length 1).
#' @param p Data frame specifying the parameters' default values as well as the
#'   test values. There must be four columns: 'name', 'default', 'min', and
#'   'max'.
#' @param sort If \code{TRUE}, the result data frames are sorted to show the
#'   most sensitive parameters at the top.
#' @param silent If \code{TRUE}, diagnostic messages are suppressed.
#' @param ... Additional arguments passed to function \code{fn}.
#'
#' @return A list of data frames. The list length and the element names are
#'   determined by the return vector of \code{fn}. In the common case where
#'   \code{fn} returns a (named) scalar, the list is of length 1.
#'   The data frame at position \eqn{i} has the following columns:
#' \itemize{
#'   \item{\code{name}: } Name of the parameter whose influence is being tested.
#'   \item{\code{pmin}: } Lower test value of the parameter.
#'   \item{\code{pmax}: } Upper test value of the parameter.
#'   \item{\code{pdef}: } Default value of the parameter.
#'   \item{\code{fmin}: } \eqn{i}'th element of the return value of \code{fn},
#'                          when called with the lower test value of the
#'                          particular parameter.
#'   \item{\code{fmax}: } \eqn{i}'th element of the return value of \code{fn},
#'                          when called with the upper test value of the
#'                          particular parameter.
#'   \item{\code{fdef}: } \eqn{i}'th element of the return value of \code{fn},
#'                          when called with the parameter's default. The result
#'                          is identical in all rows.
#'   \item{\code{rsMax}: } Maximum relative sensitivity. For each parameter,
#'     this is \eqn{max(abs(fmax-fdef)/fdef, abs(fmin-fdef)/fdef)}.
#'   \item{\code{mono}: } Logical. \code{TRUE} indicates that the value of
#'     \code{fn}'s \eqn{i}'th element increases or decreases monotonocally, when
#'     the parameter is set to the lower, the default, and the upper limit.
#'     Hence, the lowest and highest function values occurs at/beyond the limits
#'     of the parameter ranges. \code{FALSE} indicates missing monotonicity,
#'     i.e. the lowest or highest function value is supposed to occur within the
#'     tested parameter range.
#'     
#' }
#'
#' @note For a each parameter in \code{p}, function \code{fn} is called
#'   twice. In the two calls, the particular parameter is set to its lower and
#'   upper test value, respectively. All other parameters are fixed at their
#'   defaults. This procedure is applied to every parameter. In total, the
#'   sensitivity test requires \eqn{2 * nrow(p) + 1} evaluations of \code{fn}.
#'   The one additional evaluation is necessary to compute \code{fn} with
#'   default values for all parameters ('base scenario').
#'
#'   Lower and upper test values of the parameters can be constructed, for
#'   example, by multiplying the defaults with, say, 0.9 and 1.1. For this to
#'   work a parameter's default must not be zero.
#'
#' @author David Kneis \email{david.kneis@@tu-dresden.de}
#'
#' @export
#'
#' @examples
#' # Sensitivity of parameters of a linear model
#' obs= data.frame(x=c(1,2), y=c(1,2))
#' model= function(p, x) { p["slope"] * x + p["intercept"] }
#' objfun= function(p, obs) { c(sse= sum((obs$y - model(p, obs$x))^2)) }
#' p= data.frame(
#'   name=c("slope","intercept"),
#'   default= c(1, 0),
#'   min= c(0.5, -1),
#'   max= c(2, 1)
#' )
#' sens(fn=objfun, p=p, obs=obs)
#' 
#' # Like above but for a vector-valued objective function
#' objfun= function(p, obs) { c(sse= sum((obs$y - model(p, obs$x))^2),
#'   mae= sum(abs(obs$y - model(p, obs$x)))) }
#' sens(fn=objfun, p=p, obs=obs)

sens= function(fn, p, sort=TRUE, silent=TRUE, ...) {

  # check inputs
  if (!is.function(fn))
    stop(paste0("'fn' must be a function"))
  if (!is.data.frame(p))
    stop("parameters must be given as a data frame")
  required= c("name","default","min","max")
  if (is.null(names(p)) || (!all(required %in% colnames(p))))
    stop(paste0("missing column names in data frame of parameters,",
      " expecting '",paste(required,collapse="', '"),"'"))
  if (any(p$max < p$min))
    stop("parameter ranges not reasonable")
  if (any(p$default <= p$min) || any(p$default >= p$max))
    stop("parameter defaults must be > minimum and < maximum")

  # base run: all parameters at default values
  if (!silent)
    print("evaluating 'fn' with default parameters")
  p_current= setNames(p$default, p$name)
  tmp= fn(p_current, ...)
  if (is.null(names(tmp)) || (any(names(tmp) == "")))
    stop("'fn' must return a named vector")

  # establish output structure
  out= vector(mode="list", length=length(tmp))
  for (k in 1:length(out)) {
    out[[k]]= data.frame(name=p$name, pmin=p$min, pmax=p$max, pdef=p$default,
      fdef=rep(tmp[k],nrow(p)), fmin=NA, fmax=NA,
      stringsAsFactors=FALSE)
    names(out)[k]= names(tmp)[k]
  }

  # evaluate function for varied parameters
  for (i in 1:nrow(p)) {
    if (!silent)
      print(paste0("testing sensitivity of 'fn' with resp. to '",p$name[i],"'"))
    p_current= setNames(p$default, p$name)

    # test minimum of current parameter and register results
    p_current[i]= p$min[i]
    tmp= fn(p_current, ...)
    for (k in 1:length(out)) {
      out[[k]]$fmin[i]= tmp[k]
    }
    
    # test maximum of current parameter and register results
    p_current[i]= p$max[i]
    tmp= fn(p_current, ...)
    for (k in 1:length(out)) {
      out[[k]]$fmax[i]= tmp[k]
    }
  }

  # compute maximum relative sensitivity
  out= lapply(X=out, FUN=function(x){cbind(x, rsMax=
    pmax( abs((x$fmax - x$fdef) / x$fdef), abs((x$fmin - x$fdef) / x$fdef)))})
  # compute whether the function values changes monotonically over the range
  out= lapply(X=out, FUN=function(x){cbind(x, mono=
    sign(x$fmax - x$fdef) == sign(x$fdef - x$fmin))})

  # sort by sensitivity
  if (sort) {
    out= lapply(X=out, FUN=function(x){x[order(x$rsMax, decreasing=TRUE),]})
  }

  return(out)
}

