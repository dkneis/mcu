#' Dynamically dimensioned search
#'
#' An inplementation of the dynamically dimensioned search (DDS) algorithm for
#' stochastic optimization of computationally expensive objective functions. The
#' algorithm gradually shifts from global to local search.
#'
#' @param fn Objective function. This function must accept as its first argument
#'   a named numeric vector of parameters. It must return a named numeric vector
#'   (which can be of length 1). The algorithm attempts to minimize \code{fn}
#'   with respect to the \emph{first element} of its return vector, i.e.
#'   additional elements are fro diagnostic purposes only.
#' @param p Data frame with initial values and bounds for all parameters. The
#'   expected column names are 'name', 'initial', 'min', and 'max'.
#' @param m Requested number of function evaluations (integer). This is the
#'   only stopping criterion of the algorithm.
#' @param r Numeric algorithm parameter (default 0.2). It controls the jump
#'   length during generation of new parameters sets.
#' @param plot Logical. It \code{TRUE} the value of \code{fn} (first element
#'   only) is plotted after each evaluation.
#' @param ... Additional arguments passed to function \code{fn}.
#'
#' @return A list with the following elements
#' \itemize{
#'   \item{f_best} Named numeric vector. Value of the objective function at the
#'     found minimum (with respect to the vector's first element).
#'   \item{p_best} Named numeric vector. Parameters corresponding to \code{f_best}.
#'   \item{f_trace} Data frame with \code{m} rows. Holds the 'best' value of
#'     the objective function for all iterations.
#'   \item{p_trace} Data frame with \code{m} rows. Holds the 'best' parameters
#'     for all iterations.
#' }
#'
#' @author David Kneis \email{david.kneis@@tu-dresden.de}
#'
#' @references The original algorithm is described in
#'   Tolson, B. A. and Shoemaker, C. A. (2007): Dynamically
#'   dimensioned search algorithm for computationally efficient watershed model
#'   calibration, Water Resources Research, 43, W01413,
#'   doi:10.1029/2005WR004723.
#'
#' @export
#'
#' @examples
#' # Fitting coefficients of an ordinary linear model
#' nobs= 100
#' obs= data.frame(x=1:nobs, y=1:nobs + (1:nobs)*(runif(nobs)-0.5))
#' model= function(p, x) { p["slope"] * x + p["intercept"] }
#' objfun= function(p, obs) { c(sse= sum((obs$y - model(p, obs$x))^2)) }
#' p= data.frame(
#'   name=c("slope","intercept"),
#'   min= c(-1, -50),
#'   max= c(4, 50),
#'   initial= 0, 20
#' )
#' opt= dds(fn=objfun, p=p, m=250, r=0.2, obs=obs)
#' layout(matrix(c(1,2,3,4), ncol=2, byrow=TRUE))
#' plot(obs)
#' lines(range(obs$x), model(p=opt$p_best, x=range(obs$x)))
#' # trace of function
#' plot(0:(nrow(opt$f_trace)-1), opt$f_trace$sse, type="l",
#'   xlab="iter",ylab="sse")
#' # trace of pars
#' plot(0:(nrow(opt$p_trace)-1), opt$p_trace$slope, type="l",
#'   xlab="iter",ylab="slope")
#' plot(0:(nrow(opt$p_trace)-1), opt$p_trace$intercept, type="l",
#'   xlab="iter",ylab="intercept")
#' layout(matrix(1))


dds= function(fn, p, m=10, r=0.2, plot=FALSE, ...) {

  # Step 1: Check inputs
  if (!is.function(fn))
    stop(paste0("'fn' must be a function"))
  if (!is.data.frame(p))
    stop("parameters must be given as a data frame")
  required= c("name","min","max","initial")
  if (is.null(names(p)) || (!all(required %in% colnames(p))))
    stop(paste0("missing column names in data frame of parameters,",
      " expecting '",paste(required,collapse="', '"),"'"))
  if (any(p$max < p$min) || any(p$max < p$default) || any(p$min > p$default))
    stop("parameter ranges not reasonable")


  # Step 2: Initial evaluation of objective function
  i= 1
  p_best= setNames(p$initial, p$name)
  f_best= fn(p_best, ...)
  if (any(is.na(f_best)))
    stop("'fn' returns NA for initial parameter set")
  if (is.null(names(f_best)) || (any(names(f_best) == "")))
    stop("'fn' must return a named vector")

  # Save initial state
  p_trace= data.frame(matrix(p_best, nrow=1))
  names(p_trace)= names(p_best)
  f_trace= data.frame(matrix(f_best, nrow=1))
  names(f_trace)= names(f_best)

  # Iteration loop
  while (i < (m+1)) {

    # Step 3: Select parameters for perturbation
    inds_perturb= which( runif(nrow(p)) <= (1 - log(i) / log(m)) )
    if (length(inds_perturb) == 0)
      inds_perturb= sample.int(n=nrow(p), size=1)

    # Step 4: Perturbation with reflection
    p_new= p_best
    p_new[inds_perturb]= p_best[inds_perturb] + 
      r * (p$max[inds_perturb] - p$min[inds_perturb]) *
      rnorm(n=length(inds_perturb), mean=0, sd=1)
    # Reflection at min
    inds_out= inds_perturb[ p_new[inds_perturb] < p$min[inds_perturb] ]
    if (length(inds_out) > 0) {
      p_new[inds_out]= p$min[inds_out] + (p$min[inds_out] - p_new[inds_out])
      inds_out= inds_out[ p_new[inds_out] > p$max[inds_out] ]
      if (length(inds_out) > 0) {
        p_new[inds_out]= p$min[inds_out]
      }
    }
    # Reflection at max
    inds_out= inds_perturb[ p_new[inds_perturb] > p$max[inds_perturb] ]
    if (length(inds_out) > 0) {
      p_new[inds_out]= p$max[inds_out] - (p_new[inds_out] - p$max[inds_out])
      inds_out= inds_out[ p_new[inds_out] < p$min[inds_out] ]
      if (length(inds_out) > 0) {
        p_new[inds_out]= p$max[inds_out]
      }
    }

    # Step 5: Evaluate objective function and possibly update best solution
    f_new= fn(p_new, ...)
    if (!is.na(f_new[1])) {
      if (f_new[1] <= f_best[1]) {
        f_best= f_new
        p_best= p_new
      }
    }

    # Step 6: Update counter
    i= i + 1

    # Save current state
    p_trace= rbind(p_trace, p_best)
    f_trace= rbind(f_trace, f_best)

    if (plot) {
      plot(1:nrow(f_trace), f_trace[,1], type="l",
        xlab="Iteration",ylab="Obj. function (1st elem.)")
    }

  }
  return(list(
    f_best=f_best, p_best=p_best, f_trace=f_trace, p_trace=p_trace
  ))
}

