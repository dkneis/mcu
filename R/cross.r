#' Cross-tabulation
#'
#' Converts a data frame in long format (e.g. output from a data base) into
#' wide format (cross-table) based on criteria from a \emph{single} column.
#'
#' @param x Data frame to be transformed.
#' @param colsRowHead Columns in \code{x} forming the row headers of the result.
#'   These columns keep their position. They usually contain factor data such
#'   as strings, integers, or logical values.
#' @param colColHead Single column in \code{x} defining the newly created column
#'   headers in the result. This column usually contains factor data such
#'   as strings, integers, or logical values.
#' @param colValues Column with the actual data used to populate the resulting
#'   cross table.
#' @param prefix Character string used as a prefix when constructing
#'   the names of the newly created columns. This is usually required to get
#'   valid column names if the data in column \code{colColHead} are not strings.
#' @param suffix Character string used as a suffix when constructing
#'   the names of the newly created columns.
#' @param sortRows Numeric. If zero, no sorting is performed. Values greater
#'   (less than) zero sort the result in ascending (descending) order based on
#'   the combined columns in \code{colsRowHead}.
#' @param sortCols Numeric. If zero, no sorting is performed. Values greater
#'   (less than) zero sort the newly created columns in ascending (descending)
#'   order based on their names.
#' @param sep A character string used as a separator when collapsing information
#'   from different columns. This string must not appear in any of the involved
#'   columns. An error is generated if this is the case.
#'
#' @return A data frame with \eqn{m + n} columns where \eqn{m} is the 
#'   length of \code{colsRowHead} and \eqn{n} equals the number of unique
#'   values in the column specified as \code{colColHead}.
#'
#' @author david.kneis@@tu-dresden.de
#'
#' @export
#'
#' @examples
#' # Data base format: Meteorological variables at different times and locations
#' times= ISOdate(2000,6,1:10)
#' locs= c("Berlin","London","Rome")
#' vars= c("rain", "temp", "wind")
#' d= data.frame(stringsAsFactors=FALSE,
#'   time= rep(times, length(locs)*length(vars)),
#'   loc= rep(rep(locs, length(vars)), each=length(times)),
#'   var= rep(vars, each=length(times)*length(vars)),
#'   val= NA
#' )
#' d$val[d$var == "rain"]= runif(sum(d$var == "rain")) > 0.5
#' d$val[d$var == "temp"]= rnorm(n=sum(d$var == "temp"), mean=20, sd=5)
#' d$val[d$var == "wind"]= 10^runif(sum(d$var == "wind"), min=-3, max=2)
#' 
#' # Cross table with observed variables in columns
#' cross(x=d, colsRowHead=c("time","loc"), colColHead="var", colValues="val",
#'   prefix="", suffix="", sortRows=0, sortCols=0, sep="\t")
#' 
#' # Cross table with locations in columns
#' cross(x=d, colsRowHead=c("time","var"), colColHead="loc", colValues="val",
#'   prefix="", suffix="", sortRows=0, sortCols=0, sep="\t")
#' 
#' # A common mistake
#' \dontrun{
#'   cross(x=d, colsRowHead=c("time"), colColHead="loc", colValues="val",
#'     prefix="", suffix="", sortRows=0, sortCols=0, sep="\t")
#' }

cross= function(x, colsRowHead, colColHead, colValues, prefix="", suffix="",
  sortRows=0, sortCols=0, sep="\t") {

  # Check arguments
  if (!is.data.frame(x))
    stop("expecting a data frame for 'x'")

  if ((!is.character(colsRowHead)) || (length(colsRowHead) < 1))
    stop("expecting a non-empty character vector for 'colsRowHead'")
  if (!all(colsRowHead %in% names(x)))
    stop("all elements of 'colsRowHead' must match column names of 'x'")

  if ((!is.character(colColHead)) || (length(colColHead) != 1))
    stop("expecting a character string for 'colColHead'")
  if (!colColHead %in% names(x))
    stop("'colColHead' must match a column name of 'x'")

  if ((!is.character(colValues)) || (length(colValues) != 1))
    stop("expecting a character string for 'colValues'")
  if (!colValues %in% names(x))
    stop("'colValues' must match a column name of 'x'")

  if ((!is.character(prefix)) || (length(prefix) != 1))
    stop("expecting a character string for 'prefix'")

  if ((!is.character(suffix)) || (length(suffix) != 1))
    stop("expecting a character string for 'suffix'")

  if ((!is.numeric(sortRows)) || (length(sortRows) != 1))
    stop("expecting a numeric value for 'sortRows'")

  if ((!is.numeric(sortCols)) || (length(sortCols) != 1))
    stop("expecting a numeric value for 'sortCols'")

  if ((!is.character(sep)) || (length(sep) != 1))
    stop("expecting a character string for 'sep'")
  if (any(grepl(pattern=sep, x=unlist(x[,c(colsRowHead, colColHead)]), fixed=TRUE)))
    stop("'sep' detected in data; please chose another value for 'sep'")

  # Generate temporary ID column 
  idCol= paste(names(x), collapse=".")
  check= apply(x[,c(colsRowHead, colColHead)], 1, paste, collapse=sep)
  if (any(duplicated(check)))
    stop(paste0("combination of selected columns ('",
      paste(c(colsRowHead, colColHead), collapse="', '"),"') does not",
      " produce unique record identifiers; you need to add",
      " another column to 'colsRowHead' or remove duplicate cases first",
      " (e.g. by aggregation)"))
  if (length(colsRowHead) == 1) {
    x= cbind(x, x[,colsRowHead])
  } else {
    x= cbind(x, apply(x[,colsRowHead], 1, paste, collapse=sep))
  }
  names(x)[ncol(x)]= idCol

  # Initialize result
  res= x[!duplicated(x[,idCol]), c(idCol, colsRowHead)]

  # Populate new columns
  newCols= unique(x[,colColHead])
  if (sortCols > 0) newCols= sort(newCols)
  if (sortCols < 0) newCols= sort(newCols, decreasing=TRUE)
  for (nc in newCols) {
    res= merge(x=res, y=x[x[,colColHead] == nc, c(idCol, colValues)],
      by=idCol, all=TRUE)
    names(res)[match(colValues, names(res))]= paste0(prefix,nc,suffix)
  }

  # Sort and clean-up
  if (sortRows > 0) res= res[order(res[,idCol]),]
  if (sortCols < 0) res= res[order(res[,idCol], decreasing=TRUE),]
  res[,idCol]= NULL

  return(res)
}

