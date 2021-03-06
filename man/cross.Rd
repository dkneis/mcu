% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/cross.r
\name{cross}
\alias{cross}
\title{Cross-tabulation}
\usage{
cross(x, colsRowHead, colColHead, colValues, prefix = "", suffix = "",
  sortRows = 0, sortCols = 0, sep = "\\t")
}
\arguments{
\item{x}{Data frame to be transformed.}

\item{colsRowHead}{Columns in \code{x} forming the row headers of the result.
These columns keep their position. They usually contain factor data such
as strings, integers, or logical values.}

\item{colColHead}{Single column in \code{x} defining the newly created column
headers in the result. This column usually contains factor data such
as strings, integers, or logical values.}

\item{colValues}{Column with the actual data used to populate the resulting
cross table.}

\item{prefix}{Character string used as a prefix when constructing
the names of the newly created columns. This is usually required to get
valid column names if the data in column \code{colColHead} are not strings.}

\item{suffix}{Character string used as a suffix when constructing
the names of the newly created columns.}

\item{sortRows}{Numeric. If zero, no sorting is performed. Values greater
(less than) zero sort the result in ascending (descending) order based on
the combined columns in \code{colsRowHead}.}

\item{sortCols}{Numeric. If zero, no sorting is performed. Values greater
(less than) zero sort the newly created columns in ascending (descending)
order based on their names.}

\item{sep}{A character string used as a separator when collapsing information
  from different columns. This string must not appear in any of the involved
  columns. An error is generated if this is the case.}
}
\value{
A data frame with \eqn{m + n} columns where \eqn{m} is the
  length of \code{colsRowHead} and \eqn{n} equals the number of unique
  values in the column specified as \code{colColHead}.
}
\description{
Converts a data frame in long format (e.g. output from a data base) into
wide format (cross-table) based on criteria from a \emph{single} column.
}
\examples{
# Data base format: Meteorological variables at different times and locations
times= ISOdate(2000,6,1:10)
locs= c("Berlin","London","Rome")
vars= c("rain", "temp", "wind")
d= data.frame(stringsAsFactors=FALSE,
  time= rep(times, length(locs)*length(vars)),
  loc= rep(rep(locs, length(vars)), each=length(times)),
  var= rep(vars, each=length(times)*length(vars)),
  val= NA
)
d$val[d$var == "rain"]= runif(sum(d$var == "rain")) > 0.5
d$val[d$var == "temp"]= rnorm(n=sum(d$var == "temp"), mean=20, sd=5)
d$val[d$var == "wind"]= 10^runif(sum(d$var == "wind"), min=-3, max=2)

# Cross table with observed variables in columns
cross(x=d, colsRowHead=c("time","loc"), colColHead="var", colValues="val",
  prefix="", suffix="", sortRows=0, sortCols=0, sep="\\t")

# Cross table with locations in columns
cross(x=d, colsRowHead=c("time","var"), colColHead="loc", colValues="val",
  prefix="", suffix="", sortRows=0, sortCols=0, sep="\\t")

# A common mistake
\dontrun{
  cross(x=d, colsRowHead=c("time"), colColHead="loc", colValues="val",
    prefix="", suffix="", sortRows=0, sortCols=0, sep="\\t")
}
}
\author{
david.kneis@tu-dresden.de
}

