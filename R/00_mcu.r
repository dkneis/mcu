#' Model calibration utilities
#'
#' Utilities to aid in the calibration of models, especially those
#' with many parameters and significant computation times.
#'
#'
#' @name mcu-package
#' @aliases mcu
#' @docType package
{}

################################################################################
# Required packages

# Package for latin hypercube sampling
#' @import lhs
if (!("lhs" %in% installed.packages()[,1]))
  install.packages("lhs")
library("lhs")

# Packages for parallel execution
#' @import foreach
if (!("foreach" %in% installed.packages()[,1]))
  install.packages("foreach")
library("foreach")

# Packages for parallel execution
#' @import doParallel
if (!("doParallel" %in% installed.packages()[,1]))
  install.packages("doParallel")
library("doParallel")

