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

