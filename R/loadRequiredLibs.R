# make sure required libraries are installed and loaded
#
# devtools::use_package("dplyr")
# devtools::use_package("gdata")
# devtools::use_package("ggplot2")
# devtools::use_package("multtest")
devtools::use_package("ggrepel")
devtools::use_package("gtools")
devtools::use_package("readr")


#  #' @import multtest 
#' @import ggplot2 ggrepel readr
#' @importFrom gdata startsWith
#' @importFrom gtools mixedorder mixedsort 
#' @importFrom dplyr mutate filter
#' @importFrom biomaRt useMart
#' @importFrom biomaRt getBM
#' @importFrom matrixStats rowSums2 rowMins
#' @importFrom graphics plot
#'
NULL
# line containing NULL is needed for roxygen2
