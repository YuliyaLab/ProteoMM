# make sure required libraries are installed and loaded
#
# devtools::use_package("dplyr")
devtools::use_package("gdata")
devtools::use_package("ggplot2")
devtools::use_package("ggrepel")
devtools::use_package("gtools")
devtools::use_package("multtest")
devtools::use_package("readr")


#' @import gdata ggplot2 ggrepel gtools multtest readr
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom biomaRt useMart
#' @importFrom biomaRt getBM
#'
NULL
# line containing NULL is needed for roxygen2
