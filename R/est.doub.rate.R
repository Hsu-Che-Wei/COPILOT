#' Estimated doublet rate according to 10X genomics
#' @param ncell number of high quality cells remain after quality filtering

est.doub.rate <- function (ncell) {
  round(ncell*7.354e-04 + 1.260e-01,2)
}
