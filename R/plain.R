#' PrettyNum
#' @param x A numeric number

plain <- function(x) {
  prettyNum(x,big.mark = ",", scientific = FALSE)
}
