#' Order strings by numbers in them
#' @param str A vector of character trings
#' @param ... Passed to \code{order}, by default \code{decreasing} is
#' \code{TRUE}, i.e. the descending order is reported.
#' @importFrom readr parse_number
#' @examples
#' orderByNumberInStr(c("D1", "D10", "D15", "D3.5"))
#' @seealso \code{\link{factorByNumberInStr}}, which makes factors with levels
#' ordered by numbers in the string
#' @export
orderByNumberInStr <- function(str, ...) {
  num <- readr::parse_number(as.character(str))
  res <- order(num, ...)
  return(res)
}

#' Make a factor vector from a character vector by the order of the parsed numbers
#' @param str Strings
#' @param decreasing Logical, whether decreasing or increasing order is desied, passed to \code{order}.
#' @examples
#' factorByNumberInStr(c("D1", "D10", "D15", "D3.5"))
#' factorByNumberInStr(c("D1", "D10", "D15", "D3.5"), decreasing=FALSE)
#' @seealso \code{\link{orderByNumberInStr}}, which returns the order of strings
#' by numbers in them
#' @export
factorByNumberInStr <- function(str, decreasing=TRUE) {
  ord <- orderByNumberInStr(str, decreasing=decreasing)
  res <- factor(str, levels=unique(str[ord]))
  return(res)
}
