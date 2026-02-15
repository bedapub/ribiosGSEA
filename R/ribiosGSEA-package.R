#' @keywords internal
#' @useDynLib ribiosGSEA, .registration=TRUE, .fixes=`_ribiosGSEA_`
"_PACKAGE"

#' @importFrom Rcpp evalCpp
#' @importFrom ribiosExpression DesignContrast
#' @importFrom ribiosNGS EdgeObject
#' @importFrom edgeR DGEList
#' @importFrom utils read.table write.table
#' @importFrom methods as callGeneric is new
#' @importFrom stats ave cor fisher.test kmeans median p.adjust
#'             pchisq pt qt sd var model.matrix
NULL

## following lines make it not necessary to document
## re-exported functions

#' @export
ribiosExpression::DesignContrast

#' @export
ribiosExpression::designMatrix

#' @export
ribiosNGS::EdgeObject

#' @export
edgeR::DGEList

#' @export
limma::camera
