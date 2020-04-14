#' The ribiosGSEA package
#'
#' The ribiosGSEA package supports gene-set analysis
#'
#' @docType package
#' @useDynLib ribiosGSEA, .registration=TRUE, .fixes=`_ribiosGSEA_`
#' @name ribiosGSEA-package
NULL

#' @importFrom ribiosExpression DesignContrast
#' @importFrom ribiosNGS EdgeObject
#' @importFrom edgeR DGEList
#' @importFrom stats model.matrix
#' @importFrom utils write.table
NULL

#' @export
ribiosExpression::DesignContrast

#' @export
ribiosNGS::EdgeObject

#' @export
edgeR::DGEList
