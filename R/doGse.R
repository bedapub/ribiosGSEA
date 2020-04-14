#' @include gse.R
NULL

#' @describeIn doGse Gene-set enrichment with logFC gage
#' @export
gseWithLogFCgage <- function(edgeResult, gmtList) {
  gseRes <- logFCgage(edgeResult, gmtList)
  return(gseRes)
}

#' @describeIn doGse Gene-set enrichment with camera
#' @export
gseWithCamera <- function(edgeResult, gmtList) {
  gseRes <- camera.EdgeResult(edgeResult, gmtList)
  return(gseRes)
}

#' Perform gene-set enrichment (GSE) analysis
#' 
#' @param edgeResult An object of the class \code{EdgeObject}
#' @param gmtList An object of the class \code{GmtList}
#' 
#' The function performs gene-set enrichment analysis. By default,the CAMERA
#' method is applied. In case this is not successful, for instance because of
#' lack of biological replicates, the GAGE method (Generally Applicable
#' Gene-set Enrichment for pathway analysis) is applied.
#' @return An \code{EdgeGSE} object containing all information required to
#' reproduce the gene-set enrichment analysis results, as well as the
#' enrichment table. Apply \code{fullEnrichTable} to the object to extract a
#' \code{data.frame} containing results of the gene-set enrichment analysis.
#' @seealso \code{gseWithLogFCgage} and \code{gseWithCamera} are wrapped by
#' this function to perform analysis with GAGE and CAMERA, respectively.
#' \code{logFCgage} and \code{camera.EdgeResult} implements the logic, and
#' returns an object of the \code{EdgeGSE} class, which contains all relevant
#' information required to reproduce the analysis results.
#'
#' @examples
#' exMat <- matrix(rpois(120, 10), nrow=20, ncol=6)
#' exGroups <- gl(2,3, labels=c("Group1", "Group2"))
#' exDesign <- model.matrix(~0+exGroups)
#' exContrast <- matrix(c(-1,1), ncol=1, dimnames=list(c("Group1", "Group2"), c("Group2.vs.Group1")))
#' exDescon <- DesignContrast(exDesign, exContrast, groups=exGroups)
#' exFdata <- data.frame(GeneSymbol=sprintf("Gene%d", 1:nrow(exMat)))
#' exPdata <- data.frame(Name=sprintf("Sample%d", 1:ncol(exMat)),
#'                      Group=exGroups)
#' exObj <- EdgeObject(exMat, exDescon, 
#'                      fData=exFdata, pData=exPdata)
#' exDgeRes <- ribiosNGS::dgeWithEdgeR(exObj)
#' 
#' exGeneSets <- BioQC::GmtList(list(
#'     list(name="Set1", desc="set 1", genes=c("Gene1", "Gene2", "Gene3"), namespace="default"),
#'     list(name="Set2", desc="set 2", genes=c("Gene18", "Gene6", "Gene4"), namespace="default")
#' ))
#' exGse <- doGse(exDgeRes, exGeneSets)
#' fullEnrichTable(exGse)
#' 
#' exGseWithGage <- gseWithLogFCgage(exDgeRes, exGeneSets)
#' fullEnrichTable(exGseWithGage)
#' 
#' exGseWithCamera <- gseWithCamera(exDgeRes, exGeneSets)
#' fullEnrichTable(exGseWithCamera)
#' @importClassesFrom ribiosNGS EdgeResult
#' @importFrom ribiosNGS dgeWithEdgeR
#' @export doGse
doGse <- function(edgeResult, gmtList) {
  res <- try(gseWithCamera(edgeResult, gmtList))
  if(class(res)=="try-error") {
    res <- gseWithLogFCgage(edgeResult, gmtList)
  }
  return(res)
}
