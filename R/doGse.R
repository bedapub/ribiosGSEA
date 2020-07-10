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
gseWithCamera <- function(edgeResult, gmtList, doParallel=TRUE) {
  gseRes <- camera.EdgeResult(edgeResult, gmtList, doParallel=doParallel)
  return(gseRes)
}

#' Perform gene-set enrichment (GSE) analysis
#' 
#' @param edgeResult An object of the class \code{EdgeObject}
#' @param gmtList An object of the class \code{GmtList}
#' @param doParallel Logical, whether \code{parallel::mclapply} should be used. Since at the current setting it makes a job running forever, use \code{TRUE} only if you are debugging the code.
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
#' 
#' \dontrun{
#'   exMat <- matrix(rpois(120000, 10), nrow=20000, ncol=12)
#'   exGroups <- gl(4,3, labels=c("Group1", "Group2", "Group3", "Group4"))
#'   exDesign <- model.matrix(~0+exGroups)
#'   exContrast <- matrix(c(-1,1,0,0, 0,0,-1,1),
#'      ncol=2, byrow=FALSE,
#'      dimnames=list(c("Group1", "Group2", "Group3", "Group4"), 
#'        c("Group2.vs.Group1", "Group4.vs.Group3")))
#'   exDescon <- DesignContrast(exDesign, exContrast, groups=exGroups)
#'   exFdata <- data.frame(GeneSymbol=sprintf("Gene%d", 1:nrow(exMat)))
#'   exPdata <- data.frame(Name=sprintf("Sample%d", 1:ncol(exMat)),
#'                        Group=exGroups)
#'   exObj <- EdgeObject(exMat, exDescon, 
#'                        fData=exFdata, pData=exPdata)
#'   exDgeRes <- ribiosNGS::dgeWithEdgeR(exObj)
#'   
#'   ngeneset <- 1000
#'   genesetSizes <- round(runif(ngeneset)*100)+1
#'   exGeneSets <- BioQC::GmtList(lapply(seq(1:ngeneset), function(i) {
#'     name <- paste0("GeneSet", i)
#'     desc <- paste0("GeneSet", i)
#'     genes <- sample(exFdata$GeneSymbol, genesetSizes[i])
#'     res <- list(name=name, desc=desc, genes=genes, namespace="default")
#'   }))
#'   exGse <- doGse(exDgeRes, exGeneSets)
#' }
#' @importClassesFrom ribiosNGS EdgeResult
#' @importFrom ribiosNGS dgeWithEdgeR
#' @export doGse
doGse <- function(edgeResult, gmtList, doParallel=FALSE) {
  res <- try(gseWithCamera(edgeResult, gmtList, doParallel=doParallel))
  if(class(res)=="try-error") {
    res <- gseWithLogFCgage(edgeResult, gmtList)
  }
  return(res)
}
