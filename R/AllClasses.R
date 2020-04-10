setClass("gseaResItem",
         representation=list(geneset="character",
           "es"="numeric",
           "nes"="numeric",
           "np"="numeric",
           "fdr"="numeric",
           "fwer"="numeric",
           "geneIndices"="integer",
           "esProfile"="numeric",
           "coreEnrichThr"="numeric"))

setClass("annoGseaResItem",
         representation=list("gsGenes"="character",
           "gsGeneValues"="numeric"),
         contains="gseaResItem")

setClass("annoGseaRes", contains="list")
setClass("annoGseaResList", contains="list") #3 a list of annoGseaRes objects

setClass("GeneSetResult",
         representation=list(
             gsNamespace="character",
             gsName="character",
             gsEffSize="integer",
             p="numeric",
             fdr="numeric"),
         contains="VIRTUAL")
            
## Fisher's exact test
setClass("FisherResult",
         representation=list(hits="character"),
         contains="GeneSetResult")

setClass("FisherResultList",
         representation=list(
             inputName="character",
             input="character",
             universe="character"),
         contain="list")


## Camera result
setClass("CameraResult",
         representation=list(
             correlation="numeric",
             hits="character",
             score="numeric"),
         contains="GeneSetResult")
setClass("CameraResultList",
         representation=list(inputName="character"),
         contain="list")
            
##----------------------------------------##
## migrated from ribiosNGS
##----------------------------------------##

#' An EdgeGSE object contains gene-sets, enrichment method, and 
#' enrichment tables besides EdgeResult
#' @slot geneSets A GmtList
#' @slot method Gene-set enrichment method
#' @slot enrichTables A data.frame
#' @importClassesFrom ribiosNGS EdgeResult
#' @export
setClass("EdgeGSE",
         representation=list(geneSets="GmtList",
           method="character",
           enrichTables="data.frame"),
         contains="EdgeResult")

#' Build an EdgeGSE object
#' @param edgeObj An \code{EdgeObject}
#' @param gmtList A \code{GmtList} object
#' @returns An EdgeGSE object, with enrichTables as \code{NULL}
#' @importClassesFrom BioQC GmtList
#' @export
EdgeGSE <- function(edgeObj, gmtList) {
  haltifnot(all(c("GeneID", "GeneSymbol") %in% colnames(dgeList(edgeObj)$genes)),
            msg="Gene annotation of the edgeObj must contain columns 'GeneID' with EntrezGeneIDs and 'GeneSymbol' with official gene symbols")
  egse <- as(edgeObj,"EdgeGSE")
  egse@geneSets <- gmtList
  return(egse)
}

