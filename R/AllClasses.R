#' @include ribiosGSEA-package.R
NULL

##----------------------------------------------## 
## Classes for generic gene-set analysis results
##---------------------------------------------## 

#' A generic, virtual S4 class for gene-set analysis result
#' @export
setClass("GeneSetResult",
         representation=list(
             gsNamespace="character",
             gsName="character",
             gsEffSize="integer",
             p="numeric",
             fdr="numeric"),
         contains="VIRTUAL")

#' Result of Fisher's exact test
#' @export
setClass("FisherResult",
         representation=list(hits="character"),
         contains="GeneSetResult")

#' Print a FisherResult object
#' @param x A FisherResult object
#' @param ... Not used
#' @export
print.FisherResult <- function(x, ...) {
  if(!is.na(gsNamespace(x)))
      cat("Namespace:", gsNamespace(x), "\n")
  if(!is.na(gsName(x)))
      cat("Name:", gsName(x), "\n")
  cat("Gene set size:", gsEffSize(x), "\n")
  cat(sprintf("Hits (%d):", length(hits(x))),
      paste(hits(x), collapse=","), "\n")
  cat("Fisher's exact p value:", pValue(x), "\n")
  cat("BH FDR value:", fdrValue(x), "\n")
}

#' A list of results of Fisher's exact test
#' @export
setClass("FisherResultList",
         representation=list(
             inputName="character",
             input="character",
             universe="character"),
         contain="list")


#' Print a FisherResultList object
#'
#' @param x A \code{FisherResultList} object
#' @param ... Not used
#' @export
print.FisherResultList <- function(x, ...) {
   cat("--- One-sided Fisher's exact tests for gene sets ---\n")
   cat(sprintf("Total input genes: %d\n", length(x@input)))
   cat(sprintf("Gene universe: %d\n", length(x@universe)))
   cat(sprintf("Total gene sets: %d\n", length(x)))
   cat(sprintf("Minimal P-value: %e\n", minPValue(x)))
   cat(sprintf("Minimal FDR-value: %e\n", minFdrValue(x)))
}

## ## not used so far, to be deleted
## ## TODO: delete
## ## Camera result
## setClass("CameraResult",
##          representation=list(
##              correlation="numeric",
##              hits="character",
##              score="numeric"),
##          contains="GeneSetResult")
## setClass("CameraResultList",
##          representation=list(inputName="character"),
##          contain="list")

##---------------------------------------## 
## Classes for BROAD GSEA tools
##---------------------------------------## 

#' A S4 class representing the atom structure of results of the BROAD GSEA tool
#' @slot geneset Character, gene-set name
#' @slot es Numeric, enrichment score
#' @slot nes Numeric, normalised enrichment score
#' @slot np Numeric
#' @slot fdr Numeric, false discovery rate
#' @slot fwer Numeric, family-wise error rate
#' @slot geneIndices Integer vector, gene indices
#' @slot esProfile Numeric, enrichment score profile
#' @slot coreEnrichThr Numeric
#'
#' @export
setClass("broadGseaResItem",
         representation=list(geneset="character",
           "es"="numeric",
           "nes"="numeric",
           "np"="numeric",
           "fdr"="numeric",
           "fwer"="numeric",
           "geneIndices"="integer",
           "esProfile"="numeric",
           "coreEnrichThr"="numeric"))

#' Annotated BROAD GSEA result item
#'
#' @slot gsGenes Vector of character strings, gene-set genes
#' @slot gsGeneValues Vector of numeric values, statistics of gene-set genes
#'
#' @export
setClass("annoBroadGseaResItem",
         representation=list("gsGenes"="character",
           "gsGeneValues"="numeric"),
         contains="broadGseaResItem")

#' Annotated BROAD GSEA Results for one contrast
#' @export
setClass("annoBroadGseaRes", contains="list")

#' A list of annoBroadGseaRes objects
#' @export
setClass("annoBroadGseaResList", contains="list") #3 a list of annoBroadGseaRes objects
            
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
#' @return An EdgeGSE object, with enrichTables as \code{NULL}
#' @importFrom ribiosNGS dgeList
#' @importClassesFrom BioQC GmtList
#' @export
EdgeGSE <- function(edgeObj, gmtList) {
  haltifnot(all(c("GeneID", "GeneSymbol") %in% colnames(dgeList(edgeObj)$genes)),
            msg="Gene annotation of the edgeObj must contain columns 'GeneID' with EntrezGeneIDs and 'GeneSymbol' with official gene symbols")
  egse <- as(edgeObj,"EdgeGSE")
  egse@geneSets <- gmtList
  return(egse)
}

