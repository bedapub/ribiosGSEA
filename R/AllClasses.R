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
             gsEffectiveSize="integer",
             p="numeric",
             fdr="numeric"),
         contains="VIRTUAL")

#' Result of Fisher's exact test
#' @export
setClass("FisherResult",
         representation=list(hits="character"),
         contains="GeneSetResult")


#' A list of results of Fisher's exact test
#' @export
setClass("FisherResultList",
         representation=list(
             inputName="character",
             input="character",
             universe="character"),
         contain="list")


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
setClass("BroadGseaResItem",
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
setClass("AnnoBroadGseaResItem",
         representation=list("gsGenes"="character",
           "gsGeneValues"="numeric"),
         contains="BroadGseaResItem")

#' Annotated BROAD GSEA Results for one contrast
#' @export
setClass("AnnoBroadGseaRes", contains="list")

#' A list of AnnoBroadGseaRes objects
#' @export
setClass("AnnoBroadGseaResList", contains="list") #3 a list of AnnoBroadGseaRes objects
            
