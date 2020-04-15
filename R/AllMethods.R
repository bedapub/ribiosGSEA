#' @include AllClasses.R AllGenerics.R
NULL

##----------------------------------------##
## gsName
##----------------------------------------##

#' @describeIn gsName Get gene-set name from a broadGseaResItem object
#' @export
setMethod("gsName", "broadGseaResItem", function(object) return(object@geneset))

#' @describeIn gsName Get gene-set name from an annoBroadGseaRes object
#' @export
setMethod("gsName", "annoBroadGseaRes", function(object) sapply(object, gsName))

#' @describeIn gsName Get gene-set name from a FisherResult object
#' @export
setMethod("gsName", "FisherResult", function(object) object@gsName)

#' @describeIn gsName Get gene-set name from a FisherResultList object
#' @export
setMethod("gsName", "FisherResultList", function(object) names(object@.Data))

#' @describeIn gsName Get gene-set name from a GmtList object
#' @export
setMethod("gsName", "GmtList", function(object) BioQC::gsName(object))

#' @describeIn gsName Get gene-set name from a FisherResultList object
#' @export
setMethod("gsName", "FisherResultList", function(object,...) {
              names(object)
          })

##----------------------------------------##
## gseaES
##----------------------------------------##

#' @describeIn gseaES Get GSEA enrichment score from a broadGseaResItem object
#' @export
setMethod("gseaES", "broadGseaResItem", function(object) return(object@es))

#' @describeIn gseaES Get GSEA enrichment score from an annoBroadGseaRes object
#' @export
setMethod("gseaES", "annoBroadGseaRes", function(object) {
  res <- sapply(object, gseaES)
  names(res) <- gsName(object)
  return(res)
})

#' @describeIn gseaES Get GSEA enrichment score from an annoBroadGseaResList object
#' @export
setMethod("gseaES", "annoBroadGseaResList", function(object) {
  es <- lapply(object, gseaES)
  res <- vec2mat(es, sort.by="mean", decreasing=FALSE)
  return(res)
})

##----------------------------------------##
## gseaNES
##----------------------------------------##

#' @describeIn gseaNES Get GSEA normalized enrichment score 
#'     from a broadGseaResItem object
#' @export
setMethod("gseaNES", "broadGseaResItem", function(object) return(object@nes))

#' @describeIn gseaNES Get GSEA normalized enrichment score 
#'     from an annoBroadGseaRes object
#' @export
setMethod("gseaNES", "annoBroadGseaRes", function(object) {
  res <- sapply(object, gseaNES)
  names(res) <- gsName(object)
  return(res)
})

#' @describeIn gseaNES Get GSEA normalized enrichment score 
#'     from an annoBroadGseaResList object
#' @export
setMethod("gseaNES", "annoBroadGseaResList", function(object) {
  nes <- lapply(object, gseaNES)
  res <- vec2mat(nes, sort.by="mean", decreasing=FALSE)
  return(res)
})

##----------------------------------------##
## gseaNP
##----------------------------------------##

#' @describeIn gseaNP Get GSEA number of permutations
#'     from a broadGseaResItem object
#' @export
setMethod("gseaNP", "broadGseaResItem", function(object) return(object@np))

#' @describeIn gseaNP Get GSEA number of permutations
#'     from an annoBroadGseaRes object
#' @export
setMethod("gseaNP", "annoBroadGseaRes", function(object) {
  res <- sapply(object, gseaNP)
  names(res) <- gsName(object)
  return(res)
})

#' @describeIn gseaNP Get GSEA number of permutations
#'     from an annoBroadGseaResList object
#' @export
setMethod("gseaNP", "annoBroadGseaResList", function(object) {
  nps <- lapply(object, gseaNP)
  res <- vec2mat(nps, sort.by="mean", decreasing=FALSE)
  return(res)
})

##----------------------------------------##
## gseaES
##----------------------------------------##

#' @describeIn gseaFDR Get GSEA FDR values
#'     from a broadGseaResItem object
#' @export
setMethod("gseaFDR", "broadGseaResItem", function(object) return(object@fdr))

#' @describeIn gseaFDR Get GSEA FDR values
#'     from an annoBroadGseaRes object
#' @export
setMethod("gseaFDR", "annoBroadGseaRes", function(object) {
  res <- sapply(object, gseaFDR)
  names(res) <- gsName(object)
  return(res)
})

#' @describeIn gseaFDR Get GSEA FDR values
#'     from an annoBroadGseaResList object
#' @export
setMethod("gseaFDR", "annoBroadGseaResList", function(object) {
  fdrs <- lapply(object, gseaFDR)
  res <- vec2mat(fdrs, sort.by="mean", decreasing=FALSE)
  return(res)
})

##----------------------------------------##
## gseaFWER
##----------------------------------------##

#' @describeIn gseaFWER Get GSEA FWER values
#'     from a broadGseaResItem object
#' @export
setMethod("gseaFWER", "broadGseaResItem", function(object) {return(object@fwer)})

#' @describeIn gseaFWER Get GSEA FWER values
#'     from an annoBroadGseaRes object
#' @export
setMethod("gseaFWER", "annoBroadGseaRes", function(object) {
  res <- sapply(object, gseaFWER)
  names(res) <- gsName(object)
  return(res)
})

#' @describeIn gseaFWER Get GSEA FWER values
#'     from an annoBroadGseaResList object
#' @export
setMethod("gseaFWER", "annoBroadGseaResList", function(object) {
  fwers <- lapply(object, gseaFWER)
  res <- vec2mat(fwers, sort.by="mean", decreasing=FALSE)
  return(res)
})


##----------------------------------------##
## gsGeneIndices
##----------------------------------------##

#' @describeIn gsGeneIndices Get gene-set gene indices
#'     from a broadGseaResItem object, returning a vector of integers.
#' @export
setMethod("gsGeneIndices", "broadGseaResItem", function(object) return(object@geneIndices))

##----------------------------------------##
## gseaESprofile
##----------------------------------------##

#' @describeIn gseaESprofile Get GSEA enrichment profile
#'     from a broadGseaResItem object
#' @export
setMethod("gseaESprofile", "broadGseaResItem", function(object) return(object@esProfile))

##----------------------------------------##
## gseaCoreEnrichThr
##----------------------------------------##

#' @describeIn gseaCoreEnrichThr Get the threshold value of GSEA core enrichment
#'     from a broadGseaResItem object
#' @export
setMethod("gseaCoreEnrichThr", "broadGseaResItem", function(object) return(object@coreEnrichThr))

#' @describeIn gseaCoreEnrichThr Get the threshold value of GSEA core enrichment
#'     from an annoBroadGseaRes object
#' @export
setMethod("gseaCoreEnrichThr", "annoBroadGseaRes", function(object) {
  res <- sapply(object, gseaCoreEnrichThr)
  names(res) <- gsName(object)
  return(res)
})

##----------------------------------------##
## gsGenes
##----------------------------------------##

#' @describeIn gsGenes Get gene-set genes
#'     from a broadGseaResItem object, returning a character string vector.
#' @export
setMethod("gsGenes", "annoBroadGseaResItem", function(object) return(object@gsGenes))

#' @describeIn gsGenes Get gene-set genes from an annoBroadGseaRes object,
#' returning a list of character string vectors.
#' @export
setMethod("gsGenes", "annoBroadGseaRes", function(object) {
  res <- lapply(object@.Data, gsGenes)
  names(res) <- gsName(object)
  return(res)
})

#' @describeIn gsGenes Get gene-set genes from a GmtList object, returning a
#' list of character string vector. It uses the implementation in BioQC.
# '@export
setMethod("gsGenes", "GmtList", function(object) return(BioQC::gsGenes(object)))

#' @describeIn gsGenes-set Assign gene-set genes to annoBroadGseaResItem
#' @export
setMethod("gsGenes<-", c("annoBroadGseaResItem", "character"), function(object,value) {
  object@gsGenes <- value
  return(object)
})

##----------------------------------------##
## gsNamespace
##----------------------------------------##

#' @describeIn gsNamespace Return gene-set namespace from a GmtList object
#' @export
setMethod("gsNamespace", "GmtList", function(object) BioQC::gsNamespace(object))

#' @describeIn gsNamespace Return gene-set namespace from a FisherResult object
#' @export
setMethod("gsNamespace", "FisherResult", function(object) return(object@gsNamespace))

#' @describeIn gsNamespace Return gene-set namespace from a FisherResultList 
#' object.
#' @export
setMethod("gsNamespace", "FisherResultList", function(object) sapply(object@.Data, gsNamespace))


##----------------------------------------##
## gsGeneValues
##----------------------------------------##

#' @describeIn gsGeneValues Return values associated with the genes in a
#' gene-set in an annoBroadGseaResItem object in a numeric vector.
#' @export
setMethod("gsGeneValues", "annoBroadGseaResItem", function(object) return(object@gsGeneValues))

#' @describeIn gsGeneValues Return values associated with the genes in a
#' gene-set in an annoBroadGseaRes object in a list of numeric vectors.
#' @export
setMethod("gsGeneValues", "annoBroadGseaRes", function(object) {
  res <- lapply(object, gsGeneValues)
  names(res) <- gsName(object)
  return(res)
})

#' @describeIn gsGeneValues-set Assign values associated with gene-set genes to
#' an annoBraoadGseaResItem object
#' @export
setMethod("gsGeneValues<-", c("annoBroadGseaResItem", "numeric"),
	  function(object, value) { 
		  object@gsGeneValues <- value 
		  return(object)
	  })

##----------------------------------------##
## gseaCoreEnrichGenes
##----------------------------------------##

#' Return a vector of logical values, indicating whether genes belong to core
#' enrichment or not
#' @param object An \code{annoBroadGseaResItem} object
#' @return A logical vector
isGseaCoreEnrich <- function(object) {
  nes <- gseaNES(object)
  thr <- gseaCoreEnrichThr(object)
  value <- gsGeneValues(object)
  if(nes<0) {
    value <= thr
  } else {
    value >= thr
  }
}

#' @describeIn gseaCoreEnrichGenes Return core enriched genes (also known as
#' leading-edge genes) in an annoBroadGseaResItem object as a character string
#' vector.
#' @export
setMethod("gseaCoreEnrichGenes", "annoBroadGseaResItem", function(object) {
  gsGenes(object)[isGseaCoreEnrich(object)]
})

#' @describeIn gseaCoreEnrichGenes Return core enriched genes (also known as
#' leading-edge genes) in an annoBroadGseaRes object as a list of character 
#' string vectors.
#' @export
setMethod("gseaCoreEnrichGenes", "annoBroadGseaRes", function(object) {
  res <- lapply(object, gseaCoreEnrichGenes)
  names(res) <- gsName(object)
  return(res)
})

## #' @describeIn gseaCoreEnrichGenes gseaLeadingEdgeGenes is synonymous with
## #' gseaCoreEnrichGenes
## #' @export
## gseaLeadingEdgeGenes <- gseaCoreEnrichGenes

##----------------------------------------##
## setAs
##----------------------------------------##

#' Convert a list of annoBroadGseaResItem objects to an annoBroadGseaRes object
#' @param from A list of annoBroadGseaResItem objects
#' @param to An annoBroadGseaRes objecct
#' @export
setAs(from="list", to="annoBroadGseaRes", def=function(from,to) {
  haltifnot(all(sapply(from, function(x) is(x, "annoBroadGseaResItem"))),
            msg="Input list must be of annoBroadGseaResItem objects")
  res <- new("annoBroadGseaRes", from)
  return(res)
})

#' Convert a list of annoBroadGseaRes to an annoBroadGseaResList object
#' @param from A list of annoBroadGseaRes objects
#' @param to An annoBroadGseaResList objecct
#' @export
setAs(from="list", to="annoBroadGseaResList", def=function(from,to) {
  haltifnot(all(sapply(from, function(x) is(x, "annoBroadGseaRes"))),
            msg="Input list must be of annoBroadGseaRes objects")
  res <- new("annoBroadGseaResList", from)

  return(res)
})

setMethod("annoBroadGseaRes", "list", function(object) {
  return(as(object, "annoBroadGseaRes"))
})


##----------------------------------------##
## [
##----------------------------------------##

#' Subset an annoBroadGseaRes object
#'
#' @param x An annoBroadGseaRes object
#' @param i An integer or logical subsetting index
#' @param j Not used
#' @param ... Not used
#' @param drop Not used
#' @return A subset of the original data as an annoBroadGseaRes object
#'
#' @export
setMethod("[", "annoBroadGseaRes", function(x, i, j, ..., drop=FALSE) {
  if(all(is.character(i)))
    i <- match(i, gsName(x))
  res <- callGeneric(x@.Data, i)
  return(as(res, "annoBroadGseaRes"))
})


#' Subset a FisherResultList object by indexing
#'
#' @param x A FisherResultList object
#' @param i An integer or logical subsetting index
#' @param j Not used
#' @param ... Not used
#' @param drop Not used
#' @return A subset of the original data as an FisherResultList object
#'
#' @export
setMethod("[", c("FisherResultList", "ANY", "missing", "missing"), function(x, i, j, ..., drop=FALSE) {
              resList <- x@.Data[i]
              res <- new("FisherResultList", resList, input=x@input, universe=x@universe)
              return(res)
          })

#' Subset a FisherResultList object by namespace and name
#'
#' @param x A FisherResultList object
#' @param i Character string, gene-set namespace
#' @param j Character string, gene-set name
#' @param ... Not used
#' @param drop Not used
#'
#  @return If more than one elements are found, a \code{FisherResultList} is
#  returned. If only one element is found, a \code{FisherResult} object is
#  returned.
#' @export
setMethod("[", c("FisherResultList", "character", "character", "missing"),
          function(x, i,j, ..., drop) {
              isNamespace <- gsNamespace(x) %in% i
              isName <- names(x) %in% j
              isSel <- isNamespace & isName
              if(sum(isSel)==1) {
                  return(x[[which(isSel)]])
              } else if (sum(isSel)>1) {
                  return(x[isSel])
              } else {
                  stop(sprintf("No element found for namespace %s and gene set %s!\n",
                               i, j))
              }
          })

##----------------------------------------##
## show
##----------------------------------------##

#' Show a broadGseaResItem object
#' @param object A broadGseaResItem object 
#' export
setMethod("show", "broadGseaResItem", function(object) {
  gInd <- gsGeneIndices(object)
  fmt <- "GeneSet \"%s\" [%d genes]\nES=%1.3f; NES=%1.3f; \
Nominal P-value(NP)=%1.3f; FDR=%1.3f; FWER=%1.3f\
Indices:%s\nEnrichment Score (ES) profile:%s\n\
Core enrichment threshold of input value:%1.3f\n"
  str <- sprintf(fmt,
                 gsName(object),
                 length(gInd),
                 gseaES(object),
                 gseaNES(object),
                 gseaNP(object),
                 gseaFDR(object),
                 gseaFWER(object),
                 paste(gInd, collapse=","),
                 paste(gseaESprofile(object), collapse=","),
                 gseaCoreEnrichThr(object))
  cat(str)
})

#' Show an annoBroadGseaResItem object
#' @param object An annobroadGseaResItem object 
#' export
setMethod("show", "annoBroadGseaResItem", function(object) {
  gInd <- gsGeneIndices(object)
  fmt <- "AnnotatedGeneSet \"%s\" [%d genes]\nES=%1.3f; NES=%1.3f; \
Nominal P-value(NP)=%1.3f; FDR=%1.3f; FWER=%1.3f\
Enrichment Score (ES) profile:%s\n\
Core enrichment threshold of input value:%1.4f\
GeneNames:%s\nGene input values:%s\n"
  str <- sprintf(fmt,
                 gsName(object),
                 length(gInd),
                 gseaES(object),
                 gseaNES(object),
                 gseaNP(object),
                 gseaFDR(object),
                 gseaFWER(object),
                 paste(gseaESprofile(object), collapse=","),
                 gseaCoreEnrichThr(object),
                 paste(gsGenes(object),collapse=","),
                 paste(gsGeneValues(object), collapse=","))
  cat(str)
})

#' Show a anonBroadGseaRes object
#' @param object A annoBroadGseaRes object 
#' export
setMethod("show", "annoBroadGseaRes", function(object) {
  str <- sprintf("Annotated GSEA Results with %d gene sets\n",
                 length(object))
  cat(str)
})

setMethod("show", "FisherResult", function(object) {
              print(object)
          })

setMethod("show", "FisherResultList", function(object) {
              print(object)
          })

##----------------------------------------##
## gene-set effect size
##----------------------------------------##

#' @describeIn gsEffSize Gene-set effect size, returning a numeric value.
#' @return
setMethod("gsEffSize", "FisherResult", function(object) return(object@gsEffSize))

#' @describeIn gsEffSize Gene-set effect size, returning a numeric vector.
#' @return
setMethod("gsEffSize", "FisherResultList", function(object) {
    return(sapply(object@.Data, gsEffSize))
})


##----------------------------------------##
## as.data.frame
##----------------------------------------##

#' Convert an FisherResultList object into a data.frame
#' @param x An FisherResultList object
#' @param row.names Character strings.
#' @return A \code{data.frame}
#' @export
setMethod("as.data.frame", "FisherResultList", function(x, row.names=NULL) {
              categories <- sapply(x, gsNamespace)
              genesets <- sapply(x, gsName) ## TODO: gsName
              ps <- sapply(x, pValue)
              fdrs <- sapply(x, fdrValue)
              hits <- lapply(x, hits)
              hitCounts <- sapply(hits, length)
              gsEffSize <- sapply(x, gsEffSize)
              inputSize <- length(x@input)
              universeSize <- length(x@universe)
              hitPrint <- sapply(hits, paste, collapse=",")
              data.frame(Namespace=categories,
                         GeneSet=genesets,
                         Pvalue=ps,
                         FDRvalue=fdrs,
                         HitCount=hitCounts,
                         InputSize=inputSize,
                         GeneSetEffectiveSize=gsEffSize,
                         UniverseSize=universeSize,
                         Hits=hitPrint,
                         row.names=row.names)
          })


##----------------------------------------##
## hits
##----------------------------------------##

#' @describeIn hits Return hits from a FisherResult object
#' @export
setMethod("hits", "FisherResult", function(object) {
              object@hits
          })

#' @describeIn hits Return hits from a FisherResultList object, returning a list
#' if \code{geneset} is missing, or gene-set genes if \code{geneset} is present.
#' @param geneset Character string, gene-set name
#' @export
setMethod("hits", "FisherResultList", function(object, geneset) {
              if(missing(geneset)) {
                  res <- lapply(object, function(x) x@hits)
              } else {
                  res <- gsGenes(object[[geneset]])
              }
              return(res)
          })

##----------------------------------------##
## pValue
##----------------------------------------##

#' @describeIn pValue Return the p-value from a FisherResult
#' @export
setMethod("pValue", "FisherResult", function(object) {return(object@p)})

#' @describeIn pValue Return the p-values from a FisherResultList. If \code{ind}
#' is missing, all p-values are returned; otherwise, the subset indicated by
#' \code{ind} is returned.
#' @param ind An integer or logical vector for subsetting
#' @export
setMethod("pValue", "FisherResultList", function(object, ind, ...) {
              res <- sapply(object@.Data, pValue)
              if(!missing(ind)) {
                  return(res[ind])
              } else {
                  return(res)
              }
    
          })

#' @describeIn pValue Return the FDR-value from a FisherResult
#' @export
setMethod("fdrValue", "FisherResult", function(object) {return(object@fdr)})

#' @describeIn pValue Return the FDR-values from a FisherResultList. If \code{ind}
#' is missing, all FDR-values are returned; otherwise, the subset indicated by
#' \code{ind} is returned.
#' @export
setMethod("fdrValue", "FisherResultList", function(object, ind, ...) {
              res <- sapply(object@.Data, fdrValue)
              if(!missing(ind)) {
                  return(res[ind])
              } else {
                  return(res)
              }
          })

##----------------------------------------##
## Exported from ribiosExpression
##----------------------------------------##

#' Extract contrast names from an EdgeGSE object
#' @param object An \code{EdgeGSE} object
#' @importMethodsFrom ribiosExpression contrastNames
#' @export
setMethod("contrastNames",
          "EdgeGSE",
          function(object) contrastNames(as(object, "EdgeObject")))

##----------------------------------------##
## Fisher's exact test
##----------------------------------------##

#' Perform Fisher's exact test on a gene set
#'
#' @param genes a collection of genes of which over-representation of the gene set is tested
#' @param genesets A vector of character strings, genes belonging to one gene
#' set.
#' @param universe universe of genes
#' @param gsName gene set name, can be left missing
#' @param gsNamespace gene set namespace name, can be left missing
#' @param makeUniqueNonNA Logical, whether genes, geneSetGenes, and universe should be filtered to remove NA and made unique. The default is set to \code{TRUE}. When the uniqueness and absence of NA is ensured, this flag can be set to \code{FALSE} to accelerate the operation.
#' @param checkUniverse Logical, if \code{TRUE}, then genes that are in \code{genes} but are not in \code{universe} are appended to \code{universe}
#' @param useEASE Logical, whether to use the EASE method to report the p-value. 
#'
#' This function performs one-sided Fisher's exact test to test the over-representation of gene set genes in the input gene list.
#' 
#' If \code{useEASE} is \code{TRUE}, one gene is penalized (removed) within \code{geneSetGenes} that are in \code{genes} and calculating the resulting Fisher exact probability for that namespace. The theoretical basis of the EASE score lies in the concept of jackknifing a probability. See Hosack \emph{et al.} for details.
#'
#' @references 
#' \describe{
#'   \item{Hosack \emph{et al.}}{Hosack, Douglas A., Glynn Dennis, Brad T. Sherman, H. Clifford Lane, and Richard A. Lempicki. Identifying Biological Themes within Lists of Genes with EASE. Genome Biology 4 (2003): R70. \url{https://doi.org/10.1186/gb-2003-4-10-r70}}
#' }
#' 
#' @note Duplicated items in genes, genesets' genes, and the universe are per default removed
#' 
#' @examples
#' myGenes <- LETTERS[1:3]
#' myGeneSet1 <- LETTERS[1:6]
#' myGeneSet2 <- LETTERS[4:7]
#' myUniverse <- LETTERS
#' fisherTest(genes=myGenes, genesets=myGeneSet1, universe=myUniverse)
#' fisherTest(genes=myGenes, genesets=myGeneSet2, universe=myUniverse)
#' fisherTest(genes=myGenes, genesets=myGeneSet1, universe=myUniverse, 
#'            gsName="My gene set1", gsNamespace="Letters")
#'
#' ## note that duplicated items are removed by default
#' resWoRp <- fisherTest(genes=rep(myGenes,2), genesets=myGeneSet1, 
#'                       universe=myUniverse)
#' resWithRp <- fisherTest(genes=rep(myGenes,2), genesets=myGeneSet1, 
#'                       universe=rep(myUniverse,2))
#' identical(resWoRp, resWithRp)
#' 
#' resWithRpNoUnique <- fisherTest(genes=rep(myGenes,2), genesets=myGeneSet1, 
#'            universe=rep(myUniverse,2), makeUniqueNonNA=FALSE)
#' identical(resWoRp, resWithRpNoUnique)
#' @export
setMethod("fisherTest", c("character", "character", "character"),
          function(genes, genesets, universe, gsName, gsNamespace,
                   makeUniqueNonNA=TRUE, 
                   checkUniverse=TRUE,
                   useEASE=FALSE) {
              if(missing(gsName))
                  gsName <- as.character(NA)
              if(missing(gsNamespace))
                  gsNamespace <- as.character(NA)
              coreRes <- gsFisherTestCore(genes = genes, 
                           geneSetGenes = genesets,
                           universe = universe,
                           makeUniqueNonNA = makeUniqueNonNA,
                           checkUniverse = checkUniverse,
                           useEASE = useEASE)
              new("FisherResult",
                  gsNamespace=gsNamespace,
                  gsName=gsName,
                  gsEffSize=coreRes$gsEffSize,
                  hits=coreRes$hits,
                  p=coreRes$p,
                  fdr=coreRes$p)
          })

#' Perform Fisher's exact test on a GeneSet object
#'
#' @param genes a collection of genes of which over-representation of the gene set is tested
#' @param genesets A \code{GmtList} object.
#' @param universe universe of genes
#' @param makeUniqueNonNA Logical, whether genes and universe should be filtered to remove NA and made unique. The default is set to \code{TRUE}. When the uniqueness and absence of NA is ensured, this flag can be set to \code{FALSE} to accelerate the operation.
#' @param checkUniverse Logical, if \code{TRUE}, then genes that are in \code{genes} but are not in \code{universe} are appended to \code{universe}
#' @param useEASE Logical, whether to use the EASE method to report the p-value. 
#'
#' This function performs one-sided Fisher's exact test to test the over-representation of gene set genes in the input gene list.
#'
#' @importClassesFrom BioQC GmtList
#' @examples
#' myGenes <- LETTERS[1:3]
#' myS4GeneSet1 <- list(name="GeneSet1", desc="GeneSet", 
#'     genes=LETTERS[1:6], namespace="My namespace 1")
#' myS4GeneSet2 <- list(name="GeneSet1", desc="GeneSet", 
#'     genes=LETTERS[2:7], namespace="My namespace 2")
#' myUniverse <- LETTERS
#' fisherTest(myGenes, myS4GeneSet1, myUniverse)
#' fisherTest(myGenes, myS4GeneSet2, myUniverse)
#' @export
setMethod("fisherTest", c("character", "list", "character"),
          function(genes, genesets, universe,
                   makeUniqueNonNA=TRUE, 
                   checkUniverse=TRUE,
                   useEASE=FALSE) {
            if(makeUniqueNonNA) {
              genes <- uniqueNonNA(genes)
              universe <- uniqueNonNA(universe)
            }
            ## gsGenes(genesets) are garanteed to be unique and non-NA
            ## therefore fisherTest now takes makeUniqueNonNA=FALSE
            fisherTest(genes=genes, genesets=genesets$genes,
                       universe=universe,
                       gsName=genesets$name, 
                       gsNamespace=genesets$namespace,
                       makeUniqueNonNA=FALSE,
                       checkUniverse=checkUniverse,
                       useEASE=useEASE)
          })


##----------------------------------------##
## obsolete functions
##----------------------------------------##

## estimateFdr <- function(object) {
##       ps <- sapply(object, pValue)
##       fdrs <- rep(NA, length(ps))
##       categories <- gsNamespace(object)
##       categories[is.na(categories)] <- "NA"
##       categories <- factor(categories)
##       for(i in 1:nlevels(categories)) {
## 	  isCurr <- as.integer(categories)==i
## 	  fdrs[isCurr] <- p.adjust(ps[isCurr], "fdr")
##       }
##       for(i in seq(along=object)) {
## 	  object@.Data[[i]]@fdr <- fdrs[[i]]
##       }
##       return(object)
## }
