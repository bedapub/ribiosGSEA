#' @include AllClasses.R AllGenerics.R
NULL

setMethod("gsName", "broadGseaResItem", function(object) return(object@geneset))
setMethod("gsName", "annoBroadGseaRes", function(object) sapply(object, gsName))
setMethod("gsName", "FisherResult", function(object) object@gsName)
setMethod("gsName", "FisherResultList", function(object) names(object@.Data))
setMethod("gsName", "GmtList", function(object) BioQC::gsName(object))

setMethod("gseaES", "broadGseaResItem", function(object) return(object@es))
setMethod("gseaES", "annoBroadGseaRes", function(object) {
  res <- sapply(object, gseaES)
  names(res) <- gsName(object)
  return(res)
})
setMethod("gseaES", "annoBroadGseaResList", function(object) {
  es <- lapply(object, gseaES)
  res <- vec2mat(es, sort.by="mean", decreasing=FALSE)
  return(res)
})

setMethod("gseaNES", "broadGseaResItem", function(object) return(object@nes))
setMethod("gseaNES", "annoBroadGseaRes", function(object) {
  res <- sapply(object, gseaNES)
  names(res) <- gsName(object)
  return(res)
})
setMethod("gseaNES", "annoBroadGseaResList", function(object) {
  nes <- lapply(object, gseaNES)
  res <- vec2mat(nes, sort.by="mean", decreasing=FALSE)
  return(res)
})

setMethod("gseaNP", "broadGseaResItem", function(object) return(object@np))
setMethod("gseaNP", "annoBroadGseaRes", function(object) {
  res <- sapply(object, gseaNP)
  names(res) <- gsName(object)
  return(res)
})
setMethod("gseaNP", "annoBroadGseaResList", function(object) {
  nps <- lapply(object, gseaNP)
  res <- vec2mat(nps, sort.by="mean", decreasing=FALSE)
  return(res)
})

setMethod("gseaFDR", "broadGseaResItem", function(object) return(object@fdr))
setMethod("gseaFDR", "annoBroadGseaRes", function(object) {
  res <- sapply(object, gseaFDR)
  names(res) <- gsName(object)
  return(res)
})
setMethod("gseaFDR", "annoBroadGseaResList", function(object) {
  fdrs <- lapply(object, gseaFDR)
  res <- vec2mat(fdrs, sort.by="mean", decreasing=FALSE)
  return(res)
})

setMethod("gseaFWER", "broadGseaResItem", function(object) {return(object@fwer)})
setMethod("gseaFWER", "annoBroadGseaRes", function(object) {
  res <- sapply(object, gseaFWER)
  names(res) <- gsName(object)
  return(res)
})
setMethod("gseaFWER", "annoBroadGseaResList", function(object) {
  fwers <- lapply(object, gseaFWER)
  res <- vec2mat(fwers, sort.by="mean", decreasing=FALSE)
  return(res)
})


setMethod("gsGeneIndices", "broadGseaResItem", function(object) return(object@geneIndices))
setMethod("gseaESprofile", "broadGseaResItem", function(object) return(object@esProfile))

setMethod("gseaCoreEnrichThr", "broadGseaResItem", function(object) return(object@coreEnrichThr))
setMethod("gseaCoreEnrichThr", "annoBroadGseaRes", function(object) {
  res <- sapply(object, gseaCoreEnrichThr)
  names(res) <- gsName(object)
  return(res)
})

setMethod("gsGenes", "annoBroadGseaResItem", function(object) return(object@gsGenes))
setMethod("gsGenes", "annoBroadGseaRes", function(object) {
  res <- lapply(object@.Data, gsGenes)
  names(res) <- gsName(object)
  return(res)
})
setMethod("gsGenes", "GmtList", function(object) return(BioQC::gsGenes(object)))


#' @importFrom ribiosUtils ulen
setMethod("gsSize", "GmtList", function(object) BioQC::gsSize(object))

#' @export
setMethod("gsNamespace", "GmtList", function(object) BioQC::gsNamespace(object))

#' @export
setMethod("gsDesc", "GmtList", function(object) BioQC::gsDesc(object))

setMethod("gsGeneValues", "annoBroadGseaResItem", function(object) return(object@gsGeneValues))
setMethod("gsGeneValues", "annoBroadGseaRes", function(object) {
  res <- lapply(object, gsGeneValues)
  names(res) <- gsName(object)
  return(res)
})

setMethod("gseaCoreEnrichGenes", "annoBroadGseaResItem", function(object) {
  gsGenes(object)[isGseaCoreEnrich(object)]
})
setMethod("gseaCoreEnrichGenes", "annoBroadGseaRes", function(object) {
  res <- lapply(object, gseaCoreEnrichGenes)
  names(res) <- gsName(object)
  return(res)
})
gseaLeadingEdgeGenes <- gseaCoreEnrichGenes

setMethod("gsGenes<-", c("annoBroadGseaResItem", "character"), function(object,value) {
  object@gsGenes <- value
  return(object)
})
setMethod("gsGeneValues<-", c("annoBroadGseaResItem", "numeric"), function(object, value) {
  object@gsGeneValues <- value
  return(object)
})
setMethod("isGseaCoreEnrich", "annoBroadGseaResItem", function(object) {
  nes <- gseaNES(object)
  thr <- gseaCoreEnrichThr(object)
  value <- gsGeneValues(object)
  if(nes<0) {
    value <= thr
  } else {
    value >= thr
  }
})

##setAs(from="list", to="gseaRes", def=function(from,to) {
##  haltifnot(all(sapply(from, function(x) is(x, "broadGseaResItem"))),
##            msg="Input list must be of broadGseaResItem objects")
##  res <- new("gseaRes", from)
##  return(res)
##})
setAs(from="list", to="annoBroadGseaRes", def=function(from,to) {
  haltifnot(all(sapply(from, function(x) is(x, "annoBroadGseaResItem"))),
            msg="Input list must be of annoBroadGseaResItem objects")
  res <- new("annoBroadGseaRes", from)
  return(res)
})
setAs(from="list", to="annoBroadGseaResList", def=function(from,to) {
  haltifnot(all(sapply(from, function(x) is(x, "annoBroadGseaRes"))),
            msg="Input list must be of annoBroadGseaRes objects")
  res <- new("annoBroadGseaResList", from)

  return(res)
})

##setMethod("gseaRes", "list", function(object) {
##  as(object, "gseaRes")
##})

setMethod("annoBroadGseaRes", "list", function(object) {
  return(as(object, "annoBroadGseaRes"))
})
##setMethod("[", "gseaRes", function(x, i) {
##  res <- callGeneric(x@.Data, i)
##  return(as(res, "gseaRes"))
##})
setMethod("[", "annoBroadGseaRes", function(x, i,...) {
  if(all(is.character(i)))
    i <- match(i, gsName(x))
  res <- callGeneric(x@.Data, i)
  return(as(res, "annoBroadGseaRes"))
})

setMethod("show", "broadGseaResItem", function(object) {
  gInd <- gsGeneIndices(object)
  fmt <- "GeneSet \"%s\" [%d genes]\nES=%1.3f; NES=%1.3f; \
Nominal P-value(NP)=%1.3f; FDR=%1.3f; FWER=%1.3f\
Indices:%s\nEnrichment Score (ES) profile:%s\nCore enrichment threshold of input value:%1.3f\n"
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
setMethod("show", "annoBroadGseaResItem", function(object) {
  gInd <- gsGeneIndices(object)
  fmt <- "AnnotatedGeneSet \"%s\" [%d genes]\nES=%1.3f; NES=%1.3f; \
Nominal P-value(NP)=%1.3f; FDR=%1.3f; FWER=%1.3f\
Enrichment Score (ES) profile:%s\nCore enrichment threshold of input value:%1.4f\
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

setMethod("annoBroadGseaResItem", "broadGseaResItem", function(object, genes, genevalues) {
  res <- as(object, "annoBroadGseaResItem")
  if(!(length(genes)==length(genevalues) && length(genes)==length(gsGeneIndices(object))))
    stop(sprintf("genes (#=%d) and genevalues (#=%d) must be of the same length as the gsGeneIndices (#=%d)",
                 length(genes), length(genevalues), length(gsGeneIndices(object))))
  gsGenes(res) <- genes
  gsGeneValues(res) <- genevalues
  return(res)
})

setMethod("show", "annoBroadGseaRes", function(object) {
  str <- sprintf("Annotated GSEA Results with %d gene sets\n",
                 length(object))
  cat(str)
})

## extending functions


#' Extract scores from GSEA results
#' 
#' One way to score GSEA results is to multiple the absolute value of log10
#' transformed p-values (nominal p-value, FDR, or FWER) with the sign of the
#' enrichment scores. This score is intuitive since it combines statistical
#' significance and the sign of regulation.
#' 
#' \code{gseaScores} takes care of the situation where some gene sets are
#' missing in one or more conditions.
#' 
#' @aliases gseaScore gseaScores
#' @param x An \code{annoBroadGseaRes} object
#' @param type Character string, the type of p-value used to calculate the
#' score.
#' @param ... Objects of \code{annoBroadGseaRes} to be compared
#' @param names Character strings, names given to the result score sets. See
#' examples below.
#' @return \code{gseaScore} returns a double vector of scores with gene set
#' names.
#' 
#' \code{gseaScores} returns a data frame of scores, with gene set names as row
#' names.
#' @author Jitao David Zhang <jitao_david.zhang@@roche.com>
#' @seealso \code{\link{gseaNP}}, \code{\link{gseaFDR}}, \code{\link{gseaFWER}}
#' to get p-values.
#' 
#' \code{\link{parseGSEAdir}}.
#' @export gseaScore
gseaScore <- function(x, type=c("fdr", "p", "fwer")) {
  type <- match.arg(type)
  if(type=="fdr") {
    val <- gseaFDR(x)
  } else if (type=="p") {
    val <- gseaNP(x)
  } else if (type=="fwer") {
    val <- gseaFWER(x)
  }
  val[val==0] <- min(val[val!=0], na.rm=TRUE)
  res <- -log10(val) * sign(gseaES(x))
  return(res)
}

gseaScores <- function(..., names=NULL, type=c("fdr", "p", "fwer")) {
  ll <- list(...)
  scores <- lapply(ll, gseaScore, type=type)
  setnames <- munion(lapply(scores, names))
  res <- as.data.frame(sapply(scores, function(x) x[match(setnames, names(x))]))
  rownames(res) <- setnames
  if(!is.null(names)) {
    haltifnot(length(names)==length(ll), msg="names length must match the input list")
    colnames(res) <- names
  }
  return(res)
}

## back compatibility
gsNames <- gsName
gsDescs <- gsDesc


##----------------------------------------##
## Fisher's exact test
##----------------------------------------##
setMethod("gsEffSize", "FisherResult", function(object) return(object@gsEffSize))
setMethod("gsEffSize", "FisherResultList", function(object) {
    return(sapply(object@.Data, gsEffSize))
})
setMethod("hits", "FisherResult", function(object) return(object@hits))

setMethod("gsNamespace", "FisherResult", function(object) return(object@gsNamespace))
setMethod("gsNamespace", "FisherResultList", function(object) sapply(object@.Data, gsNamespace))


setMethod("as.data.frame", "FisherResultList", function(x, row.names) {
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
## Fisher's exact test
##----------------------------------------##
setMethod("gsName", "FisherResultList", function(object,...) {
              names(object)
          })

setMethod("[", c("FisherResultList", "ANY", "missing", "missing"), function(x, i) {
              resList <- x@.Data[i]
              res <- new("FisherResultList", resList, input=x@input, universe=x@universe)
              return(res)
          })
setMethod("[", c("FisherResultList", "character", "character", "missing"),
          function(x, i,j, drop) {
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


setMethod("hits", "FisherResult", function(object) {
              object@hits
          })
setMethod("hits", "FisherResultList", function(object, geneset) {
              if(missing(geneset)) {
                  res <- lapply(object, function(x) x@hits)
              } else {
                  res <- genes(object[[geneset]])
              }
              return(res)
          })
setMethod("pValue", "FisherResult", function(object) {return(object@p)})
setMethod("pValue", "FisherResultList", function(object, ind, ...) {
              res <- sapply(object@.Data, pValue)
              if(!missing(ind)) {
                  return(res[ind])
              } else {
                  return(res)
              }
    
          })
setMethod("fdrValue", "FisherResult", function(object) {return(object@fdr)})
setMethod("fdrValue", "FisherResultList", function(object, ind, ...) {
              res <- sapply(object@.Data, fdrValue)
              if(!missing(ind)) {
                  return(res[ind])
              } else {
                  return(res)
              }
          })
setMethod("minPValue", "FisherResultList", function(object,...) {
              min(pValue(object))
          })
setMethod("minFdrValue", "FisherResultList", function(object,...) {
              min(fdrValue(object))
          })
setMethod("isSigGeneSet", c("FisherResultList", "numeric"),function(object,fdr) {
              fdrValue(object)<fdr
          })
setMethod("sigGeneSet", c("FisherResultList", "numeric"),function(object,fdr) {
              gsName(object)[isSigGeneSet(object, fdr)]
          })
setMethod("sigGeneSetTable", c("FisherResultList", "numeric"),function(object,fdr,...) {
              as.data.frame(obj[isSigGeneSet(object, fdr)])
          })
setMethod("topGeneSetTable", c("FisherResultList", "numeric"),function(object,N,...) {
              ps <- pValue(object)
              pOrd <- order(ps, decreasing=FALSE)[1:pmin(N, length(ps))]
              sub <- object[pOrd]
              return(as.data.frame(sub))
          })
setMethod("topOrSigGeneSetTable", c("FisherResultList", "numeric", "numeric"), function(object, N, fdr) {
              fdrV <- fdrValue(object)
              N <- pmax(N, sum(fdrV<fdr))
              return(topGeneSetTable(object, N, fdr))
              
          })
setMethod("topOrSigGeneSetTable", c("FisherResultList", "numeric", "missing"), function(object, N, fdr) {
              topOrSigGeneSetTable(object, N, 0.05)
          })
setMethod("topOrSigGeneSetTable", c("FisherResultList", "missing", "missing"), function(object, N, fdr) {
              topOrSigGeneSetTable(object, 10, 0.05)
          })

setMethod("print", "FisherResult", function(x, ...) {
              if(!is.na(gsNamespace(x)))
                  cat("Namespace:", gsNamespace(x), "\n")
              if(!is.na(gsName(x)))
                  cat("Name:", gsName(x), "\n")
              cat("Gene set size:", gsEffSize(x), "\n")
              cat(sprintf("Hits (%d):", length(hits(x))),
                  paste(hits(x), collapse=","), "\n")
              cat("Fisher's exact p value:", pValue(x), "\n")
              cat("BH FDR value:", fdrValue(x), "\n")
          })
setMethod("show", "FisherResult", function(object) {
              print(object)
          })
setMethod("print", "FisherResultList", function(x,...) {
              cat("--- One-sided Fisher's exact tests for gene sets ---\n")
              cat(sprintf("Total input genes: %d\n", length(x@input)))
              cat(sprintf("Gene universe: %d\n", length(x@universe)))
              cat(sprintf("Total gene sets: %d\n", length(x)))
              cat(sprintf("Minimal P-value: %e\n", minPValue(x)))
              cat(sprintf("Minimal FDR-value: %e\n", minFdrValue(x)))
          })

setMethod("show", "FisherResultList", function(object) {
              print(object)
          })

estimateFdr <- function(object) {
      system.time(ps <- sapply(object, pValue))
      fdrs <- rep(NA, length(ps))
      categories <- gsNamespace(object)
      categories[is.na(categories)] <- "NA"
      categories <- factor(categories)
      for(i in 1:nlevels(categories)) {
	  isCurr <- as.integer(categories)==i
	  fdrs[isCurr] <- p.adjust(ps[isCurr], "fdr")
      }
      for(i in seq(along=object)) {
	  object@.Data[[i]]@fdr <- fdrs[[i]]
      }
      return(object)
}

##----------------------------------------##
## migrated from ribiosNGS
##----------------------------------------##

#' Extract contrast names from an EdgeGSE object
#' @param object An \code{EdgeGSE} object
#' @importMethodsFrom ribiosExpression contrastNames
#' @export
setMethod("contrastNames",
          "EdgeGSE",
          function(object) contrastNames(as(object, "EdgeObject")))
