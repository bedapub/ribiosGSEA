#' @include AllGenerics.R AllMethods.R ribiosGSEA-package.R
NULL


#' Return a logical vector indicating whether a gene-set is significantly
#' enriched or not, given the FDR threshold
#' @param object A FisherResultList object
#' @param fdr Numeric, FDR value threshold
#' @return A logical vector
#' @export
isSigGeneSet <- function(object, fdr=0.05) { 
	return(fdrValue(object)<=fdr) 
}

#' Return names of gene-sets that are significantly enriched given the 
#' FDR threshold
#' @param object A FisherResultList object
#' @param fdr Numeric, FDR value threshold
#' @return A character vector
#' @export
sigGeneSet <- function(object,fdr) {
    gsName(object)[isSigGeneSet(object, fdr)]
}

#' Return a data.frame of significantly enriched gene-sets
#' @param object A FisherResultList object
#' @param fdr Numeric, FDR value threshold
#' @return A \code{data.frame}
#' @export
sigGeneSetTable <- function(object, fdr) {
  as.data.frame(object[isSigGeneSet(object, fdr)])
}

#' Return the minimal p-value from a FisherResultList
#' @param object A FisherResultList object
#' @return A numeric value
#' @export
minPvalue <- function(object) {
	min(pValue(object))
}

#' Return the minimal FDR value from a FisherResultList
#' @param object A FisherResultList object
#' @return A numeric value
#' @export
minFDRvalue <- function(object) {
              min(fdrValue(object))
}

#' Return a data.frame of top gene-sets with the lowest p-values
#' @param object An FisherResultList object
#' @param N Integer, the number of returned gene-sets
# '@return A data.frame
#' @export
topGeneSetTable <- function(object, N) {
      ps <- pValue(object)
      pOrd <- order(ps, decreasing=FALSE)[1:pmin(N, length(ps))]
      sub <- object[pOrd]
      return(as.data.frame(sub))
}

#' Return a data.frame of significantly enriched gene-sets with a minimum number
#' @param object An FisherResultList object
#' @param fdr Numeric, the treshold of FDR value
#' @param N Integer, the number of returned gene-sets
#' The total number of returned gene-sets are determined by the maximum of
#' \code{N} and the counts of gene-sets that have \code{FDR} lower than
#' \code{fdr}.
#' @return A \code{data.frame}.
#' @export
topOrSigGeneSetTable <- function(object, fdr=0.05, N=10) {
      fdrV <- fdrValue(object)
      N <- pmax(N, sum(fdrV<=fdr))
      return(topGeneSetTable(object, N))
}

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

#' The core algorithm to perform Fisher's exact test on a gene set
#'
#' @param genes Character vector, a collection of genes of which over-representation of the gene set is tested
#' @param geneSetGenes Character vector, genes belonging to a gene set
#' @param universe Character vector, universe of genes
#' @param makeUniqueNonNA Logical, whether genes, geneSetGenes, and universe should be filtered to remove NA and made unique. The default is set to \code{TRUE}. When the uniqueness and absence of NA is ensured, this flag can be set to \code{FALSE} to accelerate the operation.
#' @param checkUniverse Logical, if \code{TRUE}, then genes that are in \code{genes} but are not in \code{universe} are appended to \code{universe}
#' @param useEASE Logical, whether to use the EASE method to report the p-value. 
#'
#' This function performs one-sided Fisher's exact test to test the over-representation of the genes given as \code{geneSetGenes} in the input \code{genes} list.
#' 
#' If \code{useEASE} is \code{TRUE}, one gene is penalized (removed) within \code{geneSetGenes} that are in \code{genes} and calculating the resulting Fisher exact probability for that namespace. The theoretical basis of the EASE score lies in the concept of jackknifing a probability. See Hosack \emph{et al.} for details.
#'
#' @return A list of three elements
#' \enumerate{
#'   \item p The p-value of one-sided (over-representation of the Fisher's test)
#'   \item gsEffSize Gene-set's effective size, namely number of genes that are in the universe
#'   \item hits Character vector, genes that are found in the gene sets
#' }
#' @references 
#' \describe{
#'   \item{Hosack \emph{et al.}}{Hosack, Douglas A., Glynn Dennis, Brad T. Sherman, H. Clifford Lane, and Richard A. Lempicki. Identifying Biological Themes within Lists of Genes with EASE. Genome Biology 4 (2003): R70. \url{https://doi.org/10.1186/gb-2003-4-10-r70}}
#' }
#' @examples
#' myGenes <- LETTERS[1:3]
#' myGeneSet1 <- LETTERS[1:6]
#' myGeneSet2 <- LETTERS[4:7]
#' myUniverse <- LETTERS
#' gsFisherTestCore(myGenes, myGeneSet1, myUniverse)
#' gsFisherTestCore(myGenes, myGeneSet2, myUniverse)
#' 
#' ## use EASE for conservative estimating
#' gsFisherTestCore(myGenes, myGeneSet1, myUniverse, useEASE=FALSE)
#' gsFisherTestCore(myGenes, myGeneSet1, myUniverse, useEASE=TRUE)
#' 
#' ## checkUniverse will make sure that \code{univese} contains all element in \code{genes}
#' gsFisherTestCore(c("OutOfUniverse", myGenes), myGeneSet1, myUniverse, checkUniverse=FALSE)
#' gsFisherTestCore(c("OutOfUniverse", myGenes), myGeneSet1, myUniverse, checkUniverse=TRUE)
#' @importFrom ribiosUtils uniqueNonNA
#' @export
gsFisherTestCore <- function(genes, geneSetGenes, universe, 
                             makeUniqueNonNA=TRUE,
                             checkUniverse=TRUE,
                             useEASE=FALSE) {
  if(makeUniqueNonNA) {
    genes <- uniqueNonNA(genes)
    geneSetGenes <- uniqueNonNA(geneSetGenes)
    universe <- uniqueNonNA(universe)
  }
  if(checkUniverse) {
    geneNotInUniverse <- !genes %in% universe
    if(any(geneNotInUniverse)) {
      universe <- union(universe, geneNotInUniverse)
    }
  }
  
  is.pos.geneset <- genes %in% geneSetGenes
  
  geneSetGenes <- intersect(geneSetGenes, universe)
  
  pos.geneset <- sum(is.pos.geneset)
  pos.nonGeneset <- length(genes)-pos.geneset
  neg.geneset <- length(geneSetGenes)-pos.geneset
  neg.nonGeneset <- length(universe)-pos.geneset-pos.nonGeneset-neg.geneset
  
  if(useEASE) {
    pos.geneset <- pos.geneset - 1
  }
  
  mat <- matrix(c(pos.geneset, pos.nonGeneset, neg.geneset, neg.nonGeneset),
                byrow=TRUE, nrow=2)
  fisher.p <- fisher.test(mat, alternative="greater")$p.value
  hits <- genes[is.pos.geneset]
  
  res <- list(p=fisher.p,
             gsEffSize=length(geneSetGenes),
             hits=hits)
  return(res)
}

#' Core algorithm to perform Fisher's exact test on a list of gene set
#'
#' @param genes Character vector, a collection of genes of which over-representation of the gene set is tested
#' @param geneSetGenesList A list of character vector, genes belonging to each gene set
#' @param universe Character vector, universe of genes
#' @param makeUniqueNonNA Logical, whether genes, geneSetGenes, and universe should be filtered to remove NA and made unique. The default is set to \code{TRUE}. When the uniqueness and absence of NA is ensured, this flag can be set to \code{FALSE} to accelerate the operation.
#' @param checkUniverse Logical, if \code{TRUE}, then genes that are in \code{genes} but are not in \code{universe} are appended to \code{universe}
#' @param useEASE Logical, whether to use the EASE method to report the p-value. 
#'
#' This function performs one-sided Fisher's exact test to test the over-representation of the genes given as \code{geneSetGenes} in the input \code{genes} list.
#' 
#' If \code{useEASE} is \code{TRUE}, one gene is penalized (removed) within \code{geneSetGenes} that are in \code{genes} and calculating the resulting Fisher exact probability for that namespace. The theoretical basis of the EASE score lies in the concept of jackknifing a probability. See Hosack \emph{et al.} for details.
#'
#' @return A list of lists, of the same length as the input geneSetGenesList, each list consisting of three elements
#' \enumerate{
#'   \item p The p-value of one-sided (over-representation of the Fisher's test)
#'   \item gsEffSize Gene-set's effective size, namely number of genes that are in the universe
#'   \item hits Character vector, genes that are found in the gene sets
#' }
#' @references 
#' \describe{
#'   \item{Hosack \emph{et al.}}{Hosack, Douglas A., Glynn Dennis, Brad T. Sherman, H. Clifford Lane, and Richard A. Lempicki. Identifying Biological Themes within Lists of Genes with EASE. Genome Biology 4 (2003): R70. \url{https://doi.org/10.1186/gb-2003-4-10-r70}}
#' }
#' 
#' @seealso \code{\link{gsFisherTestCore}}
#' @examples
#' myGenes <- LETTERS[1:3]
#' myGeneSet1 <- LETTERS[1:6]
#' myGeneSet2 <- LETTERS[4:7]
#' myUniverse <- LETTERS
#' gsListFisherTestCore(myGenes, list(myGeneSet1, myGeneSet2), myUniverse)
#' @export
gsListFisherTestCore <- function(genes, geneSetGenesList, universe, 
                             makeUniqueNonNA=TRUE,
                             checkUniverse=TRUE,
                             useEASE=FALSE) {
  if(makeUniqueNonNA) {
    genes <- uniqueNonNA(genes)
    geneSetGenes <- lapply(geneSetGenesList, uniqueNonNA)
    universe <- uniqueNonNA(universe)
  }
  if(checkUniverse) {
    geneNotInUniverse <- !genes %in% universe
    if(any(geneNotInUniverse)) {
      universe <- union(universe, geneNotInUniverse)
    }
  }
  
  res <- lapply(geneSetGenesList, function(geneSetGenes) {
    res <- gsFisherTestCore(genes=genes,
                            geneSetGenes=geneSetGenes,
                            universe=universe,
                            makeUniqueNonNA=FALSE,
                            checkUniverse=FALSE,
                            useEASE=useEASE)
  })
  names(res) <- names(geneSetGenesList)
  return(res)
}

#' Append NewHitsProp to the result \code{data.table} returned by \code{fisherTest}
#' @param fisherTestResults \code{data.table} returned by \code{fisherTest}
#' @return A new \code{data.table} containing all columns of the input and \code{NewHitsProp}, a new column including the proportion of new hits in the gene-set
#' @importFrom dplyr `%>%` arrange
fisherTestResultNewHitsProp <- function(fisherTestResults) {
  fisherTestResults <- fisherTestResults %>% arrange(FDR)
  hits <- strsplit(as.character(fisherTestResults$Hits), ",")
  cumOC <- ribiosUtils::cumOverlapDistance(hits)
  fisherTestResults$NewHitsProp <- cumOC
  return(fisherTestResults)
}

#' Perform Fisher's exact test on a GmtList object
#' @param genes character strings of gene list to be tested
#' @param genesets An GmtList object
#' @param universe Universe (background) gene list
#' @param gsNamespace Character string, gene-set namespace(s)
#' @param makeUniqueNonNA Logical, whether genes and universe should be filtered to remove NA and made unique. The default is set to \code{TRUE}. When the uniqueness and absence of NA is ensured, this flag can be set to \code{FALSE} to accelerate the operation.
#' @param checkUniverse Logical, if \code{TRUE}, then genes that are in \code{genes} but are not in \code{universe} are appended to \code{universe}
#' @param useEASE Logical, whether to use the EASE method to report the p-value. 
#' 
#' @return A \code{data.table} containing Fisher's exact test results of all gene-sets, in the same order as the input gene-sets, with following columns:
#' \enumerate{
#'   \item GeneSetNamespace
#'   \item GeneSetName
#'   \item GeneSetEffectiveSize, the count of genes in the gene-set that are found in the universe
#'   \item HitCount, the count of genes in the \code{genes} input that are in the gene-set
#'   \item Hits, a vector of character string, representing hits
#'   \item PValue
#'   \item FDR, PValue adjusted by the Benjamini-Hochberg method. If more than one gene-set categories are provided, the FDR correction is performed per namespace
#' }
#' @examples
#' gs1 <- list(name="GeneSet1", desc="desc", genes=LETTERS[1:4], namespace="A")
#' gs2 <- list(name="GeneSet2", desc="desc", genes=LETTERS[5:8], namespace="A")
#' gs3 <- list(name="GeneSet3", desc="desc", genes=LETTERS[seq(2,8,2)], namespace="A")
#' gs4 <- list(name="GeneSet3", desc="desc", genes=LETTERS[seq(1,7,2)], namespace="B")
#' gmtList <- BioQC::GmtList(list(gs1, gs2, gs3, gs4))
#' myInput <- LETTERS[2:6]
#' myUniverse <- LETTERS
#' myFisherRes <- fisherTest(myInput, gmtList, myUniverse)
#' @export
setMethod("fisherTest", 
          c("character", "GmtList", "character"),
          function(genes, genesets, universe,
                   gsNamespace,
                   makeUniqueNonNA = TRUE,
                   checkUniverse = TRUE, useEASE = FALSE) {
            if(makeUniqueNonNA) {
              genes <- uniqueNonNA(genes)
              universe <- uniqueNonNA(universe)
            }
            if(checkUniverse) {
              geneNotInUniverse <- !genes %in% universe
              if(any(geneNotInUniverse)) {
                universe <- union(universe, geneNotInUniverse)
              }
            }
            if(missing(gsNamespace) || is.null(gsNamespace)) {
              gsNamespace <- unlist(sapply(genesets, function(x) x$namespace))
              if(is.null(gsNamespace))
                gsNamespace <- NA
            } else {
              haltifnot(length(gsNamespace) == length(genesets),
                        msg=sprintf("Length of namespace (%d) does not match length of the gene-sets (%d)", 
                                   length(gsNamespace), length(genesets)))
            }
            res <- lapply(genesets, function(x) {
              gsFisherTestCore(genes = genes, 
                               geneSetGenes = x$genes,
                               universe = universe,
                               makeUniqueNonNA = makeUniqueNonNA,
                               checkUniverse = checkUniverse,
                               useEASE = useEASE)
            })

            res <- data.table::data.table(GeneSetNamespace=gsNamespace,
                              GeneSetName=sapply(genesets, function(x) x$name),
                              GeneSetEffectiveSize=sapply(res, function(x) x$gsEffSize),
                              HitCount=sapply(res, function(x) length(x$hits)),
                              Hits=sapply(res, function(x) x$hits),
                              PValue=sapply(res, function(x) x$p))

            if(all(is.na(gsNamespace))) {
              res$FDR <- p.adjust(res$PValue, method="fdr")
            } else {
              res$FDR <- ave(res$PValue, gsNamespace, 
                             FUN=function(x) p.adjust(x, "fdr"))
            }
            return(res)
          })

utils::globalVariables(c("logFC", "Contrast", "FDR", 
			 "GeneSymbol", "GeneSetEffectiveSize",
			 "NGenes", "PValue", "GeneID", "GeneSymbol", 
			 "GeneSet", "Regulation"))

#' Run Fisher's exact test on an EdgeResult object
#' 
#' @param edgeResult An \code{EdgeResult} object
#' @param gmtList A \code{GmtList} or \code{GeneSets} object
#' @param contrast Character, the contrast of interest
#' @param thr.abs.logFC Numeric, threshold of absolute log2 fold-change to 
#'     define positively and negatively regulated genes
#' @param thr.FDR Numeric, threshold of FDR values 
#' @param minGeneSetEffectiveSize Integer, minimal number of genes of a 
#'   geneset that are quantified
#' @param maxGeneSetEffectiveSize Integer, maximal number of genes of a 
#'    geneset that are quantified
#' @param ... Passed to \code{filter} to further filter the differential gene
#' expression table (\code{dgeTbl}).
#' 
#' @importFrom ribiosNGS dgeTable
#' @importFrom dplyr filter pull `%>%` ungroup arrange group_by
#' @return 
#' ## TODO: example
fisherTestEdgeResult <- function(edgeResult,
                             gmtList,
                             contrast, 
                             thr.abs.logFC=1, thr.FDR=0.05, 
                             minGeneSetEffectiveSize=5,
                             maxGeneSetEffectiveSize=500, ...) {
  dgeTbl <- ribiosNGS::dgeTable(edgeResult) %>% filter(Contrast %in% contrast)
  dgeBg <- unique(as.character(dgeTbl$GeneSymbol))
  posDgeTbl <- dgeTbl %>% filter(logFC>=thr.abs.logFC, FDR<thr.FDR, ...)
  negDgeTbl <- dgeTbl %>% filter(logFC<=(-thr.abs.logFC), FDR<thr.FDR, ...)
  posGenes <- posDgeTbl %>% pull(GeneSymbol) %>% as.character
  negGenes <- negDgeTbl %>% pull(GeneSymbol) %>% as.character
  posFisher <- fisherTest(posGenes, gmtList, dgeBg)
  negFisher <- fisherTest(negGenes, gmtList, dgeBg)
  
  posFisherHits <- posFisher %>% filter(FDR<=thr.FDR,
                                        GeneSetEffectiveSize>=minGeneSetEffectiveSize,
                                        GeneSetEffectiveSize<=maxGeneSetEffectiveSize)
  negFisherHits <- negFisher %>% filter(FDR<=0.05,
                                        GeneSetEffectiveSize>=minGeneSetEffectiveSize,
                                        GeneSetEffectiveSize<=maxGeneSetEffectiveSize)
  
  res <- data.table::data.table(Regulation=factor(rep(c("Positive", "Negative"), 
                                          c(nrow(posFisherHits), nrow(negFisherHits))),
                                      c("Positive", "Negative")),
                    rbind(posFisherHits, negFisherHits)) %>% 
    group_by(Regulation) %>% 
    fisherTestResultNewHitsProp %>%
    ungroup %>%
    arrange(Regulation, FDR)
  return(res)
}
