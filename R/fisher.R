#' @include AllGenerics.R AllMethods.R ribiosGSEA-package.R
NULL

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


#' Perform Fisher's exact test on a gene set
#'
#' @param genes a collection of genes of which over-representation of the gene set is tested
#' @param geneSetGenes genes belonging to a gene set
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
#' fisherTest(myGenes, myGeneSet1, myUniverse)
#' fisherTest(myGenes, myGeneSet2, myUniverse)
#' fisherTest(myGenes, myGeneSet1, myUniverse, gsName="My gene set1", gsNamespace="Letters")
#'
#' ## note that duplicated items are removed by default
#' resWoRp <- fisherTest(rep(myGenes,2), myGeneSet1, myUniverse)
#' resWithRp <- fisherTest(rep(myGenes,2), myGeneSet1, rep(myUniverse,2))
#' identical(resWoRp, resWithRp)
#' 
#' resWithRpNoUnique <- fisherTest(rep(myGenes,2), myGeneSet1, rep(myUniverse,2), makeUniqueNonNA=FALSE)
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
#' @param GmtList A \code{GmtList} object
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
#' myS4GeneSet1 <- list(name="GeneSet1", desc="GeneSet", genes=LETTERS[1:6], namespace="My namespace 1")
#' myS4GeneSet2 <- list(name="GeneSet1", desc="GeneSet", genes=LETTERS[2:7], namespace="My namespace 2")
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


#' Append NewHitsProp to the result \code{data.table} returned by \code{fisherTest}
#' @param fisherTestResults \code{data.table} returned by \code{fisherTest}
#' @return A new \code{data.table} containing all columns of the input and \code{NewHitsProp}, a new column including the proportion of new hits in the gene-set
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



#' Run Fisher's exact test on an EdgeResult object
#' 
#' @param edgeResult An \code{EdgeResult} object
#' @param gmtList A \code{GmtList} or \code{GeneSets} object
#' @param contrast Character, the contrast of interest
#' @param thr.abs.logFC Numeric, threshold of absolute log2 fold-change to define positively and negatively regulated genes
#' @param thr.FDR Numeric, threshold of FDR values 
#' @param minGeneSetEffectiveSize Integer, minimal number of genes of a geneset that are quantified
#' @param maxGeneSetEffectiveSize Integer, maximal number of genes of a geneset that are quantified
#' 
#' @return 
#' ## TODO: example
fisherTestEdgeResult <- function(edgeResult,
                             gmtList,
                             contrast, 
                             thr.abs.logFC=1, thr.FDR=0.05, 
                             minGeneSetEffectiveSize=5,
                             maxGeneSetEffectiveSize=500, ...) {
  dgeTbl <- dgeTable(edgeResult) %>% filter(Contrast %in% contrast)
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
