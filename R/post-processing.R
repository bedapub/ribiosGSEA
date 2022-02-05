## for post-processing

#' Parse contributing genes from the CAMERA output file
#'
#' @param str Character string, containing contributing genes
#' @return A list of \code{data.frames}, each containing two columns,
#' \code{Gene} and \code{Stat}
#' @examples
#' 
#' parseContributingGenes("AKR1C4(-1.25), AKR1D1(-1.11)")
#' parseContributingGenes(c("AKR1C4(-1.25), AKR1D1(-1.11)",
#'                          "AKT1(1.24), AKT2(1.11), AKT3(1.05)"))
#' 
#' @export parseContributingGenes
parseContributingGenes <- function(str) {
  ss <- strsplit(str, ",")
  res <- lapply(ss, function(x) {
    genes <- gsub("\\(.*\\)", "", x)
    stat <- as.numeric(as.character(gsub(".*\\((.*)\\)", "\\1", x)))
    return(data.frame(Gene=genes, Stat=stat))
  })
  return(res)
}



#' Parse contributing genes by genesets
#' 
#' @param str Character strings, containing contributing genes
#' @param genesets Character strings, geneset labels. Its length must match the
#' length of \code{str}
#' @return A \code{data.frame} containing genesets, genes, and statistics
#' @examples
#' 
#' parseGenesetsContributingGenes("AKR1C4(-1.25), AKR1D1(-1.11)", "Metabolism")
#' parseGenesetsContributingGenes(c("AKR1C4(-1.25), AKR1D1(-1.11)",
#'                          "AKT1(1.24), AKT2(1.11), AKT3(1.05)"),
#'                          c("Metabolism", "AKTs"))
#' 
#' @export parseGenesetsContributingGenes
parseGenesetsContributingGenes <- function(str, genesets) {
  stopifnot(length(str)==length(genesets))
  genes <- parseContributingGenes(str)
  res <- cbind(GeneSet=rep(genesets, sapply(genes, nrow)),
               do.call(rbind, genes))
  return(res)
}

#' Parse contributing genes by genesets from the result data.frame of the \code{CAMERA} method
#' 
#' @param cameraResTbl A \code{tibble} or \code{data.frame} holding results of \code{CAMERA}
#' @param genesets Character strings, geneset labels
#' 
#' @return A list of gene symbols, indexed by geneset names that are found in the results.
#' 
#' @export parseCameraContributingGenes
parseCameraContributingGenes <- function(cameraResTbl, genesets) {
  res <- cameraResTbl %>% dplyr::filter(GeneSet %in% genesets) %>% 
    (function(x) parseGenesetsContributingGenes(x$ContributingGenes, x$GeneSet)) %>%
    (function(x) split(as.character(x$Gene), x$GeneSet))
  res <- lapply(res[genesets], unique)
  return(res)
}


#' Read CAMERA results into a tibble object
#' @param file CAMERA results file
#' @param minNGenes NULL or integer, genesets with fewer genes are filtered out
#' @param maxNGenes NULL or integer, genesets with more genes are filtered out
#' @importFrom readr read_tsv
#' @export readCameraResults
readCameraResults <- function(file, minNGenes=3, maxNGenes=1000) {
  res <- readr::read_tsv(file, col_types = "cccicddddddddc")
  if(!is.null(minNGenes)) {
    res <- subset(res, NGenes>=minNGenes)
  }
  if(!is.null(maxNGenes)) {
    res <- subset(res, NGenes<=maxNGenes)
  }
  return(res)
}

#' Expand genes in the CAMERA result table 
#' @param tbl A \code{data.frame}
#' @return A longer \code{data.frame}, with each row one gene.
#' @export
expandCameraTableGenes <- function(tbl) {
    genes <- tbl$ContributingGenes
    expGenes <- parseContributingGenes(genes)
    isRemove <- colnames(tbl) %in% "ContributingGenes"
    nGenes <- sapply(expGenes, nrow)
    pathTbl <- tbl[rep(1:nrow(tbl),  nGenes), !isRemove]
    geneTbl <- do.call(rbind, expGenes)
    res <- cbind(pathTbl, geneTbl)
    rownames(res) <- NULL
    return(res)
}

#' Convert a CAMERA table into a graph
#'
#' @param df Data.frame, CAMERA results
#' @param jacThr Numeric, between 0 and 1, Jaccard Index threshold
#' @param plot Logical, whether plotting the results
#' @param ... Passed to \code{plot}
#'
#' @importFrom igraph make_empty_graph graph_from_adjacency_matrix 
#'     layout_in_circle V
#' @importFrom ribiosUtils pairwiseJaccardIndex matchColumn boundNorm
#' @importFrom ribiosPlot royalbluegrayred
#' @export
cameraTable2graph <- function(df, jacThr=0.25, plot=TRUE, ...) {
    retObj <- list(graph=igraph::make_empty_graph(),
                   resTbl=data.frame(Namespace=character(0),
                        GeneSet=character(0),
                       score=numeric(0)))
    if(nrow(df)==0) {
        return(retObj)
    }
    
    scores <- df$Score
    labels <- df$label <- with(df, paste(Namespace,":",GeneSet,sep=""))
                               
    geneLists <- with(df, split(as.character(Gene), labels))
    gJacSim <- ribiosUtils::pairwiseJaccardIndex(geneLists)

    gJacBin <- gJacSim
    gJacBin[gJacBin>=jacThr] <- 1L
    gJacBin[gJacBin<jacThr] <- 0L

    gLabels <- names(geneLists)
    rownames(gJacBin) <- colnames(gJacBin) <- gLabels

    gNamespace <- sapply(strsplit(rownames(gJacBin), ":"), "[[", 1L)
    gGeneSet <- sapply(strsplit(rownames(gJacBin), ":"), function(x) paste(x[2:length(x)], collapse=" "))
    ugNamespace <- unique(gNamespace)
    for(ugc in ugNamespace) {
        if(grepl("^MPS",ugc)) {  ## MPS specific!
            isUgc <- gNamespace==ugc
            gJacBin[isUgc, isUgc] <- 1L
        }
    }
        
    graph <- igraph::graph_from_adjacency_matrix(gJacBin, mode="undirected", diag=FALSE)
    
    graphVscores <- matchColumn(igraph::V(graph)$name, df, "label")$Score
    
    absMaxScore <- max(abs(graphVscores))
    score2range <- as.integer(ribiosUtils::boundNorm(graphVscores)*100,
                              low=-absMaxScore, high=absMaxScore)
    scoreColor <- royalbluegrayred(101)[score2range+1]

    gOrder <- order(graphVscores)
    layout <- layout_in_circle(graph, order=gOrder)
    if(plot) {
        plot(graph,
             layout=layout,
             vertex.color=scoreColor,
             vertex.label=gsub(":", "\n", gsub("^MPS ", "", gLabels)),  ## MPS specific!
             vertex.label.cex=0.8,
             edge.arrow.size=2, ...)
    }
    retObj$graph <- graph
    retObj$resTbl <- data.frame(Namespace=gNamespace,
                                GeneSet=gGeneSet,
                                Score=graphVscores)
    return(retObj)
}

#' Read significant CAMERA results into a tibble
#' @param file A tsv file, output of \code{\link{biosCamera}}
#' @param returnAllContrasts Logical, if TRUE, results of all contrasts for gene-sets that are significant in at least one contrast are returned.
#' @param maxPValue Numeric, max unadjusted P-value of CAMERA that is considered significant
#' @param minAbsEffectSize Numeric, minimal absolute effect size
#' @param minNGenes Integer, size of the smallest gene set that is considered
#' @param maxNGenes Integer, size of the largest gene set that is considered
#' @param excludeNamespace Character, vector of namespaces to be excluded
#' @importFrom dplyr `%>%` filter
#' @export
readSigCameraResults <- function(file, 
                                 returnAllContrasts=TRUE,
                                 maxPValue=0.01, 
                                 minAbsEffectSize=0.5,
                                 minNGenes=5, 
                                 maxNGenes=200,
                                 excludeNamespace=c("goslim", "immunespace", 
                                                    "immunomics", "mbdisease",
                                                    "mbpathology", "mbtoxicity",
                                                    "msigdbC7", "msigdbC2",
                                                    "MolecularPhenotyping")) {
  EffectSize <- Namespace <- NULL
  cameraRes <- readCameraResults(file)
  cameraSig <- cameraRes %>% filter(PValue<maxPValue,
                                    NGenes>=minNGenes,
                                    NGenes<=maxNGenes,
                                    abs(EffectSize)>=minAbsEffectSize,
                                    ! Namespace %in% excludeNamespace)
  if(returnAllContrasts) {
    res <- cameraRes %>% filter(GeneSet %in% cameraSig$GeneSet) 
  } else {
    res <- cameraSig
  }
  return(res)
}

#' Read significant CAMERA results into a matrix
#' @param file A tsv file, output of \code{\link{biosCamera}}
#' @param ... passed to \code{readSigCameraResults}
#' @importFrom ribiosUtils longdf2matrix
#' @export
readSigCameraScoreMatrix <- function(file, ...) {
  sigRes <- readSigCameraResults(file, returnAllContrasts = TRUE, 
                                 ...)
  res <- sigRes %>% 
    longdf2matrix(row.col="GeneSet", column.col="Contrast", value.col="Score")
  return(res)
}

## #' @export
## expandSigCameraResults <- function(cameraTable,
##                                    pVal=0.01, minNGenes=3, includeAllFromSigGenesets=FALSE) {
##     if(includeAllFromSigGenesets) {
##         sigGeneSets <- subset(cameraTable, PValue<=pVal)$GeneSet
##         sigCameraTable <- subset(cameraTable, GeneSet %in% sigGeneSets & NGenes >= minNGenes)
##     } else {
##         sigCameraTable <- subset(cameraTable, PValue<=pVal & NGenes >= minNGenes)
##     }
##     expSigCameraTable <- expandCameraTableGenes(sigCameraTable)
##     return(expSigCameraTable)
## }
## 
## #' @export
## visualizeCameraGraphsByContrast <- function(expCameraTable, ...) {
##     contrasts <- sort(unique(expCameraTable$Contrast))
##     nres <- lapply(as.character(contrasts), function(x) {
##                        subexpCameraTable <- subset(expCameraTable, Contrast==x)
##                        cameraTable2graph(subexpCameraTable, plot=TRUE, main=x, ...)
##                   })
##     tbls <- lapply(nres, function(x) x$resTbl)
##     res <- cbind(Contrast=rep(contrasts, sapply(tbls, nrow)),
##                  do.call(rbind, tbls))
##                       
##     return(res)
## }
