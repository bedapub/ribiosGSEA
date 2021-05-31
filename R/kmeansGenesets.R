#' Cluster gene-sets by enrichment profiles with k-means clustering, and select representative gene-sets by gene-set composition
#' 
#' @param enrichProfMatrix A numeric matrix representing gene-set enrichment profile. Each row represent one gene-set and each column represent one enrichment profile, for instance a contrast in differential gene expression analysis. The values of the matrix represent enrichment of gene-sets, for instance enrichment score or absolute log10-transform p-values can be used. The row names are gene-set names.
#' @param genesetGenes A list of character strings, each element being genes of a gene-set in the \code{enrichProfMatrix}. The names of the list must exactly match the row-names of \code{enrichProfMatrix}, namely the names of gene-sets in the same order.
#' @param optK Integer, the number of initial clusters of gene-sets. Because one or more gene-sets may be selected from each gene-set cluster, the number of finally selected gene-sets is equal to or larger than \code{optK}.
#' @param iter.max Integer, the maximum numbers of iterations allowed. This parameter is passed to \code{\link[stats]{kmeans}}.
#' @param nstart Integer, how many random sets should be chosen to initialize cluster centers. This parameter is passed to \code{\link[stats]{kmeans}}.
#' @param thrCumJaccardIndex Numeric, between 0 and 1, the threshold of cumulative Jaccard Index. The larger the value is, the more gene-sets will be selected from each cluster
#' @param maxRepPerCluster Integer, maximum number of representative genesets per cluster. If NULL or NA, no limit is set.
#' @param metaClusterColumns Columns used to cluster the clusters by their average enrichment profile. By default, all columns are used.
#' 
#' This function performs \code{k-means} clustering of enrichment profiles of gene-sets. Within each cluster, we first identify the union set of unique genes covered any gene-set in the cluster, and then calculate Jaccard Index between genes in each gene-set and the union set. Gene-sets are sorted descendingly by the Jaccard Index, and the cumulative Jaccard Index is calculated. Among the sorted gene-sets, the gene-sets up to the position when the cumulative Jaccard Index exceeds \code{thrCumJaccardIndex} are selected (excluding redundant gene-sets).
#' 
#' The geneset clusters are ordered by their average profiles - similar clusters are near to each other.
#' 
#' @return A list:
#' \itemize{
#'  \item kmeans Result object returned by \code{kmeans}. 
#'  \item genesetClusterData A \code{data.frame} with following columns: \code{GenesetCluster}, \code{GenesetInd}, \code{GenesetName}, \code{JaccardIndex}, \code{CumJaccardIndex}, \code{IsRepresentative}.
#'  \item repGenesets Character vector, gene-set names that are selected as representative gene-sets from each gene-set clsuter.
#'  \item gsCompOverlapSelInd Factor vector, indicating the gene-set clusters represented by each representative gene-set.
#' }
#' 
#' @importFrom ribiosUtils munion jaccardIndex
#' @importFrom stats hclust dist
#' 
#' @examples
#' set.seed(1887)
#' profMat <- matrix(rnorm(100), nrow=20, 
#'     dimnames=list(sprintf("geneset%d", 1:20), sprintf("contrast%d", 1:5)))
#' gsGenes <- lapply(1:nrow(profMat), function(x) 
#'     unique(sample(LETTERS, 10, replace=TRUE)))
#' names(gsGenes) <- rownames(profMat)
#' kmeansGeneset(profMat, gsGenes, optK=5)
#' 
#' @export
kmeansGeneset <- function(enrichProfMatrix, genesetGenes,
                          optK=pmin(25,floor(nrow(enrichProfMatrix)/2)),
                          iter.max=15, nstart=50, 
                          thrCumJaccardIndex=0.5,
                          maxRepPerCluster=10,
                          metaClusterColumns=1:ncol(enrichProfMatrix)) {
  stopifnot(nrow(enrichProfMatrix)==length(genesetGenes) &&
              identical(rownames(enrichProfMatrix),
                        names(genesetGenes)))
  gsNames <- as.character(rownames(enrichProfMatrix))
  maxOptK <- floor(nrow(enrichProfMatrix)/2)
  if(optK>maxOptK) {
    optK <- maxOptK
  }
  kmeansObj <- kmeans(enrichProfMatrix, optK, nstart=nstart, iter.max=iter.max)
  pathCluster <- kmeansObj$cluster
  pathClusterFac <- factor(pathCluster); levels(pathClusterFac) <- sprintf("GenesetCluster%s", levels(pathClusterFac))
  gsComps <- split(seq(along=genesetGenes), pathClusterFac)
  gsCompOverlap <- lapply(seq(along=gsComps), function(i) {
    ind <- gsComps[[i]]
    clusterUnionGenes <- ribiosUtils::munion(genesetGenes[ind])
    coefs <- sapply(genesetGenes[ind], function(genes) jaccardIndex(genes, clusterUnionGenes))
    coefOrd <- order(coefs, decreasing = TRUE)
    ind <- ind[coefOrd]
    coefs <- coefs[coefOrd]
    genesNewList <- genesetGenes[ind]
    cumJac <- sapply(seq(along=genesNewList), function(j) jaccardIndex(ribiosUtils::munion(genesNewList[1:j]), 
                                                                       clusterUnionGenes))

    res <- data.frame(GenesetCluster=names(gsComps)[i],
                      GenesetInd=ind, 
                      GenesetName=gsNames[ind],
                      JaccardIndex=coefs, 
                      CumJaccardIndex=cumJac,
                      row.names = NULL)
    
    hasDiff <- c(TRUE, diff(res$CumJaccardIndex)!=0) ## ignore redudant genesets
    indFirstOverThr <- min(which(res$CumJaccardIndex>thrCumJaccardIndex &
                                   hasDiff))
    if(!is.null(maxRepPerCluster) && !is.na(maxRepPerCluster))
      indFirstOverThr <- pmin(indFirstOverThr, as.integer(maxRepPerCluster))
    sel <- c(rep(TRUE, indFirstOverThr),
             rep(FALSE, nrow(res)-indFirstOverThr))
    res$IsRepresentative <- sel
    
    return(res)
  })

  gsCompDf <- do.call(rbind, gsCompOverlap)
  gsCompOverlapSels <- with(gsCompDf,
                            as.character(GenesetName[IsRepresentative]))
  gsCompOverlapSelInd <- with(gsCompDf,
                              GenesetCluster[IsRepresentative])
  ## re-arrange the levels of the clusters so clusters with similar profiles
  ## cluster together
  repRows <- which(gsCompDf$IsRepresentative)
  clusterAvgProfList <- tapply(repRows,
                               gsCompOverlapSelInd,
                               function(i) {
                                 gsInd <- gsCompDf$GenesetInd[i]
                                 res <- colMeans(enrichProfMatrix[gsInd,
                                                                  metaClusterColumns,
                                                                  drop=FALSE],
                                                 na.rm=TRUE)
                               })
  if(is.list(clusterAvgProfList)) {
    clusterAvgProf <- do.call(rbind, clusterAvgProfList)
  } else {
    clusterAvgProf <- matrix(clusterAvgProfList, ncol=1,
                             dimnames = list(names(clusterAvgProfList), NULL))
  }
  clusterHclust <- stats::hclust(stats::dist(clusterAvgProf), 
                                 "ward.D2")
  clusterOrder <- clusterHclust$label[clusterHclust$order]
  
  gsClusters <- factor(gsCompOverlapSelInd,
                       levels=clusterOrder)
  levels(gsClusters) <- sprintf("GenesetCluster%s", 1:nlevels(gsClusters))
  gsClusterOrder <- order(gsClusters)
  
  orderedRepGenesets <- gsCompOverlapSels[gsClusterOrder]
  orderedRepGenesetClusters <- gsClusters[gsClusterOrder]
  
  res <- list(kmeans=kmeansObj,
              genesetClusterData=gsCompDf,
              repGenesets=orderedRepGenesets,
              repGenesetClusters=orderedRepGenesetClusters)
  return(res)
}

assignInNamespace("kmeansGeneset", kmeansGeneset, "ribiosGSEA")
