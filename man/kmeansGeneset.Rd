% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kmeansGenesets.R
\name{kmeansGeneset}
\alias{kmeansGeneset}
\title{Cluster gene-sets by enrichment profiles with k-means clustering, and select representative gene-sets by gene-set composition}
\usage{
kmeansGeneset(
  enrichProfMatrix,
  genesetGenes,
  optK = pmin(25, floor(nrow(enrichProfMatrix)/2)),
  iter.max = 15,
  nstart = 50,
  thrCumJaccardIndex = 0.5,
  maxRepPerCluster = 10,
  metaClusterColumns = 1:ncol(enrichProfMatrix)
)
}
\arguments{
\item{enrichProfMatrix}{A numeric matrix representing gene-set enrichment profile. Each row represent one gene-set and each column represent one enrichment profile, for instance a contrast in differential gene expression analysis. The values of the matrix represent enrichment of gene-sets, for instance enrichment score or absolute log10-transform p-values can be used. The row names are gene-set names.}

\item{genesetGenes}{A list of character strings, each element being genes of a gene-set in the \code{enrichProfMatrix}. The names of the list must exactly match the row-names of \code{enrichProfMatrix}, namely the names of gene-sets in the same order.}

\item{optK}{Integer, the number of initial clusters of gene-sets. Because one or more gene-sets may be selected from each gene-set cluster, the number of finally selected gene-sets is equal to or larger than \code{optK}.}

\item{iter.max}{Integer, the maximum numbers of iterations allowed. This parameter is passed to \code{\link[stats]{kmeans}}.}

\item{nstart}{Integer, how many random sets should be chosen to initialize cluster centers. This parameter is passed to \code{\link[stats]{kmeans}}.}

\item{thrCumJaccardIndex}{Numeric, between 0 and 1, the threshold of cumulative Jaccard Index. The larger the value is, the more gene-sets will be selected from each cluster}

\item{maxRepPerCluster}{Integer, maximum number of representative genesets per cluster. If NULL or NA, no limit is set.}

\item{metaClusterColumns}{Columns used to cluster the clusters by their average enrichment profile. By default, all columns are used.

This function performs \code{k-means} clustering of enrichment profiles of gene-sets. Within each cluster, we first identify the union set of unique genes covered any gene-set in the cluster, and then calculate Jaccard Index between genes in each gene-set and the union set. Gene-sets are sorted descendingly by the Jaccard Index, and the cumulative Jaccard Index is calculated. Among the sorted gene-sets, the gene-sets up to the position when the cumulative Jaccard Index exceeds \code{thrCumJaccardIndex} are selected (excluding redundant gene-sets).

The geneset clusters are ordered by their average profiles - similar clusters are near to each other.}
}
\value{
A list:
\itemize{
 \item kmeans Result object returned by \code{kmeans}. 
 \item genesetClusterData A \code{data.frame} with following columns: \code{GenesetCluster}, \code{GenesetInd}, \code{GenesetName}, \code{JaccardIndex}, \code{CumJaccardIndex}, \code{IsRepresentative}.
 \item repGenesets Character vector, gene-set names that are selected as representative gene-sets from each gene-set clsuter.
 \item gsCompOverlapSelInd Factor vector, indicating the gene-set clusters represented by each representative gene-set.
}
}
\description{
Cluster gene-sets by enrichment profiles with k-means clustering, and select representative gene-sets by gene-set composition
}
\examples{
set.seed(1887)
profMat <- matrix(rnorm(100), nrow=20, 
    dimnames=list(sprintf("geneset\%d", 1:20), sprintf("contrast\%d", 1:5)))
gsGenes <- lapply(1:nrow(profMat), function(x) 
    unique(sample(LETTERS, 10, replace=TRUE)))
names(gsGenes) <- rownames(profMat)
kmeansGeneset(profMat, gsGenes, optK=5)

}
