% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fisher.R
\name{gsListFisherTestCore}
\alias{gsListFisherTestCore}
\title{Core algorithm to perform Fisher's exact test on a list of gene set}
\usage{
gsListFisherTestCore(
  genes,
  geneSetGenesList,
  universe,
  makeUniqueNonNA = TRUE,
  checkUniverse = TRUE,
  useEASE = FALSE
)
}
\arguments{
\item{genes}{Character vector, a collection of genes of which over-representation of the gene set is tested}

\item{geneSetGenesList}{A list of character vector, genes belonging to each gene set}

\item{universe}{Character vector, universe of genes}

\item{makeUniqueNonNA}{Logical, whether genes, geneSetGenes, and universe should be filtered to remove NA and made unique. The default is set to \code{TRUE}. When the uniqueness and absence of NA is ensured, this flag can be set to \code{FALSE} to accelerate the operation.}

\item{checkUniverse}{Logical, if \code{TRUE}, then genes that are in \code{genes} but are not in \code{universe} are appended to \code{universe}}

\item{useEASE}{Logical, whether to use the EASE method to report the p-value. 

This function performs one-sided Fisher's exact test to test the over-representation of the genes given as \code{geneSetGenes} in the input \code{genes} list.

If \code{useEASE} is \code{TRUE}, one gene is penalized (removed) within \code{geneSetGenes} that are in \code{genes} and calculating the resulting Fisher exact probability for that namespace. The theoretical basis of the EASE score lies in the concept of jackknifing a probability. See Hosack \emph{et al.} for details.}
}
\value{
A list of lists, of the same length as the input geneSetGenesList, each list consisting of three elements
\enumerate{
  \item p The p-value of one-sided (over-representation of the Fisher's test)
  \item gsEffectiveSize Gene-set's effective size, namely number of genes that are in the universe
  \item hits Character vector, genes that are found in the gene sets
}
}
\description{
Core algorithm to perform Fisher's exact test on a list of gene set
}
\examples{
myGenes <- LETTERS[1:3]
myGeneSet1 <- LETTERS[1:6]
myGeneSet2 <- LETTERS[4:7]
myUniverse <- LETTERS
gsListFisherTestCore(myGenes, list(myGeneSet1, myGeneSet2), myUniverse)
}
\references{
\describe{
  \item{Hosack \emph{et al.}}{Hosack, Douglas A., Glynn Dennis, Brad T.
  Sherman, H. Clifford Lane, and Richard A. Lempicki. Identifying Biological
  Themes within Lists of Genes with EASE. Genome Biology 4 (2003): R70. \doi{10.1186/gb-2003-4-10-r70}}
}
}
\seealso{
\code{\link{gsFisherTestCore}}
}
