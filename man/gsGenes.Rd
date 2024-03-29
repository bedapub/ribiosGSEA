% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/AllMethods.R
\name{gsGenes}
\alias{gsGenes}
\alias{gsGenes,AnnoBroadGseaResItem-method}
\alias{gsGenes,AnnoBroadGseaRes-method}
\alias{gsGenes,GmtList-method}
\title{Return gene-set genes}
\usage{
gsGenes(object, ...)

\S4method{gsGenes}{AnnoBroadGseaResItem}(object)

\S4method{gsGenes}{AnnoBroadGseaRes}(object)

\S4method{gsGenes}{GmtList}(object)
}
\arguments{
\item{object}{An object}

\item{...}{Other parameters}
}
\description{
Return gene-set genes
}
\section{Methods (by class)}{
\itemize{
\item \code{gsGenes(AnnoBroadGseaResItem)}: Get gene-set genes
from a BroadGseaResItem object, returning a character string vector.

\item \code{gsGenes(AnnoBroadGseaRes)}: Get gene-set genes from an AnnoBroadGseaRes object,
returning a list of character string vectors.

\item \code{gsGenes(GmtList)}: Get gene-set genes from a GmtList object, returning a
list of character string vector. It uses the implementation in BioQC.

}}
