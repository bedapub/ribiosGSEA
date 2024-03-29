% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/AllMethods.R
\name{gseaCoreEnrichGenes}
\alias{gseaCoreEnrichGenes}
\alias{gseaLeadingEdgeGenes}
\alias{gseaCoreEnrichGenes,AnnoBroadGseaResItem-method}
\alias{gseaCoreEnrichGenes,AnnoBroadGseaRes-method}
\title{Return GSEA core enrichment genes (also known as leading-edge genes)}
\usage{
gseaCoreEnrichGenes(object)

\S4method{gseaCoreEnrichGenes}{AnnoBroadGseaResItem}(object)

\S4method{gseaCoreEnrichGenes}{AnnoBroadGseaRes}(object)
}
\arguments{
\item{object}{An object}
}
\description{
Return GSEA core enrichment genes (also known as leading-edge genes)
}
\section{Methods (by class)}{
\itemize{
\item \code{gseaCoreEnrichGenes(AnnoBroadGseaResItem)}: Return core enriched genes (also known as
leading-edge genes) in an AnnoBroadGseaResItem object as a character string
vector.

\item \code{gseaCoreEnrichGenes(AnnoBroadGseaRes)}: Return core enriched genes (also known as
leading-edge genes) in an AnnoBroadGseaRes object as a list of character 
string vectors.

}}
