% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/post-processing.R
\name{readSigCameraResults}
\alias{readSigCameraResults}
\title{Read significant CAMERA results into a tibble}
\usage{
readSigCameraResults(
  file,
  returnAllContrasts = TRUE,
  maxPValue = 0.01,
  minAbsEffectSize = 0.5,
  minNGenes = 5,
  maxNGenes = 200,
  excludeNamespace = c("goslim", "immunespace", "immunomics", "mbdisease", "mbpathology",
    "mbtoxicity", "msigdbC7", "msigdbC2", "MolecularPhenotyping")
)
}
\arguments{
\item{file}{A tsv file, output of \code{\link{biosCamera}}}

\item{returnAllContrasts}{Logical, if TRUE, results of all contrasts for gene-sets that are significant in at least one contrast are returned.}

\item{maxPValue}{Numeric, max unadjusted P-value of CAMERA that is considered significant}

\item{minAbsEffectSize}{Numeric, minimal absolute effect size}

\item{minNGenes}{Integer, size of the smallest gene set that is considered}

\item{maxNGenes}{Integer, size of the largest gene set that is considered}

\item{excludeNamespace}{Character, vector of namespaces to be excluded}
}
\description{
Read significant CAMERA results into a tibble
}
