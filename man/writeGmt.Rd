% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gmt.R
\name{writeGmt}
\alias{writeGmt}
\title{Write an GeneSet object into a file}
\usage{
writeGmt(geneSet, file)
}
\arguments{
\item{geneSet}{A GeneSet object}

\item{file}{Character string, output file name}
}
\description{
Write an GeneSet object into a file
}
\examples{
gmtFile <- system.file("extdata", "example.gmt", package="ribiosGSEA")
mySet <- readGmt(gmtFile)[1:5]
myTempFile <- tempfile()
writeGmt(mySet, file=myTempFile)
readLines(myTempFile)
}
