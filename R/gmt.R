#' Write an GmtList object into a file
#'
#' @param gmtList A \code{GmtList} object
#' @param file Character string, output file name
#'
#' @note 
#' The function will be moved to BioQC once the ribiosIO is reposited in CRAn
#' 
#' @importFrom ribiosIO write_gmt
#' @importFrom BioQC readGmt
#' @examples
#' gmtFile <- system.file("extdata", "example.gmt", package="ribiosGSEA")
#' mySet <- BioQC::readGmt(gmtFile)[1:5]
#' myTempFile <- tempfile()
#' writeGmt(mySet, file=myTempFile)
#' readLines(myTempFile)
writeGmt <- function(gmtList, file) {
  ribiosIO::write_gmt(gmtList, file=file)
}

#' Read molecular-phenotyping genesets
#'
#' @param file GMT file which stores default molecular-phenotyping genesets
#'
#' @importFrom BioQC readGmt
#' @return A \code{GmtList} object containing molecular-phenotypic screening (MPS) categories and genes
readMPSGmt <- function(file) {
  gs <- BioQC::readGmt(file)
  namespace <- sprintf("MPS %s", BioQC::gsDesc(gs))
  gsNamespace(gs) <- namespace
  return(gs)
}

  
#' Read default genesets for gene-set enrichment analysis
#'
#' In Roche Bioinformatics we use a default collection of gene-sets for gene-set enrichment analysis. This function loads this collection.
#'
#' @param path Character, path to the directory where the gmt files are stored
#' @param mps Logical,  whether molecular-phenotypic screening (MPS) genesets should be read in as pathway-centric namespaces (\code{TRUE}) or as one namespace named \code{MolecularPhenotyping} (\code{FALSE}).
#' 
#' @details
#'
#' The default collection includes both publicly available genesets as well as proprietary genesets, and therefore they are not included as part of the ribios package.
#'
#' Publicly available genesets include
#' \itemize{
#'  \item{MSigDB: collections C2, C7 and Hallmark}
#'  \item{RONET: which is a collection of publicly available pathway databases including REACTOME and NCI-Nature}
#'  \item{goslim}
#' }
#' @importFrom BioQC appendGmtList readGmt
#' @importFrom ribiosUtils assertFile
#' @examples 
#' \dontrun{
#'   readDefaultGenesets("/tmp/defaultGmts")
#' }
readDefaultGenesets <- function(path,
                                mps=FALSE) {
  assertFile(msigdb.c2.file <- file.path(path, "msigdb.c2.all.symbols.gmt"))
  assertFile(msigdb.c7.file <- file.path(path, "msigdb.c7.all.symbols.gmt"))
  assertFile(msigdb.hallmark.file <- file.path(path, "msigdb.hallmark.all.symbols.gmt"))
  
  assertFile(mps.pathway <- file.path(path, "MolecularPhenotyping-genesets.gmt"))
  
  assertFile(ronet.file <- file.path(path, "path.ronet.roche.symbols.gmt"))
  assertFile(goslim.file <- file.path(path, "go_slim.bp.roche.symbols.gmt"))

  assertFile(upstream.file <- file.path(path, "MetaBase.downstream.expression.gmt"))  
  assertFile(mbdisease.file <- file.path(path, "MetaBase.DiseaseBiomarker.gmt"))
  assertFile(mbmetabolic.file <- file.path(path, "MetaBase.Metabolic.gmt"))
  assertFile(mbpath.file <- file.path(path, "MetaBase.PathwayMap.gmt"))
  assertFile(mbtoxicity.file <- file.path(path, "MetaBase.Toxicity.gmt"))
  assertFile(mbpathology.file <- file.path(path, "MetaBase.ToxicPathology.gmt"))
  assertFile(mbprocess.file <- file.path(path, "MetaBase.TRprocesses.gmt"))

  assertFile(immunomics.file <- file.path(path, "exp.immune.roche.symbols.gmt"))
  assertFile(immunespace.file <- file.path(path, "immunespace.gmt"))
  
  if(mps) {
    mpsGSCs <- readMPSGmt(mps.pathway)
    gscs <- BioQC::readGmt(msigdbC2=msigdb.c2.file,
                    msigdbC7=msigdb.c7.file,
                    msigdbHallmark = msigdb.hallmark.file,
                    immunomics = immunomics.file,
                    immunespace = immunespace.file)
    gscs <- appendGmtList(mpsGSCs, gscs)
  } else {
    gscs <- BioQC::readGmt(MolecularPhenotyping=mps.pathway,
                    upstream=upstream.file,
                    ronet=ronet.file,
                    goslim=goslim.file,
                    mbdisease=mbdisease.file,
                    mbmetabolic=mbmetabolic.file,
                    mbpath=mbpath.file,
                    mbprocess=mbprocess.file,
                    mbtoxicity=mbtoxicity.file,		
                    mbpathology=mbpathology.file,
                    msigdbC2=msigdb.c2.file,
                    msigdbC7=msigdb.c7.file,
                    msigdbHallmark = msigdb.hallmark.file,
                    immunomics = immunomics.file,
                    immunespace = immunespace.file)
  }
  return(gscs)
}
