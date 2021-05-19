#' @include biosCamera.R
NULL

#' Wrap the gage::gage method to report consistent results as the CAMERA method
#' 
#' @param logFC A named vector of logFC values of genes
#' @param gmtList A \code{\link[BioQC]{GmtList}} object containing gene-sets
#' @param ... Other parameters passed to \code{\link[gage]{gage}}
#' 
#' @importFrom gage gage
#' @export
myGage <- function(logFC, gmtList, ...) {
  gsnames <- gsName(gmtList)
    genes <- gsGenes(gmtList)
    gsizes <- gsSize(gmtList)
    cate <- gsNamespace(gmtList)
    gdf <- data.frame(geneset=gsnames, namespace=cate)
    gage.res <- gage::gage(logFC, gsets=genes, 
                           set.size=c(min(gsizes), max(gsizes)),
                           ref=NULL, samp=NULL, ...)
    greater <- as.data.frame(gage.res$greater[,c("stat.mean", "p.val", "q.val", "set.size")])
    less <- as.data.frame(gage.res$less[,c("p.val", "q.val")])
    greater$geneset <- gsnames
    less$geneset <- gsnames
    greater <- merge(greater, gdf, by="geneset")
    less <- merge(less, gdf, by="geneset")
    res.raw <- merge(greater, less, by=c("geneset","namespace"),
                     suffix=c(".greater", ".less"))
    direction <- with(res.raw, ifelse(p.val.less<p.val.greater, "Down", "Up"))
    pVal.pmin <- with(res.raw, ifelse(p.val.less<p.val.greater, p.val.less, p.val.greater))
    pVal <- pVal.pmin * 2; pVal[pVal>1] <- 1

    ## TODO: contributing genes

    res <- data.frame(Namespace=res.raw$namespace,
                      GeneSet=res.raw$geneset,
                      NGenes=res.raw$set.size,
                      Direction=direction,
                      PValue=pVal,
                      FDR=rep(NA,length(pVal)))
    res <- subset(res, NGenes>=1 & !is.na(PValue) & !is.nan(PValue))
    res$FDR <- p.adjust(res$PValue, "fdr")
    return(res)
}

#' Perform the GAGE analysis for EdgeResult and GmtList
#'
#' @param edgeResult An \code{EdgeResult} object.
#' @param gmtList A \code{GmtList} object.
#' @return A \code{data.frame} containing enrichment analysis results.
#'
#' @importFrom ribiosUtils putColsFirst
#' @importFrom ribiosNGS dgeTables fData
#' @export
logFCgage <- function(edgeResult, gmtList) {
    geneSymbols <- fData(edgeResult)$GeneSymbol
    logFCs <- lapply(ribiosNGS::dgeTables(edgeResult), function(x) {
                         res <- x$logFC
                         names(res) <- x$GeneSymbol
                         return(res)
                     })
    erTables <- lapply(logFCs, function(x) myGage(x, gmtList))
    
    erTable <- do.call(rbind, erTables)
    erTable$Contrast <- rep(names(logFCs),sapply(erTables, nrow))
    erTable <- putColsFirst(erTable, c("Namespace", "Contrast", "GeneSet"))
    rownames(erTable) <- NULL
    
    return(erTable)
}

##----------------------------------------##
## camera
##----------------------------------------##
#' Calculate mid-p quantile residuals
#' 
#' @param y An DGEList object
#' @param design Design matrix
#' @param contrast Contrast vector
#' 
#' The function is a carbon copy of edgeR:::.zscoreDGE, which is unfortunately
#' not exported
#' @examples
#' 
#' dgeMatrix <- matrix(rpois(1200, 10), nrow=200)
#' dgeList <- DGEList(dgeMatrix)
#' dgeList <- edgeR::estimateCommonDisp(dgeList)
#' dgeDesign <- model.matrix(~gl(2,3))
#' dgeZscore <- zscoreDGE(dgeList, dgeDesign, contrast=c(0,1))
#' head(dgeZscore)
#' 
#' @importFrom edgeR getDispersion zscoreNBinom glmFit
#' @importFrom limma contrastAsCoef
#' @export zscoreDGE
zscoreDGE <- function(y, design=NULL, contrast=ncol(design)) {
  allzero <- rowSums(y$counts > 1e-08) == 0
  if (any(allzero)) 
    warning(sum(allzero), "rows with all zero counts")
  dispersion <- getDispersion(y)
  if (is.null(dispersion)) 
    stop("Dispersion estimate not found. Please estimate the dispersion(s) before you proceed.")
  if (is.null(design))
    design <- y$design
  if (is.null(design)) {
    if (nlevels(y$samples$group) < 2) 
      stop("design not supplied and samples all belong to the same group")
    design <- model.matrix(~y$samples$group)
    rownames(design) <- colnames(y)
  }
  nbeta <- ncol(design)
  if (nbeta < 2) 
    stop("design matrix must have at least two columns")
  if (is.character(contrast)) {
    if (length(contrast) > 1) 
      stop("contrast should specify only one column of design")
    contrast <- which(contrast == colnames(design))
    if (!length(contrast)) 
      stop("contrast doesn't match any column of design")
  }
  if (length(contrast) == 1) {
    design0 <- design[, -contrast, drop = FALSE]
  }
  else {
    design <- contrastAsCoef(design, contrast = contrast, 
                             first = FALSE)$design
    design0 <- design[, -nbeta, drop = FALSE]
  }
  fit.null <- edgeR::glmFit(y, design0, prior.count = 0)
  y <- zscoreNBinom(y$counts, 
                    mu = pmax(fit.null$fitted.values, 1e-17),
                    size = 1/dispersion)
  return(y)
}

#' @importFrom parallel detectCores
getCores <- function(maxCores) {
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    cl <- 2L
  } else {
    # use all cores 
    cl <- parallel::detectCores()
  }
  cl <- pmin(cl, maxCores)
  return(cl)
}

#' Merge CAMERA results using limma default parameters and biosCamera parameters
#' @param matrix A numeric matrix, passed to \code{camera} and \code{biosCamera}
#' @param index An index vector or a list of index vectors of features.
#' @param designMatrix Design matrix.
#' @param contrast A numeric vector of the same length as the number of columns in the design matrix, coefficients of contrasts.
#' @param featureLabels A character vector of the same length as the number of rows of the matrix, feature labels, for instance gene symbols.
#' @param weights NULL or numeric matrix of precision weights, passed to \code{camera}.
#' @param use.ranks Logical, passed to \code{camera}.
#' 
#' The function merges the output of \code{camera} with default options with
#' the output of \code{biosCamera}, which appends additional information to 
#' camera methods, to return a comprehensive table as output.
#' 
#' @seealso \code{\link{camera}}, \code{\link{biosCamera}}
#'
#' @examples 
#' 
#' y <- matrix(rnorm(1000*6),1000,6)
#' features <- sprintf("Feature%d", 1:nrow(y))
#' design <- cbind(Intercept=1,Group=c(0,0,0,1,1,1))
#' # First set of 20 genes are genuinely deferentially expressed 
#' index1 <- 1:20
#' y[index1,4:6] <- y[index1,4:6]+1
#' # The second set of 20 genes are not
#' index2 <- 21:40
#' index1Res <- mergeCameraResults(y, index=index1, 
#'   designMatrix=design, contrast=c(0,1), featureLabels=features)
#' index1ListRes <- mergeCameraResults(y, index=list(index1), 
#'    designMatrix=design, contrast=c(0,1), featureLabels=features)
#' index12ListRes <- mergeCameraResults(y, index=list(index1, index2), 
#'    designMatrix=design, contrast=c(0,1), featureLabels=features)
#' @export
mergeCameraResults <- function(matrix, 
                               index,
                               designMatrix, 
                               contrast,
                               featureLabels,
                               weights=NULL, 
                               use.ranks=FALSE) {
  haltifnot(length(featureLabels)==nrow(matrix),
            msg="featureLabels must a character vector of the same length as the row count of the input matrix")
  haltifnot(length(contrast)==ncol(designMatrix),
            msg="contrast must be a numeric vector of the same length as the column count of the design matrix")
  
  tbl.priorCor <- limma::camera.default(matrix,
                                        index=index,
                                        design=designMatrix,
                                        contrast=contrast,
                                        weights=weights,
                                        use.ranks = use.ranks, 
                                        allow.neg.cor=FALSE,
                                        inter.gene.cor=0.01,
                                        trend.var=FALSE,
                                        sort=FALSE)
  tbl.estCor <- ribiosGSEA::biosCamera(matrix, 
                                       index=index, 
                                       design = designMatrix, 
                                       contrast = contrast,
                                       weights = weights, 
                                       use.ranks = use.ranks, 
                                       geneLabels = featureLabels, 
                                       allow.neg.cor = FALSE, 
                                       trend.var = FALSE, 
                                       sort = FALSE)
  if(!is.list(index) || (is.list(index) && length(index)==1)) {
    tbl.priorCor$FDR <- tbl.priorCor$PValue
    tbl.estCor$FDR <- tbl.estCor$PValue
  }
  tbl <- cbind(tbl.estCor,
               PValue.cor0.01=tbl.priorCor$PValue,
               FDR.cor0.01=tbl.priorCor$FDR,
               Score.cor0.01=ribiosUtils::pScore(tbl.priorCor$PValue,
                                                 tbl.priorCor$Direction=="Up",
                                                 method="absLog10"))
  rownames(tbl) <- NULL
  tbl <- tbl[,c("GeneSet", "NGenes","Direction",
                "Correlation", "EffectSize", "PValue", "FDR", "Score",
                "PValue.cor0.01", "FDR.cor0.01", "Score.cor0.01", 
                "ContributingGenes")]
  tbl <- ribiosUtils::sortByCol(tbl, "PValue",decreasing=FALSE)
  return(tbl)  
}

#' Apply the CAMERA method to a DGEList object and a contrast
#' 
#' @param dgeList A DGEList object, with GeneSymbol available, and dispersion must be estimated
#' @param index List of integer indices of genesets, names are names of gene
#' sets
#' @param design Design matrix
#' @param contrasts Contrast matrix
#' @param doParallel Logical, whether \code{parallel::mclapply} should be used. Since at the current setting it makes a job running forever, use \code{TRUE} only if you are debugging the code.
#' 
#' @importFrom parallel mclapply
#' @importFrom limma camera
#' @export camera
#' @importFrom edgeR estimateDisp
#' @importFrom ribiosNGS humanGeneSymbols
#' @importFrom ribiosUtils sortByCol
#' 
#' @examples 
#' exMat <- matrix(rpois(120, 10), nrow=20, ncol=6)
#' exGroups <- gl(2,3, labels=c("Group1", "Group2"))
#' exDesign <- model.matrix(~0+exGroups)
#' colnames(exDesign) <- levels(exGroups)
#' exContrast <- matrix(c(-1,1), ncol=1, dimnames=list(c("Group1", "Group2"), c("Group2.vs.Group1")))
#' exDescon <- DesignContrast(exDesign, exContrast, groups=exGroups)
#' exFdata <- data.frame(GeneSymbol=sprintf("Gene%d", 1:nrow(exMat)))
#' exPdata <- data.frame(Name=sprintf("Sample%d", 1:ncol(exMat)),
#'                      Group=exGroups)
#' exDgeList <- DGEList(exMat, genes=exFdata, samples=exPdata)
#' exDgeList <- edgeR::estimateDisp(exDgeList, exDesign)
#' cameraDGEListByContrast(exDgeList, index=1:5, design=exDesign, contrasts=exContrast)
#' cameraDGEListByContrast(exDgeList, 
#'   index=list(1:5, 6:10),
#'   design=exDesign, contrasts=exContrast)
#'   
#' @export
cameraDGEListByContrast <- function(dgeList, index, design, contrasts, doParallel=FALSE) {
  geneSymbols <- ribiosNGS::humanGeneSymbols(dgeList)
  if(is.null(geneSymbols))
    stop("EdgeResult must have 'GeneSymbol' in its fData to perform camera!")
  
  if(doParallel) {
    cl <- getCores(ncol(contrasts))
  } else {
    cl <- 1L
  }
  cameraRes <- lapply(1:ncol(contrasts),
                      function(x) {
                        contrast <- contrasts[,x]
                        zscores <- zscoreDGE(y=dgeList,
                                             design=design,
                                             contrast=contrast)
                        res <- mergeCameraResults(matrix=zscores,
                                                  index=index,
                                                  designMatrix=design,
                                                  contrast=contrast,
                                                  featureLabels = geneSymbols,
                                                  weights=NULL, 
                                                  use.ranks=FALSE)
                        return(res)
                      })
  
  cRes <- do.call(rbind, cameraRes)
  
  if(is.null(colnames(contrasts)))
    colnames(contrasts) <- sprintf("Contrast%d", 1:ncol(contrasts))
  
  bg <- data.frame(Contrast=rep(colnames(contrasts), sapply(cameraRes, nrow)))
  res <- cbind(bg, cRes)
  rownames(res) <- NULL
  res <- subset(res, NGenes>=1 & !is.na(PValue) & !is.nan(PValue))
  return(res)
}

#' Run CAMERA method using EdgeResult 
#' 
#' @param y A EdgeResult object
#' @param gmtList Gene set collections, for example read by
#' \code{\link[BioQC]{readGmt}}, with namespace.
#' @param doParallel Logical, whether \code{parallel::mclapply} should be used. Since at the current setting it makes a job running forever, use \code{TRUE} only if you are debugging the code.
#' @param ... Not used
#' 
#' Note that the EdgeResult object must have a column 'GeneSymbol' in its
#' \code{fData}.
#' 
#' @return A \code{data.frame} containing CAMERA results.
#' 
#' @importFrom ribiosExpression contrastNames designMatrix contrastMatrix
#' @importFrom ribiosNGS humanGeneSymbols dgeWithEdgeR
#' @importFrom ribiosUtils putColsFirst
#' @importFrom BioQC gsNamespace matchGenes GmtList
#' @importClassesFrom ribiosNGS EdgeResult
#' @importFrom limma camera
#' @examples 
#' exMat <- matrix(rpois(120, 10), nrow=20, ncol=6)
#' exGroups <- gl(2,3, labels=c("Group1", "Group2"))
#' exDesign <- model.matrix(~0+exGroups)
#' colnames(exDesign) <- levels(exGroups)
#' exContrast <- matrix(c(-1,1,1,-1), ncol=2, dimnames=list(c("Group1", "Group2"), 
#'   c("Group2.vs.Group1", "Group1.vs.Group2")))
#' exDescon <- DesignContrast(exDesign, exContrast, groups=exGroups)
#' exFdata <- data.frame(GeneSymbol=sprintf("GeneSymbol%d", 1:nrow(exMat)))
#' exPdata <- data.frame(Name=sprintf("Sample%d", 1:ncol(exMat)),
#'                      Group=exGroups)
#' exDgeList <- DGEList(exMat, genes=exFdata, samples=exPdata)
#' exDgeList <- edgeR::estimateDisp(exDgeList, exDesign)
#' exEdgeObject <- EdgeObject(exDgeList, exDescon)
#' exEdgeRes <- ribiosNGS::dgeWithEdgeR(exEdgeObject)
#' exGmt <- BioQC::GmtList(list(GeneSet1=sprintf("GeneSymbol%d", 1:5),
#'   GeneSet2=sprintf("GeneSymbol%d", 6:10)))
#'   
#' exCameraRes <- camera(exEdgeRes, exGmt)
#' 
#' @export
camera.EdgeResult <- function(y, gmtList, doParallel=FALSE, ...) {
  ctnames<- contrastNames(y)
  design <- designMatrix(y)
  ct <- contrastMatrix(y)
  geneSymbols <- humanGeneSymbols(y)
  if(is.null(geneSymbols))
    stop("EdgeResult must have 'GeneSymbol' in its fData to perform camera!")

  namespaces <- BioQC::gsNamespace(gmtList)
  if(!is.character(namespaces) || is.factor(namespaces)) {
    namespaces <- rep("default", length(gmtList))
  }
  gsIndex <- BioQC::matchGenes(gmtList, geneSymbols)
  names(gsIndex) <- make.unique(names(gsIndex))

  erTables <- tapply(seq(along=gsIndex), namespaces, function(i) {
    currInd <- gsIndex[i]
    tt <- cameraDGEListByContrast(y@dgeList,
                        index=currInd,
                        design=design, contrasts=ct,
                        doParallel=doParallel)
    return(tt)
  })
  
  erTable <- do.call(rbind, erTables)
  erTable$Namespace <- rep(names(erTables), sapply(erTables, nrow))
  erTable <- putColsFirst(erTable, "Namespace")
  rownames(erTable) <- NULL
  
  return(erTable)
}

#' Apply the CAMERA method to a DGEList object
#' 
#' @param limmaVoomResults A LimmaVoomResults object, with GeneSymbol available
#' @param index List of integer indices of genesets, names are names of gene
#' sets
#' @param doParallel Logical, whether \code{parallel::mclapply} should be used. Since at the current setting it makes a job running forever, use \code{TRUE} only if you are debugging the code.
#' @param ... Not used
#' 
#' @return A \code{data.frame} containing CAMERA results.
#' 
#' @importFrom parallel mclapply
#' @importFrom limma camera
#' @importFrom edgeR estimateDisp
#' @importFrom ribiosNGS humanGeneSymbols dgeWithLimmaVoom
#' @importFrom ribiosUtils sortByCol
#' 
#' @examples 
#' exMat <- matrix(rpois(120, 10), nrow=20, ncol=6)
#' exGroups <- gl(2,3, labels=c("Group1", "Group2"))
#' exDesign <- model.matrix(~0+exGroups)
#' colnames(exDesign) <- levels(exGroups)
#' exContrast <- matrix(c(-1,1), ncol=1, dimnames=list(c("Group1", "Group2"), c("Group2.vs.Group1")))
#' exDescon <- DesignContrast(exDesign, exContrast, groups=exGroups)
#' exFdata <- data.frame(GeneSymbol=sprintf("Gene%d", 1:nrow(exMat)))
#' exPdata <- data.frame(Name=sprintf("Sample%d", 1:ncol(exMat)),
#'                      Group=exGroups)
#' exDgeList <- DGEList(exMat, genes=exFdata, samples=exPdata)
#' exDgeList <- edgeR::estimateDisp(exDgeList, exDesign)
#' edgeObj <- EdgeObject(exDgeList, exDescon)
#' limmaVoomRes <- ribiosNGS::dgeWithLimmaVoom(edgeObj)
#' cameraLimmaVoomResultsByContrast(limmaVoomRes, index=c(1:5))
#' cameraLimmaVoomResultsByContrast(limmaVoomRes, index=list(GS1=1:5, GS2=6:10))
#' 
#' @export
cameraLimmaVoomResultsByContrast <- function(limmaVoomResults, index, doParallel=FALSE, ...) {
  geneSymbols <- ribiosNGS::humanGeneSymbols(limmaVoomResults)
  if(is.null(geneSymbols))
    stop("LimmaVoomResults must have 'GeneSymbol' in its fData to perform camera!")
  
  design <- designMatrix(limmaVoomResults)
  contrasts <- contrastMatrix(limmaVoomResults)
  
  if(doParallel) {
    cl <- getCores(ncol(contrasts))
  } else {
    cl <- 1L
  }
  cameraRes <- lapply(1:ncol(contrasts),
                      function(x) {
                        contrast <- contrasts[,x]
                        res <- mergeCameraResults(matrix=limmaVoomResults@voom,
                                                  index=index,
                                                  designMatrix=design,
                                                  contrast=contrast,
                                                  featureLabels = geneSymbols,
                                                  weights=NULL, 
                                                  use.ranks=FALSE)
                        return(res)
                      })
  
  cRes <- do.call(rbind, cameraRes)
  
  if(is.null(colnames(contrasts)))
    colnames(contrasts) <- sprintf("Contrast%d", 1:ncol(contrasts))
  
  bg <- data.frame(Contrast=rep(colnames(contrasts), sapply(cameraRes, nrow)))
  res <- cbind(bg, cRes)
  rownames(res) <- NULL
  res <- subset(res, NGenes>=1 & !is.na(PValue) & !is.nan(PValue))
  return(res)
}

#' Run the CAMERA method using LimmaVoomResult
#' 
#' @param y A LimmaVoomResult object
#' @param gmtList Gene set collections, for example read by
#' \code{\link[BioQC]{readGmt}}
#' @param doParallel Logical, whether \code{parallel::mclapply} should be used. Since at the current setting it makes a job running forever, use \code{TRUE} only if you are debugging the code.
#' @param ... Passed to \code{cameraLimmaVoomResultsByContrast}
#' 
#' Note that the LimmaVoomResult object must have a column 'GeneSymbol' in its
#' \code{fData}.
#' 
#' @return A \code{data.frame} containing CAMERA results.
#' 
#' @importFrom ribiosExpression contrastNames designMatrix contrastMatrix
#' @importFrom ribiosNGS humanGeneSymbols
#' @importFrom ribiosUtils putColsFirst
#' @importFrom BioQC gsNamespace matchGenes GmtList
#' @importClassesFrom ribiosNGS LimmaVoomResult
#' @examples 
#' exMat <- matrix(rpois(120, 10), nrow=20, ncol=6)
#' exGroups <- gl(2,3, labels=c("Group1", "Group2"))
#' exDesign <- model.matrix(~0+exGroups)
#' colnames(exDesign) <- levels(exGroups)
#' exContrast <- matrix(c(-1,1), ncol=1, dimnames=list(c("Group1", "Group2"), c("Group2.vs.Group1")))
#' exDescon <- DesignContrast(exDesign, exContrast, groups=exGroups)
#' exFdata <- data.frame(GeneSymbol=sprintf("Gene%d", 1:nrow(exMat)))
#' exPdata <- data.frame(Name=sprintf("Sample%d", 1:ncol(exMat)),
#'                      Group=exGroups)
#' exDgeList <- DGEList(exMat, genes=exFdata, samples=exPdata)
#' exDgeList <- edgeR::estimateDisp(exDgeList, exDesign)
#' edgeObj <- EdgeObject(exDgeList, exDescon)
#' limmaVoomRes <- ribiosNGS::dgeWithLimmaVoom(edgeObj)
#' exGmt <- BioQC::GmtList(list(GeneSet1=sprintf("GeneSymbol%d", 1:5),
#'   GeneSet2=sprintf("GeneSymbol%d", 6:10)))
#'   
#' camera(limmaVoomRes, exGmt)
#' 
#' @export
camera.LimmaVoomResult <- function(y, gmtList, doParallel=FALSE, ...) {
  ctnames<- contrastNames(y)
  design <- designMatrix(y)
  ct <- contrastMatrix(y)
  geneSymbols <- humanGeneSymbols(y)
  if(is.null(geneSymbols))
    stop("EdgeResult must have 'GeneSymbol' in its fData to perform camera!")
  
  namespaces <- BioQC::gsNamespace(gmtList)
  if(!is.character(namespaces) || is.factor(namespaces)) {
    namespaces <- rep("default", length(gmtList))
  }
  gsIndex <- BioQC::matchGenes(gmtList, geneSymbols)
  names(gsIndex) <- make.unique(names(gsIndex))
  
  erTables <- tapply(seq(along=gsIndex), namespaces, function(i) {
    currInd <- gsIndex[i]
    tt <- cameraLimmaVoomResultsByContrast(y,
                                           index=currInd,
                                           doParallel=doParallel,
                                           ...)
    return(tt)
  })
  
  erTable <- do.call(rbind, erTables)
  erTable$Namespace <- rep(names(erTables), sapply(erTables, nrow))
  erTable <- putColsFirst(erTable, "Namespace")
  rownames(erTable) <- NULL
  
  return(erTable)
}



