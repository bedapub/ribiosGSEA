#' An adapted and enhanced version of limma::camera
#' 
#' @param y 	a numeric matrix of log-expression values or log-ratios of
#' expression values, or any data object containing such a matrix. Rows
#' correspond to probes and columns to samples. Any type of object that can be
#' processed by getEAWP is acceptable.  
#' @param index an index vector or a list
#' of index vectors. Can be any vector such that y[index,] of statistic[index]
#' selects the rows corresponding to the test set. The list can be made using
#' \code{ids2indices}.  
#' @param design Design matrix 
#' @param contrast contrast of
#' the linear model coefficients for which the test is required. Can be an
#' integer specifying a column of design, or else a numeric vector of same
#' length as the number of columns of design.  
#' @param weights numeric matrix of
#' observation weights of same size as \code{y}, or a numeric vector of array
#' weights with length equal to \code{ncol(y)}, or a numeric vector of gene
#' weights with length equal to \code{nrow(y)}.  
#' @param geneLabels Labels of the features in the input matrix.
#' @param use.ranks do a rank-based test (TRUE) or
#' a parametric test (FALSE)?  
#' @param allow.neg.cor should reduced variance
#' inflation factors be allowed for negative correlations?  
#' @param trend.var
#' logical, should an empirical Bayes trend be estimated? See \code{eBayes} for
#' details.  
#' @param sort logical, should the results be sorted by p-value?
#' @param .fixed.inter.gene.cor Numeric value, vector, or \code{NULL}/\code{NA},
#' advanced parameter corresponding to \code{inter.gene.cor} in the original
#' implementation in limma. If set, gene-sets are set to have the fixed
#' inter-gene correlation; the vector will be recycled to meet the correct
#' length. If set as \code{NULL}/\code{NA}, correlations are estimated from each
#' gene-set.  
#' @param .approx.zscoreT logical, advanced parameter only used for
#' debugging purposes. If \code{TRUE}, the code is expected to return the exact
#' same results as edgeR::camera (version 3.20.9), and maybe faster in
#' execution.
#' 
#' The function was adapted from \code{\link[limma]{camera}}, with following
#' improvments \enumerate{ \item The output data.frame is more user-friendly
#' \item The column 'FDR' is always present, even when only one gene-set was
#' tested \item Scores are calculated, defined as
#' \code{log10(pValue)*I(directionality)}, where \code{I(directionality)} equals
#' \code{1} if the directionality is \code{Up} and \code{-1} if the
#' directionality is \code{Down} \item Contributing genes and statistics are
#' printed }
#' 
#' @return A \code{data.frame} with one row per set and the following columns:
#' \describe{ 
#' \item{GeneSet}{Gene set name} 
#' \item{NGenes}{Number of genes in the set}
#' \item{Correlation}{Estimated correlation} 
#' \item{EffectSize}{Estimated
#' difference between the mean values of genes in the geneset and the background
#' genes} 
#' \item{Direction}{Direction of set-wise regulation, \code{Up} or
#' \code{Down}}
#' \item{Score}{Gene-set enrichment score, defined as
#' \code{log10(pValue)*I(directionality)}, where \code{I(directionality)} equals
#' \code{1} if the directionality is \code{Up} and \code{-1} if the
#' directionality is \code{Down}}
#' \item{ContribuingGenes}{A character string,
#' containing all genes labels of genes that are in the set and regulated in the
#' same direction as the set-wise direction, and the respective statistic} }
#' 
#' @note Since limma 3.29.6, the default setting of allow.neg.cor changes from
#' TRUE to FALSE, and a new parameter, inter.gene.cor, is added with the default
#' value of 0.01, namely a prior inter-gene correlation is set for all gene
#' sets. Currently, \code{biosCamera} does not have the parameter
#' \code{inter.gene.cor}, but \code{allow.neg.cor} is set by default to
#' \code{FALSE} to be consistent with the latest camera function.
#' 
#' @importFrom ribiosUtils haltifnot
#' @importFrom limma squeezeVar zscoreT rankSumTestWithCorrelation
#'
#' @examples 
#' y <- matrix(rnorm(1000*6),1000,6)
#' design <- cbind(Intercept=1,Group=c(0,0,0,1,1,1))
#' # First set of 20 genes are genuinely differentially expressed 
#' index1 <- 1:20
#' y[index1,4:6] <- y[index1,4:6]+1
#' # The second set of 20 genes are not
#' index2 <- 21:40
#' biosCamera(y, index1, design) 
#' biosCamera(y, index2, design)
#' biosCamera(y, list(index1, index2), design)
#'
#' # compare with the output of camera: columns 'GeneSet', 'Score',
#' # 'ContributingGenes' are missing, and in case \code{inter.gene.cor} is (as
#' # default) set to a numeric value, the column 'Correlation' is also missing
#'
#' limmaDefOut <- limma::camera(y, index1, design)
#' limmaCorDefOut <-
#'     limma::camera(y, index1, design, inter.gene.cor=NA)
#'
#' \dontrun{ 
#'   # when \code{.approx.zscoreT=TRUE},  PValue reported by
#'   # \code{limma::camera(inter.gene.cor=NA)} and \code{ribiosGSEA::biosCamera}
#'   # should equal 
#'   biosCorOut <- biosCamera(y, index1, design, .approx.zscoreT=TRUE)
#'
#'   # when \code{.fixed.inter.gene.cor=0.01} and \code{.approx.zscoreT=TRUE},
#'   # PValue reported by \code{limma::camera} and \code{ribiosGSEA::biosCamera}
#'   # should equal 
#'   biosFixCorOut <- biosCamera(y, index1, design,
#'       .fixed.inter.gene.cor=0.01, .approx.zscoreT=TRUE)
#'   testthat::expect_equal(biosFixCorOut$PValue, limmaDefOut$PValue)
#' }
#'
#' @export
biosCamera <- function (y, index, design = NULL, contrast = ncol(design), 
			weights = NULL,
                        geneLabels=NULL,
                        use.ranks = FALSE, allow.neg.cor = FALSE, 
			trend.var = FALSE, 
                        sort = FALSE,
                        .fixed.inter.gene.cor = NULL,
                        .approx.zscoreT=FALSE) {
    y <- as.matrix(y)
    G <- nrow(y)
    n <- ncol(y)
    if(is.null(geneLabels)) {
        geneLabels <- rownames(y)
        if(is.null(geneLabels))
            geneLabels <- 1:nrow(y)
    } else {
        haltifnot(length(geneLabels)==nrow(y),
                  msg="geneLabels's length must equal to nrow(y)")
    }
    if(G < 3)
      stop("Too few genes in the dataset: need at least 3")
    if(!is.list(index)) 
        index <- list(set1 = index)
    if(is.null(names(index)))
      names(index) <- sprintf("set%d", seq(along=index))
    nsets <- length(index)
    if (nsets == 0L)
      stop("Geneset (index) is empty")
    if (is.null(design)) {
      stop("no design matrix")
    } else {
      design <- as.matrix(design)
      if (mode(design) != "numeric")
        stop("design must be a numeric matrix")
      if(nrow(design) != n)
        stop(paste0("row dimensions of the design matrix must ",
		    "match the column dimension of data"))
    }
    p <- ncol(design)
    df.residual <- n - p
    if (df.residual < 1)
      stop("No residual df: cannot compute t-tests")
    fixed.cor <- length(.fixed.inter.gene.cor)>1 || 
	   !(is.null(.fixed.inter.gene.cor) || 
	     is.na(.fixed.inter.gene.cor))
    if(fixed.cor) {
      df.camera <- ifelse(use.ranks, Inf, G-2)
      .fixed.inter.gene.cor <- rep_len(.fixed.inter.gene.cor, nsets)
    } else {
      df.camera <- min(df.residual, G - 2)
    }
    if (!is.null(weights)) {
        if (any(weights <= 0)) 
            stop("weights must be positive")
        if (length(weights) == n) {
            sw <- sqrt(weights)
            y <- t(t(y) * sw)
            design <- design * sw
            weights <- NULL
        }
    }
    if (!is.null(weights)) {
        if (length(weights) == G) 
            weights <- matrix(weights, G, n)
        weights <- as.matrix(weights)
        if (any(dim(weights) != dim(y))) 
            stop("weights not conformal with y")
    }
    if (is.character(contrast)) {
        contrast <- which(contrast == colnames(design))
        if (length(contrast) == 0) 
            stop("coef ", contrast, " not found")
    }
    if (length(contrast) == 1) {
        if(contrast < p) {
          ## Reorders the to-be-tested contrast to the last 
  	  ## column of the design matrix
            j <- c((1:p)[-contrast], contrast)
            design <- design[, j]
        }
    }  else {
      ## Transforms the design matrix into a new one by t(t(Q) %*% y), 
      ## where Q=qr(contrast) and y=t(design), and then reorders to-be-tested 
      ## contrast to the last column of the design matrix.
      ## The QR decomposition of the contrast matrix is used to 
      ## re-parameterize the design matrix so as to encode the desired 
      ## comparison directly in the last column.
        QR <- qr(contrast)
        design <- t(qr.qty(QR, t(design)))
        if (sign(QR$qr[1, 1] < 0)) 
            design[, 1] <- -design[, 1]
        design <- design[, c(2:p, 1)]
    }
    ## The transformed design matrix will next be used to estimate the effect 
    ## of the contrast, which involves another QR transformation
    if (is.null(weights)) {
        QR <- qr(design)
        if (QR$rank < p) 
            stop("design matrix is not of full rank")
        effects <- qr.qty(QR, t(y))
        unscaledt <- effects[p, ]
        if (QR$qr[p, p] < 0) 
            unscaledt <- -unscaledt
    }  else {
        effects <- matrix(0, n, G)
        unscaledt <- rep(0, G)
        sw <- sqrt(weights)
        yw <- y * sw
        for (g in 1:G) {
            xw <- design * sw[g, ]
            QR <- qr(xw)
            if (QR$rank < p) 
                stop("weighted design matrix not of full rank for gene ", 
                  g)
            effects[, g] <- qr.qty(QR, yw[g, ])
            unscaledt[g] <- effects[p, g]
            if (QR$qr[p, p] < 0) 
                unscaledt[g] <- -unscaledt[g]
        }
    }

    ## Effects is a n x G matrix (n=ncol(y), G=nrow(y)), and U are the 
    ## residuals removing the main effects, transposed, and normalised 
    ## by sqrt(mean(u^2)).
    U <- effects[-(1:p), , drop = FALSE]
    sigma2 <- colMeans(U^2)
    U <- t(U)/sqrt(sigma2)
    if (trend.var) 
        A <- rowMeans(y)
    else A <- NULL
    sv <- squeezeVar(sigma2, df = df.residual, covariate = A)
    modt <- unscaledt/sqrt(sv$var.post)
    if (use.ranks) {
      Stat <- modt
    } else {
      ## moderated-t is further transformed into z-score
      df.total <- min(df.residual + sv$df.prior, G * df.residual)
      Stat <- zscoreT(modt, df = df.total, approx=.approx.zscoreT)
    }
    ## meanStat and varStat are mean and variance of z-score 
    ## of the moderated t statistic of all features in the matrix
    meanStat <- mean(Stat)
    varStat <- var(Stat)
    tab <- matrix(0, nsets, 6)
    rownames(tab) <- NULL
    colnames(tab) <- c("NGenes", "Correlation", "Down", "Up", 
                       "EffectSize", "TwoSided")

    conts <- vector("character", nsets)
    for (i in 1:nsets) {
        iset <- index[[i]]
        ## TODO: also export unscaledt values (which should correspond to logFCs)
        StatInSet <- Stat[iset]
        m <- length(StatInSet)
        m2 <- G - m
        if (fixed.cor) {
          correlation <- .fixed.inter.gene.cor[i]
          vif <- 1 + (m - 1)*correlation
        } else {
          if (m > 1) {
              Uset <- U[iset, , drop = FALSE]
              vif <- m * mean(colMeans(Uset)^2)
              correlation <- (vif - 1)/(m - 1)
          } else {
            vif <- 1
            correlation <- NA
          }
        }
        tab[i, 1] <- m
        tab[i, 2] <- correlation
        effectSize <- mean(StatInSet) - meanStat
        if (use.ranks) {
            if (!allow.neg.cor) 
                correlation <- max(0, correlation)
            tab[i, 3:4] <- limma::rankSumTestWithCorrelation(iset, 
		statistics = Stat, 
                correlation = correlation, df = df.camera)
        } else {
            if (!allow.neg.cor) 
                vif <- max(1, vif)
            delta <- G/m2 * effectSize
            varStatPooled <- ((G - 1) * varStat - delta^2 * m * 
                m2/G)/(G - 2)
            two.sample.t <- delta/sqrt(varStatPooled * (vif/m + 
                1/m2))
            tab[i, 3] <- pt(two.sample.t, df = df.camera)
            tab[i, 4] <- pt(two.sample.t, df = df.camera, lower.tail = FALSE)
        }
        tab[i,5] <- effectSize
        isDown <- tab[i,3] <= tab[i,4]
        if(!is.null(isDown) && !is.na(isDown)) {
            if(isDown) { ## pDown < pUp
                contInds <- iset[StatInSet<meanStat]
            } else {
                contInds <- iset[StatInSet>meanStat]
            }
            contVals <- Stat[contInds]
            contOrd <- order(contVals, decreasing=!isDown)
            contInds <- contInds[contOrd]
            contVals <- contVals[contOrd]
            
            conts[i] <- paste(sprintf("%s(%1.2f)",
                                      geneLabels[contInds], contVals), collapse=",")
        }
    }
    
    tab[, 6] <- 2 * pmin(tab[, 3], tab[, 4])
    tab <- data.frame(tab, stringsAsFactors = FALSE)
    Direction <- rep.int("Up", nsets)
    Direction[tab$Down < tab$Up] <- "Down"
    tab$Direction <- Direction
    tab$PValue <- tab$TwoSided
    tab$Down <- tab$Up <- tab$TwoSided <- NULL
    if (nsets > 1) {
        tab$FDR <- p.adjust(tab$PValue, method = "BH")
    } else {
        tab$FDR <- tab$PValue
    }
    tab$Score <- ribiosUtils::pScore(tab$PValue,
                                     sign=Direction=="Up",
                                     method="absLog10")
    tab$ContributingGenes <- as.character(conts)
    tab$GeneSet <- names(index)
    tab <- putColsFirst(tab, "GeneSet")
    if (sort && nsets > 1) {
        o <- order(tab$PValue)
        tab <- tab[o, ]
    }
    return(tab)
}

#' Apply biosCamera to an expression matrix and GmtList
#'
#' @param matrix Expression matrix, features (genes) in rows and 
#'    samples in columns
#' @param geneSymbols Gene-symbols corresponding to the rows of matrix
#' @param gmtList A \code{\link[BioQC]{GmtList}} object 
#' @param design Design matrix
#' @param contrasts Contrast matrix
#' 
#' @return A \code{data.frame} containing CAMERA results
#' 
#' @importFrom ribiosExpression designMatrix
#' @examples 
#' mat <- matrix(c(10,3,5,9,3,5,
#'                 2,4,8,12,9,9,
#'                 3,5,7,5,4,4,
#'                 3,3,12,12,0,1), ncol=6, byrow=TRUE)
#' rownames(mat) <- sprintf("gene%d", 1:nrow(mat))
#' designMatrix <- matrix(c(rep(1,6), c(0,0,1,1,0,0), c(0,0,0,0,1,1)),
#'                        byrow=FALSE, ncol=3)
#' contrastMatrix <- matrix(c(0,1,0,0,0,1), ncol=2, byrow=FALSE)
#' gs1 <- list(name="GeneSet1", desc="", genes=c("gene1", "gene3"), namespace="default")
#' gs2 <- list(name="GeneSet2", desc="", genes=c("gene2", "gene4"), namespace="default")
#' gs3 <- list(name="GeneSet3", desc="", genes=c("gene1", "gene4"), namespace="default")
#' gmtlist <- BioQC::GmtList(list(gs1, gs2, gs3))
#' cameraOut <- gmtListCamera(mat, rownames(mat), gmtlist, designMatrix, contrastMatrix)
#' cameraOut
#' @importFrom parallel detectCores mclapply
#' @importFrom BioQC gsGenes gsName 
#' @importFrom ribiosUtils sortByCol 
#' @export
gmtListCamera <- function(matrix, geneSymbols, gmtList, design, contrasts) {
  genes <- BioQC::gsGenes(gmtList)
  genes.inds <- lapply(genes, function(x) {
    ind <- match(x, geneSymbols)
    return(ind[!is.na(ind)])
  })
  names(genes.inds) <- BioQC::gsName(gmtList)

  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    cl <- 2L
  } else {
    # use all cores 
    cl <- parallel::detectCores()
  }
  cl <- pmin(cl, ncol(contrasts))
  cameraRes <- parallel::mclapply(1:ncol(contrasts),
                        function(x) {
                          tbl <- biosCamera(matrix,
                                            design=design,
                                            index=genes.inds,
                                            contrast=contrasts[,x],
                                            geneLabels=geneSymbols,
                                            sort=FALSE)
                          if(!"FDR" %in% colnames(tbl)) {
                            ## TRUE if there is only one gene set
                            tbl$FDR <- tbl$PValue
                          }
                          tbl <- tbl[,c("GeneSet", "NGenes", "Correlation", "Direction", "EffectSize", 
                                        "PValue", "FDR", "ContributingGenes")]
                          tbl <- sortByCol(tbl, "PValue",decreasing=FALSE)
                          return(tbl)
                        }, mc.cores=cl)
  
  cRes <- do.call(rbind, cameraRes)
  
  if(is.null(colnames(contrasts)))
    colnames(contrasts) <- sprintf("Contrast%d", 1:ncol(contrasts))
  
  bg <- data.frame(Contrast=rep(colnames(contrasts), sapply(cameraRes, nrow)))
  res <- cbind(bg, cRes)
  rownames(res) <- NULL
  res <- subset(res, NGenes>=1 & !is.na(PValue) & !is.nan(PValue))
  return(res)
}
