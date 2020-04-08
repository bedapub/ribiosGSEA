## getRcrit and getSubGraphs were copied and adapted from the scorem package
getRcrit <- function (a, n) {
  t <- 0 - qt(a, n - 2)
  sqrt(t^2/(t^2 + n - 2))
}

#' @importFrom RGBL connComp
getSubGraphs <- function (object, alpha, nc, w.crit) {
  r.crit <- getRcrit(alpha, nc)
  cc <- connComp(as(new("graphAM", object > r.crit), "graphNEL"))
  sg <- NULL
  for (x in cc) {
    mx <- object[x, x]
    w <- mean(mx)
    n <- length(x)
    if (n == 2L) {
      if (w < w.crit) {
        x <- as.list(x)
      }
    }
    else if (n == 3L) {
      if (w < w.crit) {
        rs <- c(mx[1L,2L], mx[1L,3L], mx[2L,3L])
        if (max(rs) > r.crit) {
          pairs <- list(c(1L, 2L), c(1L, 3L), c(2L, 3L))
          tog <- pairs[[which.max(rs)]]
          sep <- setdiff(1:3, tog)
          x <- list(x[tog], x[sep])
        }
        else {
          x <- as.list(x)
        }
      }
    }
    else {
      if (w < w.crit) {
        x <- getSubGraphs(mx, alpha, nc - 2, w.crit)
      }
    }
    if (!is.list(x)) {
      x <- list(x)
    }
    sg <- c(sg, x)
  }
  sg
}



#' Use Kendall's W and graph theory to assign independent measurements into
#' sub-groups by correlation
#' 
#' Kendall's W, also known as Kendall's coefficient of concordance, is a
#' non-parametric statistic developed to assess agreement among raters used in
#' psychological or similar experimental settings.
#' 
#' In computational biology, the concept of associating features with similar
#' patterns while keeping outliers can be useful in many cases. See the Details
#' section for examples.
#' 
#' This function implements the Kendall's W recursively with graph theory. It
#' split grouped measurements into strongly associated sub-groups. See the
#' Details section.
#' 
#' We take a microarray experiment as an example to demonstrate how the
#' function works. In microarrays, a gene is often represented by more than one
#' probeset, and it is not rare that they do not all resemble the same
#' expression pattern. Usually a \code{one gene-one value} relation is desired.
#' Common practices including choosing the probeset with the highest average
#' signal or the highest variance, as well as taking the mean/median value of
#' all probesets mapped to one gene as the representative value.
#' 
#' Kendall's W takes a very different approach. First it tries to judge whether
#' multiple probesets of one gene are concordant. The concordance is determined
#' by a non-parametric statistic closely related to Spearman correlation
#' coefficient as well as Friedman's test. If all probesets are concordant, it
#' means that their expression patterns are closely associated with each other.
#' Any one of them, or the mean value, can be then used to represent the
#' expression level of the gene.
#' 
#' In cases where there is little concordance among probesets, we can take use
#' of graph theory to iteratively search for sub-groups of probesets resemble
#' each other's expression patterns. In the extreme case, each probeset can be
#' different from the rest, and in this case the number of sub-groups will be
#' equal to the number of probesets mapped to the gene. Such cases can appear,
#' for instance, when each probeset was designed to target a different region
#' of a transcript with splice variants. By using Kendall's W statistic with
#' graph theory, the \code{kendallWmat} function can detect sub-groups with
#' strongly correlated expression patterns, while keeping outliers on their
#' own, therefore providing help for both conventional expression analysis and
#' post-hoc analysis with the help of sequence analysis. See reference for
#' examples on this application.
#' 
#' We believe this approach is only useful for microarray, but can be also
#' interesting for other applications like next-generation sequencing (NGS) or
#' pathway/network analysis. For instance, in NGS experiments, this method can
#' help to determine which splice variants of a transcript have similar
#' expression patterns, and how different are other variants. In pathway
#' analysis, when rows indicate gene expression values and \code{row.factor}
#' indicate pathway membership, the result reveals which sub-networks are
#' regulated associatively.
#' 
#' @param mat A numeric matrix. It must contain at least 2 rows and 2 columns.
#' @param row.factor A factor indicating groups of rows. In expression
#' analysis, for instance, this can be GeneIDs indicating which probesets in
#' rows belong to the same gene.
#' @param summary Character, action to take once the sub-groups have been
#' determined. \sQuote{none} indicates no action should be taken, the original
#' data is returned with the information of sub-grouping. The option
#' \sQuote{mean} (or \sQuote{median}) will take mean/median of features in each
#' sub-group as result. On contrast, \code{max.mean.sig} or \code{max.var.sig}
#' picks the feature of the largest mean signal or the largest variance in each
#' sub-group as the representative. See details below.
#' @param na.rm Logical, should those features whose \code{row.factor} are
#' \code{NA} be left out? If set to \code{TRUE} (which is default), these
#' unannotated features will be discarded from the results.
#' @param alpha Nunmeric value, the significance level of the Kendall's W
#' statistic. The larger the value, the more abbreviations from strong
#' associations are allowed in sub-groups. Default is \code{0.01}.
#' @return Currently a matrix with one attribute slot named \code{info}.
#' @author Jitao David Zhang <jitao_david.zhang@@roche.com>
#' @references The concept of Kendall's W was introduced in the seminal paper
#' \emph{The problem of m rankings} by M.G. Kendall and B.B. Smith (The Annals
#' of Mathematical Statistics, 1939). Schneider, Smith and Hansen developed the
#' SCOREM algorithm combining this statistic with graph theory (\emph{SCOREM:
#' statistical consolidation of redundant expression measures}, Nucleic Acids
#' Research, 2011). This implementation is very much based on the SCOREM
#' algorithm. The main changes are (1) the current implementation is more
#' generic, applicable to native R data structures, therefore able to be
#' applied in other scenario than microarray analysis (2) it takes
#' not-annotated features into account as well and (3) it is possible to
#' directly calculate summary statistics from sub-groups.
#' @examples
#' 
#' ## use a mock example
#' emat <- matrix(c(2,3,5,
#'                  8,9,2,
#'                  3,4,7,
#'                  0,2,1,
#'                  NA, 3, 1.2,
#'                  5, -3,4,
#'                  5,7,11), ncol=3, byrow=TRUE,
#'                dimnames=list(paste("row", 1:7, sep=""),NULL))
#' efac <- factor(c("a", "b", "c", NA, "b", "a", "a"),
#'                levels=letters[1:5])
#' 
#' print(emat)
#' kendallWmat(emat, efac, summary="none")
#' kendallWmat(emat, efac, summary="none", na.rm=FALSE)
#' kendallWmat(emat, efac, summary="mean")
#' kendallWmat(emat, efac, summary="mean", na.rm=FALSE)
#' kendallWmat(emat, efac, summary="median")
#' kendallWmat(emat, efac, summary="median", na.rm=FALSE)
#' kendallWmat(emat, efac, summary="max.mean.sig")
#' kendallWmat(emat, efac, summary="max.mean.sig", na.rm=FALSE)
#' kendallWmat(emat, efac, summary="max.var.sig")
#' kendallWmat(emat, efac, summary="max.var.sig", na.rm=TRUE)
#' 
#' ## kendallW acts as an interface to matrix
#' kendallW(emat, efac, summary="none")
#' 
#' ## kendallW acts as an interface to ExpressionSet
#' data(ribios.ExpressionSet)
#' kendallW(ribios.ExpressionSet, 
#'   Biobase::fData(ribios.ExpressionSet)$GeneID,
#'   summary="none")
#' kendallW(ribios.ExpressionSet, 
#'   Biobase::fData(ribios.ExpressionSet)$GeneID, 
#'   summary="mean")
#' 
#' @export kendallWmat
kendallWmat <- function(mat,
                        row.factor,
                        summary=c("none", "mean", "median", "max.mean.sig", "max.var.sig"),
                        na.rm=TRUE,
                        alpha=0.01) {
  stopifnot(length(row.factor)==nrow(mat))
  if(!is.factor(row.factor))
    row.factor <- factor(row.factor)
  if(is.null(rownames(mat)) || any(duplicated(rownames(mat))))
    stop("'mat' must have unique rownames\n")
  if(missing(summary))
    summary <- "none"
  summary <- match.arg(summary)
  
  rl <- levels(row.factor)
  rt <- table(row.factor)
  
  nc <- ncol(mat)
  r.crit <- getRcrit(alpha, nc)
  w.crit <- (1L + r.crit)/2L
  
  is.mulFac <- row.factor %in% rl[rt > 1L]
  is.sglFac <- row.factor %in% rl[rt==1L]
  is.naFac <- !(is.mulFac | is.sglFac)
  sglmat <- subset(mat, is.sglFac)
  mulmat <- subset(mat, is.mulFac)
  namat <- subset(mat, is.naFac)
  mulmatFac <- factor(subset(row.factor, is.mulFac))
  
  mul2sgl <- by(mulmat, mulmatFac, function(x) {
    mcor <- cor(t(x), method="spearman")
    w <- mean(mcor, na.rm=TRUE)
    if(w >= w.crit) {
      ids <- list(rownames(x))
    } else {
      ids <- getSubGraphs(mcor, alpha, nc, w.crit)
    }
    if(summary=="none") {
      mm <- x
    } else if (summary=="mean") {
      mm <- t(sapply(ids, function(id) colMeans(x[id,,drop=FALSE], na.rm=TRUE)))
      rownames(mm) <- sapply(ids, "[[", 1)
    } else if (summary=="median") {
      mm <- t(sapply(ids, function(id) apply(x[id,,drop=FALSE],2, median, na.rm=TRUE)))
      rownames(mm) <- sapply(ids, "[[", 1)
    } else if (summary=="max.mean.sig") {
      mmp <- sapply(ids, function(id) id[which.max(rowMeans(x[id, ,drop=FALSE], na.rm=TRUE))])
      mm <- x[mmp,]
      rownames(mm) <- mmp
    } else if (summary=="max.var.sig") {
      mmp <- sapply(ids, function(id) {
        if(length(id)==1)
          return(id)
        id[which.max(apply(x[id,,drop=FALSE], 1, sd, na.rm=TRUE))]
      })
      mm <- x[mmp,]
      rownames(mm) <- mmp
    } else {
      stop("Not implemented summary method\n")
    }
    mm <- data.matrix(mm)
    return(list(matrix=mm, ids=ids))
  })
  
  sgls.names <- paste(row.factor[is.sglFac],
                      rownames(sglmat), sep="|")
  if(summary!="none") {
    muls <- unlist(lapply(mul2sgl, function(x) x$ids), use.names=F, recursive=F)
    muls.names <- paste(rep(levels(mulmatFac), sapply(mul2sgl, function(x) length(x$ids))),
                        sapply(muls, paste, collapse="|"),sep="|")
  } else {
    fp <- rep(levels(mulmatFac), sapply(mul2sgl, function(x) length(unlist(x$ids))))
    fs <- lapply(mul2sgl, function(x) lapply(x$ids,
                                             function(id) rep(paste(id, collapse="|"), length(id))))
    fsl <- unlist(fs, use.names=FALSE, recursive=TRUE)
    muls.names <- paste(fp, fsl, sep="|")
  }
  mulmats <- do.call(rbind, lapply(mul2sgl, function(x) x$matrix))
  stopifnot(length(muls.names) == nrow(mulmats))
  
  res.mat <- rbind(sglmat, mulmats)
  nrf <- factor(c(sgls.names, muls.names))
  nrf.o <- sapply(as.character(nrf),
                  function(x) strsplit(x, "\\|")[[1]][[1]])
  orf <- factor(nrf.o,
                levels=levels(row.factor))
  if(!na.rm) {
    res.mat <- rbind(res.mat, namat)
    ana <- rep(NA, sum(is.naFac))
    orf <- factor(c(as.character(orf),
                    as.character(ana)), levels=levels(orf))
    nrf <- factor(c(as.character(nrf),
                    as.character(ana)), levels=levels(nrf))
    res.df <- data.frame(id=rownames(res.mat),
                         row.factor=orf,
                         new.row.factor=nrf)
  } else {
    res.df <- data.frame(id=rownames(res.mat),
                         row.factor=orf,
                         new.row.factor=nrf)
    
  }

  attr(res.mat, "info") <- res.df
  return(res.mat)
}

#' @exportMethod kendallW
setGeneric("kendallW",
           function(object,...) standardGeneric("kendallW"))

#' @exportMethod kendallWinfo
setGeneric("kendallWinfo",
           function(object) standardGeneric("kendallWinfo"))
#' @exportMethod `kendallWinfo<-`
setGeneric("kendallWinfo<-",
           function(object, value) standardGeneric("kendallWinfo<-"))

#' @export 
setMethod("kendallWinfo", "matrix", function(object) {attr(object, "info")})

#' @export 
setReplaceMethod("kendallWinfo", c("matrix","ANY"), function(object, value) {attr(object, "info") <- value; return(object)})

#' @export 
setMethod("kendallW", "matrix", function(object,
                                         row.factor,
                                         summary=c("none", "mean", "median",
                                           "max.mean.sig", "max.var.sig"),
                                         na.rm=TRUE, alpha=0.01) {
  kendallWmat(mat=object, row.factor=row.factor,
              summary=summary, na.rm=na.rm, alpha=alpha)
})

#' @export 
setMethod("kendallW", "ExpressionSet", function(object,
                                                row.factor,
                                                summary=c("none", "mean", "median",
                                                  "max.mean.sig", "max.var.sig"),
                                                na.rm=TRUE, alpha=0.01) {
  exp <- exprs(object)
  new.exp <- kendallWmat(mat=exp, row.factor=row.factor,
                         summary=summary, na.rm=na.rm, alpha=alpha)
  new.info <- kendallWinfo(new.exp)
  kendallWinfo(new.exp) <- NULL
  new.fData <- cbind(fData(object)[match(rownames(new.exp),
                                         featureNames(object)),,drop=FALSE],
                     new.info)
  res <- new("ExpressionSet",
             exprs=new.exp,
             phenoData=phenoData(object),
             featureData=new("AnnotatedDataFrame", new.fData))
  return(res)
})

