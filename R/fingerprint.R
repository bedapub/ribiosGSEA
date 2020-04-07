## read and transform Q-values from GSEA output files
gseaResQvalue <- function(file, threshold=1E-4, log=FALSE, posLog=FALSE) {
  tbl <- read.table(file, sep="\t", header=TRUE)
  name <- tbl[,"NAME"]
  q <- tbl[,"FDR.q.val"]
  if(!is.null(threshold) && !is.na(threshold))
    q[q<threshold] <- threshold
  if(log)
    q <- log10(q)
  if(posLog)
    q <- (-q)
  return(data.frame(name=name, value=q))
}
## read enrichment scores from GSEA output files


#' Read GSEA statistic for pathway fingerprinting
#' 
#' Read GSEA statistics (log-transformed q-value [q], Enrichment Score [ES], or
#' normalized Enrichement Score [NES]) to profile pathway activitities.
#' 
#' In many cases we want to extract pathway signatures from a set of
#' experiments. Both \code{gseaResQvalue} and \code{gseaES} can read GSEA
#' output files and extract desired statistic: q-value, ES or NES.
#' 
#' See the GSEA document for definitions of the three values. For comparing a
#' few conditions to another, we recommend using \emph{q-value}. For
#' large-scale comparisons between pathways (or other gene signatures), we have
#' found \emph{ES} very useful. It is adviced to choose proper statistic to
#' extract pathway signatures only when you are sure of the aim. Using any
#' statistic without good reasoning may as always lead to wrong intepretations
#' of the data.
#' 
#' These functions are usually not directly called by end-users. See
#' \code{\link{gseaFingerprint}} and \code{link{gseaFingerprintMatrix}}
#' instead.
#' 
#' @aliases gseaResQvalue gseaResES
#' @param file GSEA output tab-delimited file, usually with the file name
#' \sQuote{gsea\_report\_for.*\_pos\_.*\.xls} or
#' \sQuote{gsea\_report\_for.*\_neg\_.*\.xls}. Located in GSEA output
#' directory.
#' @param threshold Valid for q value: what is the minimum threshold of q-value
#' (FDR)? It can be set to the number of permutation tests divided by \code{1}.
#' By default \code{1/10000}
#' @param log Valid for q value: whether the FDR q value should be transformed
#' by base-10 (log10) logarithm. By default \code{FALSE}
#' @param posLog Valid for q value: whether the logged FDR q value should be
#' negated to get positive value.This is useful when the sign of \code{q} is
#' used to distinguish between positive and negative enriched pathways. By
#' default \code{FALSE}.
#' @param normalized Valid for enrichment score: if set to \code{TRUE},
#' normalized enrichment score (nes) will be returned instead of (es). By
#' default set to \code{FALSE}
#' @return A \code{data.frame} with two columns: \code{name} and \code{value}.
#' The column \code{name} contains gene signatures (e.g. pathways), and
#' \code{value} contains the statistic.
#' @author Jitao David Zhang <jitao_david.zhang@@roche.com>, with input from
#' Martin Ebeling, Laura Badi and Isabelle Wells.
#' @seealso End-users will probably find \code{\link{gseaFingerprint}} and
#' \code{link{gseaFingerprintMatrix}} more useful, since they operate on the
#' level of GSEA result directories, instead of single output tab-delimited
#' files.
#' @references GSEA documentation
#' \url{http://www.broadinstitute.org/gsea/doc/GSEAUserGuideFrame.html}
#' @examples
#' 
#' gsea.dir <- system.file(package="ribiosGSEA",
#' "extdata/gseaDirs/VitaminA_24h_High")
#' gsea.file <- file.path(gsea.dir,
#' "gsea_report_for_na_neg_1336489010730.xls")
#' 
#' gsea.q <- gseaResQvalue(gsea.file)
#' gsea.logq <- gseaResQvalue(gsea.file, log=TRUE)
#' gsea.logq.pos <- gseaResQvalue(gsea.file, log=TRUE, posLog=TRUE)
#' 
#' gsea.es <- gseaResES(gsea.file)
#' gsea.nes <- gseaResES(gsea.file, normalized=TRUE)
#' 
#' @export gseaResES
gseaResES <- function(file, normalized=FALSE) {
  tbl <- read.table(file, sep="\t", header=TRUE)
  name <- tbl[,1L]
  vc <- ifelse(normalized, "NES", "ES")
  v <- tbl[,vc]
  return(data.frame(name=name, value=v))
}

## read finger print from GSEA directory


#' Extract pathway fingerprints from GSEA results
#' 
#' \code{gseaFingerprint} extracts pathway fingerprints from the result of one
#' GSEA result. \code{gseaFingerprintMatrix} extracts multiple signatures and
#' organizes into the form of rectangular matrix.
#' 
#' \code{gseaFingerprint} extracts pathway signature from one GSEA output
#' directory. While \code{gseaFingerprintMatrix} simultaneously extracts from
#' more than one GSEA output directories, and organizes pathway signatures in a
#' rectangular matrix form.
#' 
#' \code{gseaFingerprintMatrix} takes care of signature mapping between
#' different GSEA result sets.
#' 
#' @aliases gseaFingerprint gseaFingerprintMatrix
#' @param gseaDir Character, a GSEA output directory. Notice the directory must
#' be accessible by the R session. A common mistake is to use a relative path
#' which cannot be found.
#' @param gseaDirs Character vector, GSEA output directories
#' @param value Character, the statistic to extract, currently supporting
#' \code{q}, \code{es} and \code{nes}
#' @param threshold Numeric, minimum threshold of q-value, passed to
#' \code{gseaQvalue}
#' @param sortByName Logical, whether signatures should be sorted by name
#' @param ... Parameters passed to \code{gseaFingerprint} by
#' \code{gseaFingerprintMatrix}
#' @return \code{gseaFingerprint} returns a \code{data.frame} with two columns
#' \code{name} and \code{value}, recording gene signature (pathway) names and
#' the statistic chosen by the user.
#' 
#' \code{gseaFingerprintMatrix} returns a \code{matrix}, with the union set of
#' gene signatures from all GSEA output result sets as rows, and GSEA result
#' names as columns.
#' @author Jitao David Zhang <jitao_david.zhang@@roche.com>
#' @seealso See \code{gseaQvalue} and \code{gseaES} for how to choose the
#' statistic to produce pathway signatures.
#' @examples
#' 
#' gsea.dir <- system.file(package="ribiosGSEA","extdata/gseaDirs/")
#' gsea.dirs <- dir(gsea.dir, full.names=TRUE)
#' gsea.fp <- gseaFingerprint(gsea.dirs[1], value="q")
#' gsea.fps <- gseaFingerprintMatrix(gsea.dirs, value="q")
#' 
#' @export gseaFingerprint
gseaFingerprint <- function(gseaDir, value=c("q", "es", "nes"), threshold=1E-4, sortByName=TRUE) {
  value <- match.arg(value)
  xls <- dir(gseaDir, pattern="gsea_report_for_.*\\.xls")
  pos.xls <- grep("gsea_report_for_.*_pos_.*\\.xls", xls)
  neg.xls <- grep("gsea_report_for_.*_neg_.*\\.xls", xls)

  if(length(xls)==2 & length(pos.xls)==0 & length(neg.xls)==0) {
    keywords <- gsub("gsea_report_for_(.*)_[0-9]*\\.xls", "\\1", xls)
    kw <- tolower(keywords)
    pos.xls <- order(kw)[1L]
    neg.xls <- order(kw)[2L]
  }
  
  if(length(pos.xls)==1L) {
    if(value=="q") {
      poss <- gseaResQvalue(file.path(gseaDir, xls[pos.xls]), threshold=threshold, log=TRUE, posLog=TRUE)
    } else {
      nes <- ifelse(value=="nes", TRUE, FALSE)
      poss <- gseaResES(file.path(gseaDir, xls[pos.xls]), normalized=nes)
    }
  } else {
    poss <- NULL
  }
  if(length(neg.xls)==1L) {
    if(value=="q") {
      negs <- gseaResQvalue(file.path(gseaDir, xls[neg.xls]), threshold=threshold, log=TRUE, posLog=FALSE)
    } else {
      nes <- ifelse(value=="nes", TRUE, FALSE)
      negs <- gseaResES(file.path(gseaDir, xls[neg.xls]), normalized=nes)
    }
  } else {
    negs <- NULL
  }
  paths <- rbind(poss, negs)
  if(!is.null(paths) && sortByName) {
    paths <- sortByCol(paths, "name", decreasing=FALSE)
  }
  return(paths)
}

gseaFingerprintMatrix <- function(gseaDirs, value="es",...) {
  hs.fps <- lapply(gseaDirs,gseaFingerprint, value=value,...)
  isNull <- sapply(hs.fps, is.null)
  if(all(isNull))
    stop("No valid GSEA output directories were detected.\n",
         "Please make sure that input directories are GSEA result folders (not their parent folders)\n")
  fps <- hs.fps[!isNull]
  fps.names <- unique(unlist(lapply(fps, function(x) x[,1L])))
  
  gseaNames <- gsub("\\.GseaPreranked\\.[0-9]*$", "", basename(gseaDirs))
  gseaNames <- gseaNames[!isNull]
  
  fpsMat <- matrix(NA, nrow=length(fps.names), ncol=length(fps),
                   dimnames=list(fps.names, gseaNames))
  for(i in seq(along=fps)) {
    fpsMat[match(fps[[i]][,1L], fps.names),i] <- fps[[i]][,2L]
  }
  return(fpsMat)
}
