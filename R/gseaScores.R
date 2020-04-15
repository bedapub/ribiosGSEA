#' @include AllMethods.R
NULL

#' Extract scores from GSEA results
#' 
#' One way to score GSEA results is to multiple the absolute value of log10
#' transformed p-values (nominal p-value, FDR, or FWER) with the sign of the
#' enrichment scores. This score is intuitive since it combines statistical
#' significance and the sign of regulation.
#' 
#' \code{gseaScores} takes care of the situation where some gene sets are
#' missing in one or more conditions.
#' 
#' @aliases gseaScore gseaScores
#' @param x An \code{AnnoBroadGseaRes} object
#' @param type Character string, the type of p-value used to calculate the
#' score.
#' @return \code{gseaScore} returns a double vector of scores with gene set
#' names.
#' 
#' \code{gseaScores} returns a data frame of scores, with gene set names as row
#' names.
#' @author Jitao David Zhang <jitao_david.zhang@@roche.com>
#' @seealso \code{\link{gseaNP}}, \code{\link{gseaFDR}}, \code{\link{gseaFWER}}
#' to get p-values.
#' 
#' \code{\link{parseGSEAdir}}.
#' @export gseaScore
gseaScore <- function(x, type=c("fdr", "p", "fwer")) {
  type <- match.arg(type)
  if(type=="fdr") {
    val <- gseaFDR(x)
  } else if (type=="p") {
    val <- gseaNP(x)
  } else if (type=="fwer") {
    val <- gseaFWER(x)
  }
  val[val==0] <- min(val[val!=0], na.rm=TRUE)
  res <- -log10(val) * sign(gseaES(x))
  return(res)
}

#' @describeIn gseaScore gseaScore applied to multiple objects
#' @param ... Objects of \code{AnnoBroadGseaRes} to be compared
#' @param names Character strings, names given to the result score sets. See
#' examples below.
#' @export
gseaScores <- function(..., names=NULL, type=c("fdr", "p", "fwer")) {
  ll <- list(...)
  scores <- lapply(ll, gseaScore, type=type)
  setnames <- munion(lapply(scores, names))
  res <- as.data.frame(sapply(scores, function(x) x[match(setnames, names(x))]))
  rownames(res) <- setnames
  if(!is.null(names)) {
    haltifnot(length(names)==length(ll), msg="names length must match the input list")
    colnames(res) <- names
  }
  return(res)
}

