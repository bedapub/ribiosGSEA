#' Pretty RONET Gene-set Names
#' 
#' @param x Character strings, RONET gene-set names
#' @param nchar Integer, number of chararacters to be displayed. 
#' 
#' @return Character strings
#' 
#' @examples
#' strs <- c("ARNT_GeneID405_negativeTargets",
#'   "Neurophysiological_process_nNOS_signaling_in_neuronal_synapses",
#'   "NR5A1_GeneID2516_allTargets", 
#'   "IL4_GeneID3565_negativeTargets",
#'   "Apoptosis_REACTOME")
#' prettyRonetGenesetNames(strs)
#'   
prettyRonetGenesetNames <- function(x, nchar=50) {
  res <- gsub("_REACTOME", "", x)
  res <- gsub("MB_[A-Za-z]*_", "", res)
  res <- gsub("GeneID([0-9]*)", "(GeneID=\\1)", res)
  res <- gsub("allTargets", "targets", res)
  res <- gsub("positiveTargets", "positive targets", res)
  res <- gsub("negativeTargets", "negative targets", res)
  res <- gsub("\\(.*\\)$", "", res)
  res <- gsub("_", " ", res)
  res <- gsub(" Downstream$", "", res)
  res <- gsub("\\s+", " ", res)
  res <- ribiosUtils::shortenStr(res, nchar=nchar)
  return(res)
}
