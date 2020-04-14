#' Extract gene-set namespace from RONET GMT files
#' @param GmtList A GmtList object read from a RONET GMT file
#' @return Character vector of the same length, indicating categories
#' @export
ronetGeneSetNamespace <- function(gmtList) {
  sapply(strsplit(sapply(gmtList, function(x) x$desc), "\\|"), "[[", 1L)
}

#' Read RONET GMT files with namespace information
#' @param file A GMT file in the RONET format, where in the 'desc' field a namespace is appended at the beginning, separated from the rest of the description with a pipe
#' @return A \code{GmtList} object with an additional 'namespace' item in each list
#' @export
readRonetGmt <- function(file) {
  resRaw <- BioQC::readGmt(file)
  res <- new("GmtList")
  res@.Data <- lapply(resRaw, function(x) {
    splitDesc <- strsplit(x$desc, "\\|")[[1]]
    x$namespace <- splitDesc[[1]]
    x$desc <- paste(splitDesc[-1], collapse="|")
    return(x)
  })
  names(res) <- names(resRaw)
  return(res)
}
