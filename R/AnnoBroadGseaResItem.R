#' Convert a BroadGseaResItem object to an AnnoBroadGseaResItem object
#' @param object A BroadGseaResItem object
#' @param genes A character string vector
#' @param geneValues A numeric vector
#' @return An annoBroadGseaResItem object
#' @export
AnnoBroadGseaResItem <- function(object, genes, geneValues) {
  res <- as(object, "AnnoBroadGseaResItem")
  if(!(length(genes)==length(geneValues) && length(genes)==length(gsGeneIndices(object))))
    stop(sprintf("genes (#=%d) and geneValues (#=%d) must be of the same length as the gsGeneIndices (#=%d)",
                 length(genes), length(geneValues), length(gsGeneIndices(object))))
  gsGenes(res) <- genes
  gsGeneValues(res) <- geneValues
  return(res)
}
