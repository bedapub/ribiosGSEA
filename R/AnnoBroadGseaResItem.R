#' Convert a broadGseaResItem object to an annoBroadGseaResItem object
#' @param object A broadGseaResItem object
#' @param genes A character string vector
#' @param geneValues A numeric vector
#' @return An annobroadGseaResItem object
#' @export
annoBroadGseaResItem <- function(object, genes, geneValues) {
  res <- as(object, "annoBroadGseaResItem")
  if(!(length(genes)==length(geneValues) && length(genes)==length(gsGeneIndices(object))))
    stop(sprintf("genes (#=%d) and geneValues (#=%d) must be of the same length as the gsGeneIndices (#=%d)",
                 length(genes), length(geneValues), length(gsGeneIndices(object))))
  gsGenes(res) <- genes
  gsGeneValues(res) <- geneValues
  return(res)
}
