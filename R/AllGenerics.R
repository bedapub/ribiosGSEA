#' @include ribiosGSEA-package.R
NULL

#' Return gene-set namespace
#' @param object An object
#' @param ... Other parameters
#' @return A character vector of gene-set namespaces.
#' @export
setGeneric("gsNamespace", function(object, ...) standardGeneric("gsNamespace"))

#' Return gene-set name
#' @param object An object
#' @param ... Other parameters
#' @return A character vector of gene-set names.
#' @export
setGeneric("gsName", function(object, ...) standardGeneric("gsName"))

#' Return gene-set genes
#' @param object An object
#' @param ... Other parameters
#' @return A character vector or list of character vectors of gene-set genes.
#' @export
setGeneric("gsGenes", function(object,...) standardGeneric("gsGenes"))

#' Return gene-set gene count
#' @param object An object
#' @param ... Other parameters
#' @return An integer vector of gene counts.
#' @export
setGeneric("gsGeneCount", function(object,...) standardGeneric("gsGeneCount"))

#' Return gene-set gene values
#' @param object An object
#' @return A numeric vector or list of numeric vectors of gene values.
#' @export
setGeneric("gsGeneValues", function(object) standardGeneric("gsGeneValues"))

#' Set gene-set genes
#' @title gsGenes-set
#' @param object An object
#' @param value Value
#' @return The modified object.
#' @export
setGeneric("gsGenes<-", function(object,value) standardGeneric("gsGenes<-"))

#' Set gene-set gene statistics (values)
#' @title gsGeneValues-set
#' @param object An object
#' @param value Value
#' @return The modified object.
#' @export
setGeneric("gsGeneValues<-", function(object,value) standardGeneric("gsGeneValues<-"))

#' Return GSEA enrichment scores
#' @param object An object
#' @return A numeric vector of enrichment scores.
#' @export
setGeneric("gseaES", function(object) standardGeneric("gseaES"))

#' Return GSEA normalized enrichment scores
#' @param object An object
#' @return A numeric vector of normalized enrichment scores.
#' @export
setGeneric("gseaNES", function(object) standardGeneric("gseaNES"))

#' Return GSEA number of permutation
#' @param object An object
#' @return A numeric vector.
#' @export
setGeneric("gseaNP", function(object) standardGeneric("gseaNP"))

#' Return GSEA FDR
#' @param object An object
#' @return A numeric vector of FDR values.
#' @export
setGeneric("gseaFDR", function(object) standardGeneric("gseaFDR"))

#' Return GSEA FWER values
#' @param object An object
#' @return A numeric vector of FWER values.
#' @export
setGeneric("gseaFWER", function(object) standardGeneric("gseaFWER"))

#' Return gene-set gene indices
#' @param object An object
#' @return An integer vector of gene indices.
#' @export
setGeneric("gsGeneIndices", function(object) standardGeneric("gsGeneIndices"))

#' Return GSEA enrichment score profile
#' @param object An object
#' @return A numeric vector of enrichment score profiles.
#' @export
setGeneric("gseaESprofile", function(object) standardGeneric("gseaESprofile"))

#' Return GSEA core enrichment score threshold
#' @param object An object
#' @return A numeric value.
#' @export
setGeneric("gseaCoreEnrichThr", function(object) standardGeneric("gseaCoreEnrichThr"))

#' Return GSEA core enrichment genes (also known as leading-edge genes)
#' @aliases gseaLeadingEdgeGenes
#' @param object An object
#' @return A character vector of core enrichment genes.
#' @export
setGeneric("gseaCoreEnrichGenes", function(object) standardGeneric("gseaCoreEnrichGenes"))

##----------------------------------------##
## Fisher's exact test
##----------------------------------------##

#' Return hits
#' @param object An object
#' @param ... Other parameters
#' @return A character vector or list of hit genes.
#' @export
setGeneric("hits", function(object, ...) standardGeneric("hits"))

#' Return P-values
#' @param object An object
#' @param ... Other parameters
#' @return A numeric vector of p-values.
#' @export
setGeneric("pValue", function(object, ...) standardGeneric("pValue"))

#' Return FDR values
#' @param object An object
#' @param ... Other parameters
#' @return A numeric vector of FDR values.
#' @export
setGeneric("fdrValue", function(object, ...) standardGeneric("fdrValue"))

#' Return the effective size of gene-set
#' @param object An object
#' @param ... Other parameters
#' @return An integer vector of effective sizes.
#' @export
setGeneric("gsEffectiveSize", function(object, ...) standardGeneric("gsEffectiveSize"))

#' Perform Fisher's exact test
#' @param genes Genes
#' @param genesets Gene-sets
#' @param universe The universe of genes
#' @param ... Other parameters
#' @return A \code{FisherResult} object or a \code{data.table} of results.
#' @export
setGeneric("fisherTest", function(genes, genesets, universe,  ...) standardGeneric("fisherTest"))

#' Filter by size
#' @param object An object
#' @param min Integer, minimum size
#' @param max Integer, maximum size
#' @return The filtered object.
#' @export
setGeneric("filterBySize", function(object,min,max) standardGeneric("filterBySize"))
