#' @include ribiosGSEA-package.R
NULL

#' Return gene-set namespace
#' @param object An object
#' @param ... Other parameters
#' @export
setGeneric("gsNamespace", function(object, ...) standardGeneric("gsNamespace"))

#' Return gene-set name
#' @param object An object
#' @param ... Other parameters
#' @export
setGeneric("gsName", function(object, ...) standardGeneric("gsName"))

#' Return gene-set genes
#' @param object An object
#' @param ... Other parameters
#' @export
setGeneric("gsGenes", function(object,...) standardGeneric("gsGenes"))

#' Return gene-set gene count
#' @param object An object
#' @param ... Other parameters
#' @export
setGeneric("gsGeneCount", function(object,...) standardGeneric("gsGeneCount"))

#' Return gene-set gene values
#' @param object An object
#' @export
setGeneric("gsGeneValues", function(object) standardGeneric("gsGeneValues"))

#' Set gene-set genes
#' @title gsGenes-set
#' @param object An object
#' @param value Value
#' @export
setGeneric("gsGenes<-", function(object,value) standardGeneric("gsGenes<-"))

#' Set gene-set gene statistics (values)
#' @title gsGeneValues-set
#' @param object An object
#' @param value Value
#' @export
setGeneric("gsGeneValues<-", function(object,value) standardGeneric("gsGeneValues<-"))

#' Return GSEA enrichment scores
#' @param object An object
#' @export
setGeneric("gseaES", function(object) standardGeneric("gseaES"))

#' Return GSEA normalized enrichment scores
#' @param object An object
#' @export
setGeneric("gseaNES", function(object) standardGeneric("gseaNES"))

#' Return GSEA number of permutation
#' @param object An object
#' @export
setGeneric("gseaNP", function(object) standardGeneric("gseaNP"))

#' Return GSEA FDR
#' @param object An object
#' @export
setGeneric("gseaFDR", function(object) standardGeneric("gseaFDR"))

#' Return GSEA FWER values
#' @param object An object
#' @export
setGeneric("gseaFWER", function(object) standardGeneric("gseaFWER"))

#' Return gene-set gene indices
#' @param object An object
#' @export
setGeneric("gsGeneIndices", function(object) standardGeneric("gsGeneIndices"))

#' Return GSEA enrichment score profile
#' @param object An object
#' @export
setGeneric("gseaESprofile", function(object) standardGeneric("gseaESprofile"))

#' Return GSEA core enrichment score threshold
#' @param object An object
#' @export
setGeneric("gseaCoreEnrichThr", function(object) standardGeneric("gseaCoreEnrichThr"))

#' Return GSEA core enrichment genes (also known as leading-edge genes)
#' @aliases gseaLeadingEdgeGenes
#' @param object An object
#' @export
setGeneric("gseaCoreEnrichGenes", function(object) standardGeneric("gseaCoreEnrichGenes"))

##----------------------------------------##
## Fisher's exact test
##----------------------------------------##

#' Return hits
#' @param object An object
#' @param ... Other parameters
#' @export
setGeneric("hits", function(object, ...) standardGeneric("hits"))

#' Return P-values
#' @param object An object
#' @param ... Other parameters
#' @export
setGeneric("pValue", function(object, ...) standardGeneric("pValue"))

#' Return FDR values
#' @param object An object
#' @param ... Other parameters
#' @export
setGeneric("fdrValue", function(object, ...) standardGeneric("fdrValue"))

#' Return the effective size of gene-set(s)
#' @param object An object
#' @param ... Other parameters
#' @export
setGeneric("gsEffSize", function(object, ...) standardGeneric("gsEffSize"))

#' Perform Fisher's exact test
#' @param genes Genes
#' @param genesets Gene-sets
#' @param universe The universe of genes
#' @param ... Other parameters
#' @export
setGeneric("fisherTest", function(genes, genesets, universe,  ...) standardGeneric("fisherTest"))

#' Filter by size
#' @param object An object
#' @param min Integer, minimum size
#' @param max Integer, maximum size
#' @export
setGeneric("filterBySize", function(object,min,max) standardGeneric("filterBySize"))
