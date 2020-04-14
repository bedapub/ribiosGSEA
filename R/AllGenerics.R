#' @include ribiosGSEA-package.R
NULL

#' Return gene-set namespace
#' @param object An object
#' @export
setGeneric("gsNamespace", function(object, ...) standardGeneric("gsNamespace"))

#' Return gene-set name
#' @param object An object
#' @export
setGeneric("gsName", function(object, ...) standardGeneric("gsName"))

#' Return gene-set description
#' @param object An object
#' @export
setGeneric("gsDesc", function(object, ...) standardGeneric("gsDesc"))

#' Return gene-set genes
#' @param object An object
#' @export
setGeneric("gsGenes", function(object,...) standardGeneric("gsGenes"))

#' Return gene-set gene count
#' @param object An object
#' @export
setGeneric("gsGeneCount", function(object,...) standardGeneric("gsGeneCount"))

#' Return gene-set gene values
#' @param object An object
#' @export
setGeneric("gsGeneValues", function(object) standardGeneric("gsGeneValues"))

#' Set gene-set genes
#' @param object An object
#' @param value Value
#' @export
setGeneric("gsGenes<-", function(object,value) standardGeneric("gsGenes<-"))

#' Set gene-set gene statistics (values)
#' @param object An object
#' @param value Value
#' @export
setGeneric("gsGeneValues<-", function(object,value) standardGeneric("gsGeneValues<-"))

#' Returns a logical value whether the enrichment belongs to the core
#' @param object An object
#' @export
setGeneric("isGseaCoreEnrich", function(object) standardGeneric("isGseaCoreEnrich"))

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

#' Return GSEA core enrichment genes
#' @param object An object
#' @export
setGeneric("gseaCoreEnrichGenes", function(object) standardGeneric("gseaCoreEnrichGenes"))

#' Return annotated BROAD GSEA result item
#' @param object An object
#' @export
setGeneric("annoBroadGseaResItem", function(object, ...) standardGeneric("annoBroadGseaResItem"))

#' Return annotated BROAD GSEA result
#' @param object An object
#' @export
setGeneric("annoBroadGseaRes", function(object) standardGeneric("annoBroadGseaRes"))

##----------------------------------------##
## Fisher's exact test
##----------------------------------------##

#' Return hits
#' @param object An object
#' @export
setGeneric("hits", function(object, ...) standardGeneric("hits"))

#' Return P-values
#' @param object An object
#' @export
setGeneric("pValue", function(object, ...) standardGeneric("pValue"))

#' Return FDR values
#' @param object An object
#' @export
setGeneric("fdrValue", function(object, ...) standardGeneric("fdrValue"))

#' Return gene-set size
#' @param object An object
#' @export
setGeneric("gsSize", function(object, ...) standardGeneric("gsSize"))

#' Return gene-set effect size
#' @param object An object
#' @export
setGeneric("gsEffSize", function(object, ...) standardGeneric("gsEffSize"))

#' Return the minimum FDR value
#' @param object An object
#' @export
setGeneric("minFdrValue", function(object, ...) standardGeneric("minFdrValue"))

#' Return the minimum P value
#' @param object An object
#' @export
setGeneric("minPValue", function(object, ...) standardGeneric("minPValue"))

#' Return a logical vector indicating whether a gene-set is significantly
#' enriched
#' @param object An object
#' @param fdr Numeric, FDR threshold
#' @export
setGeneric("isSigGeneSet", function(object, fdr,...) standardGeneric("isSigGeneSet"))

#' Return a character vector indicating significantly enriched gene-sets
#' @param object An object
#' @param fdr Numeric, FDR threshold
#' @export
setGeneric("sigGeneSet", function(object, fdr,...) standardGeneric("sigGeneSet"))

#' Return a table of significantly enriched gene-ests
#' @param object An object
#' @param fdr Numeric, FDR threshold
#' @export
setGeneric("sigGeneSetTable", function(object, fdr,...) standardGeneric("sigGeneSetTable"))

#' Return top enriched gene-sets
#' @param object An object
#' @param N integer, number of top gene-sets to be returned
#' @export
setGeneric("topGeneSetTable", function(object, N,...) standardGeneric("topGeneSetTable"))

#' Return top or significantly enriched gene-sets
#' @param object An object
#' @param N integer, number of top gene-sets to be returned
#' @param fdr Numeric, FDR threshold
#' @export
setGeneric("topOrSigGeneSetTable", function(object, N, fdr, ...) standardGeneric("topOrSigGeneSetTable"))

#' Perform Fisher's exact test
#' @param genes Genes
#' @param genesets Gene-sets
#' @param universe The universe of genes
#' @export
setGeneric("fisherTest", function(genes, genesets, universe,  ...) standardGeneric("fisherTest"))

#' Filter by size
#' @param object An object
#' @param min Integer, minimum size
#' @param max Integer, maximum size
#' @export
setGeneric("filterBySize", function(object,min,max) standardGeneric("filterBySize"))
