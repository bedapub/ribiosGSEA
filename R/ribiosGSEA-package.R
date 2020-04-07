

#' Parsing GMT file for downstream analysis
#' 
#' The function reads in a GMT file, map identifiers to the given character
#' vector, and filter by size.
#' 
#' 
#' @param file A GMT file
#' @param vec A character vector of gene identifiers in the same order as in
#' the expression matrix, i.e. the row names or a column in the
#' \code{ExpressionSet} feature data.
#' @param min Integer (optional), filter gene sets with less genes than the
#' threshold
#' @param max Integer (optional), filter gene sets with more genes than the
#' threhsold
#' @return A list of filtered gene set integer indices
#' @author Jitao David Zhang <jitao_david.zhang@@roche.com> %% ~~who you are~~
#' @seealso \code{\link{readGmt}} to read the GMT file only.
#' @examples
#' 
#' sample.gmt.file <- system.file("extdata", "test.gmt",package="ribiosIO")
#' sample.genes <- c("GLB1", "GLRX", "GNAL", "GNG5", "IDH2", "HTR4", "KLF8", "KLK10")
#' (sample.gmtId <- parseGmt(sample.gmt.file, sample.genes))
#' (sample.gmtId <- parseGmt(sample.gmt.file, sample.genes, min=1, max=20))
#' 
NULL





#' Variables useful for GSEA analysis in Roche bioinformatics environment.
#' 
#' These variables ar useful for performing GSEA analysis in Roche
#' bioinformatics environment.
#' 
#' 
#' @name GSEA_DATA_DIR
#' @aliases GSEA_DATA_DIR GSEA_ANNOTATION_DIR GSEA_GENESET_DIR JAVA_BIN
#' GSEA_JAR DEFAULT_GMT DEFAULT_CHIP
#' @docType data
NULL



