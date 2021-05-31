library(ribiosGSEA)
library(testthat)

set.seed(1887)
profMat <- matrix(rnorm(100), nrow=20,
    dimnames=list(sprintf("geneset%d", 1:20), sprintf("contrast%d", 1:5)))
gsGenes <- lapply(1:nrow(profMat), function(x)
    unique(sample(LETTERS, 10, replace=TRUE)))
names(gsGenes) <- rownames(profMat)
kgsRes <- kmeansGeneset(profMat, gsGenes, optK=5)
kgsSingleRes <- kmeansGeneset(profMat[, 1L, drop=FALSE], gsGenes, optK=5)

expKgsNames <- c("kmeans", "genesetClusterData", 
                 "repGenesets", "repGenesetClusters")

testthat::test_that("kmeansGeneset returns an object as expected",  {
  expect_identical(names(kgsSingleRes), expKgsNames)
  expect_identical(names(kgsRes), expKgsNames)
})
