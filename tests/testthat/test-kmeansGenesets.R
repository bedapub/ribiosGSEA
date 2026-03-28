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

test_that("kmeansGeneset handles single gene-set (nrow=1)", {
  mat1 <- matrix(1:3, nrow=1, dimnames=list("gs1", c("c1", "c2", "c3")))
  genes1 <- list(gs1=c("A", "B", "C"))
  res <- kmeansGeneset(mat1, genes1, optK=5)
  expect_identical(names(res), expKgsNames)
  expect_null(res$kmeans)
  expect_equal(length(res$repGenesets), 1)
  expect_equal(res$repGenesets, "gs1")
  expect_equal(nrow(res$genesetClusterData), 1)
  expect_true(res$genesetClusterData$IsRepresentative[1])
})

test_that("kmeansGeneset handles two gene-sets (nrow=2)", {
  mat2 <- matrix(rnorm(6), nrow=2, dimnames=list(c("gs1", "gs2"),
                                                   c("c1", "c2", "c3")))
  genes2 <- list(gs1=c("A", "B"), gs2=c("B", "C"))
  res <- kmeansGeneset(mat2, genes2, optK=5)
  expect_identical(names(res), expKgsNames)
  expect_true(length(res$repGenesets) >= 1)
  expect_true(length(res$repGenesets) <= 2)
})

test_that("kmeansGeneset handles three gene-sets with optK=1 (single cluster)", {
  mat3 <- matrix(rnorm(9), nrow=3, dimnames=list(c("gs1", "gs2", "gs3"),
                                                   c("c1", "c2", "c3")))
  genes3 <- list(gs1=c("A", "B"), gs2=c("B", "C"), gs3=c("A", "C"))
  res <- kmeansGeneset(mat3, genes3, optK=1)
  expect_identical(names(res), expKgsNames)
  expect_true(length(res$repGenesets) >= 1)
})
