library(ribiosGSEA)
library(testthat)

expect_identical(ribiosGSEA::parseContributingGenes("RXRA(-1.42),CD36(-1.21),PPARG(-0.56),NR1H3(-0.39),CEBPA(-0.28)")[[1]],
                 data.frame(Gene=c("RXRA", "CD36", "PPARG", "NR1H3", "CEBPA"),
                            Stat=c(-1.42, -1.21, -0.56, -0.39, -0.28)))

test_that("parseGenesetsContributingGenes handles empty input", {
  res <- parseGenesetsContributingGenes(character(0), character(0))
  expect_is(res, "data.frame")
  expect_equal(nrow(res), 0)
  expect_identical(colnames(res), c("GeneSet", "Gene", "Stat"))
})

test_that("parseGenesetsContributingGenes returns data.frame for normal input", {
  res <- parseGenesetsContributingGenes("AKR1C4(-1.25), AKR1D1(-1.11)", "Metabolism")
  expect_is(res, "data.frame")
  expect_equal(nrow(res), 2)
  expect_identical(colnames(res), c("GeneSet", "Gene", "Stat"))
  expect_equal(as.character(res$GeneSet), c("Metabolism", "Metabolism"))
})

test_that("parseCameraContributingGenes handles empty genesets", {
  mockTbl <- data.frame(GeneSet=character(0),
                        ContributingGenes=character(0),
                        stringsAsFactors=FALSE)
  res <- parseCameraContributingGenes(mockTbl, character(0))
  expect_is(res, "list")
  expect_equal(length(res), 0)
})

test_that("parseCameraContributingGenes handles single geneset", {
  mockTbl <- data.frame(GeneSet="PathA",
                        ContributingGenes="AKT1(1.24),AKT2(1.11)",
                        stringsAsFactors=FALSE)
  res <- parseCameraContributingGenes(mockTbl, "PathA")
  expect_is(res, "list")
  expect_equal(length(res), 1)
  expect_true("PathA" %in% names(res))
  expect_true(all(c("AKT1", "AKT2") %in% res[["PathA"]]))
})

test_that("parseCameraContributingGenes handles missing genesets gracefully", {
  mockTbl <- data.frame(GeneSet="PathA",
                        ContributingGenes="AKT1(1.24)",
                        stringsAsFactors=FALSE)
  res <- parseCameraContributingGenes(mockTbl, c("PathA", "PathB"))
  expect_equal(length(res), 2)
  expect_true(!is.null(res[["PathA"]]))
})

