library(ribiosGSEA)
library(testthat)

context("Test strWithNumber")

days <- c("D1", "D10", "D15", "D3.5")

test_that("orderByNumberInStr works", {
    expect_identical(orderByNumberInStr(days),
		     c(3, 2, 4, 1))
    expect_identical(orderByNumberInStr(days, decreasing=FALSE),
		     c(1, 4, 2, 3))
})

test_that("factorByNumberInStr works", {
    expect_identical(factorByNumberInStr(days, decreasing=FALSE),
		     factor(days, levels=c("D1", "D3.5", "D10", "D15")))
    expect_identical(orderByNumberInStr(days, decreasing=TRUE),
		     factor(days, levels=c("D15", "D10", "D3.5", "D1")))
})
