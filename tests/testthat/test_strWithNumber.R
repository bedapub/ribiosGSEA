library(ribiosGSEA)
library(testthat)

context("Test strWithNumber")

days <- c("D1", "D10", "D15", "D3.5")

test_that("orderByNumberInStr works", {
    expect_identical(orderByNumberInStr(days, decreasing=TRUE),
		     c(3L, 2L, 4L, 1L))
    expect_identical(orderByNumberInStr(days, decreasing=FALSE),
		     c(1L, 4L, 2L, 3L))
})

test_that("factorByNumberInStr works", {
    expect_identical(factorByNumberInStr(days, decreasing=FALSE),
		     factor(days, levels=c("D1", "D3.5", "D10", "D15")))
    expect_identical(factorByNumberInStr(days, decreasing=TRUE),
		     factor(days, levels=c("D15", "D10", "D3.5", "D1")))
})
