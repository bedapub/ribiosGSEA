## test fisherTest and S4 methods
library(testthat)
library(ribiosUtils)
library(ribiosGSEA)
library(BioQC)

##----------------------------------------##
## test fisherTest, the underlying function
##----------------------------------------##

inputGenes <- LETTERS[1:3]
## geneSet1: complete inclusion
geneSet1 <- LETTERS[1:6]
## geneSet1: complete exclusion
geneSet2 <- LETTERS[4:7]
## geneSet1: partial overlap
geneSet3 <- LETTERS[2:5]
## geneSet1: genes in the gene set outside of the universe
geneSet4 <- c(LETTERS[2:5], letters[1:3])
universe <- LETTERS

badGenes <- function(vector, dup=2, naCount=3, sample=FALSE) {
    res <- c(rep(vector, dup), rep(NA, naCount))
    if(sample)
        res <- sample(res)
    return(res)
}

gsEnrich1 <- ribiosGSEA::fisherTest(inputGenes, geneSet1, universe, "testName", "testNamespace")
gsManual1 <- fisher.test(matrix(c(3,3,0,20), nrow=2, byrow=TRUE), alternative="greater")
expect_equal(gsManual1$p.value, ribiosGSEA::pValue(gsEnrich1))
expect_equal(length(geneSet1), ribiosGSEA::gsEffectiveSize(gsEnrich1))
expect_that(inputGenes, equals(ribiosGSEA::hits(gsEnrich1)))
expect_that(gsNamespace(gsEnrich1), equals("testNamespace"))
expect_that(gsName(gsEnrich1), equals("testName"))
            
gsEnrich2 <- ribiosGSEA::fisherTest(inputGenes, geneSet2, universe)
gsManual2 <- fisher.test(matrix(c(0,4,3,19), nrow=2, byrow=TRUE), alternative="greater")
expect_that(gsManual2$p.value, equals(pValue(gsEnrich2)))
expect_that(length(geneSet2), equals(ribiosGSEA::gsEffectiveSize(gsEnrich2)))
expect_that(character(), equals(ribiosGSEA::hits(gsEnrich2)))

gsEnrich3 <- ribiosGSEA::fisherTest(inputGenes, geneSet3, universe)
gsManual3 <- fisher.test(matrix(c(2,2,1,21), nrow=2, byrow=TRUE), alternative="greater")
expect_that(gsManual3$p.value, equals(pValue(gsEnrich3)))
expect_that(length(geneSet3), equals(ribiosGSEA::gsEffectiveSize(gsEnrich3)))
expect_that(c("B", "C"), equals(ribiosGSEA::hits(gsEnrich3)))


gsEnrich4 <- ribiosGSEA::fisherTest(inputGenes, geneSet4, universe)
gsManual4 <- fisher.test(matrix(c(2,2,1,21), nrow=2, byrow=TRUE), alternative="greater")
expect_that(gsManual4$p.value, equals(pValue(gsEnrich4)))
expect_that(4L, equals(ribiosGSEA::gsEffectiveSize(gsEnrich4)))
expect_that(c("B", "C"), equals(ribiosGSEA::hits(gsEnrich4)))


##----------------------------------------##
## test fisherTest for character, character, character
##----------------------------------------##
gsS4.1 <- fisherTest(inputGenes, geneSet1, universe, "testName", "testNamespace")
expect_identical(gsS4.1, gsEnrich1)

gsS4.2 <- fisherTest(inputGenes, geneSet2, universe)
expect_identical(gsS4.2, gsEnrich2)

gsS4.3 <- fisherTest(inputGenes, geneSet3, universe)
expect_identical(gsS4.3, gsEnrich3)

gsS4.4 <- fisherTest(inputGenes, geneSet4, universe)
expect_identical(gsS4.4, gsEnrich4)

#### manipulate input genes
## test that both redundant and NA genes are removed
gsS4.1.bad <- fisherTest(badGenes(inputGenes), geneSet1, universe, "testName", "testNamespace")
expect_identical(gsS4.1.bad, gsS4.1)

#### manipulate gene sets
gsS4.4.bad <- fisherTest(inputGenes, badGenes(geneSet4), universe)
expect_identical(gsS4.4.bad, gsS4.4)

#### manipulate genes, geneset genes, and universe simultaneously
gsS4.4.ubad <- fisherTest(badGenes(inputGenes), badGenes(geneSet4), badGenes(universe))
expect_identical(gsS4.4.ubad, gsS4.4)

##----------------------------------------##
## test fisherTest for character, list, character
##----------------------------------------##
myGeneSet4 <- list(name="Hamburger Sportverein", desc="some description", genes=geneSet4, namespace="Bundesliga")
myGeneSet4.fisher <- fisherTest(badGenes(inputGenes), myGeneSet4, badGenes(universe))
expect_that(pValue(myGeneSet4.fisher), equals(pValue(gsS4.4)))
expect_that(gsName(myGeneSet4.fisher), equals("Hamburger Sportverein"))
expect_that(gsNamespace(myGeneSet4.fisher), equals("Bundesliga"))

##----------------------------------------##
## test fisherTest for character, GmtList, character
##----------------------------------------##
gs1 <- list(name="GeneSet1", desc="", genes=geneSet1, namespace="A")
gs2 <- list(name="GeneSet2", desc="", genes=geneSet2, namespace="A")
gs3 <- list(name="GeneSet3", desc="", genes=geneSet3, namespace="A")
gs4 <- list(name="GeneSet4", desc="", genes=geneSet4, namespace="B")
gss <- BioQC::GmtList(list(gs1, gs2, gs3, gs4))
myFisherRes <- fisherTest(inputGenes, gss, universe)

myFisherRes.expP <- c(pValue(gsEnrich1),pValue(gsEnrich2),pValue(gsEnrich3),pValue(gsEnrich4))
expect_that(myFisherRes$PValue, equals(myFisherRes.expP))
expect_that(myFisherRes$GeneSetNamespace, equals(c(rep("A", 3), "B")))
expect_that(myFisherRes$GeneSetName, equals(sprintf("GeneSet%d", 1:4)))

myFisherRes.bad <- fisherTest(badGenes(inputGenes), gss, badGenes(universe))
expect_identical(myFisherRes.bad ,myFisherRes)
