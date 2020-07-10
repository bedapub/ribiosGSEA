library(ribiosGSEA)
library(profvis)
library(ggplot2)

  exMat <- matrix(rpois(120000, 10), nrow=20000, ncol=12)
  exGroups <- gl(4,3, labels=c("Group1", "Group2", "Group3", "Group4"))
  exDesign <- model.matrix(~0+exGroups)
  exContrast <- matrix(c(-1,1,0,0, 0,0,-1,1),
     ncol=2, byrow=FALSE,
     dimnames=list(c("Group1", "Group2", "Group3", "Group4"),
       c("Group2.vs.Group1", "Group4.vs.Group3")))
  exDescon <- DesignContrast(exDesign, exContrast, groups=exGroups)
  exFdata <- data.frame(GeneSymbol=sprintf("Gene%d", 1:nrow(exMat)))
  exPdata <- data.frame(Name=sprintf("Sample%d", 1:ncol(exMat)),
                       Group=exGroups)
  exObj <- EdgeObject(exMat, exDescon,
                       fData=exFdata, pData=exPdata)
  exDgeRes <- ribiosNGS::dgeWithEdgeR(exObj)

  ngeneset <- 100
  genesetSizes <- round(runif(ngeneset)*100)+1
  exGeneSets <- BioQC::GmtList(lapply(seq(1:ngeneset), function(i) {
    name <- paste0("GeneSet", i)
    desc <- paste0("GeneSet", i)
    genes <- sample(exFdata$GeneSymbol, genesetSizes[i])
    res <- list(name=name, desc=desc, genes=genes, namespace="default")
  }))
  debug(doGse)
  profvis({exGse <- doGse(exDgeRes, exGeneSets, doParallel=FALSE)})
  ## profvis({exGse <- doGse(exDgeRes, exGeneSets)})
  
