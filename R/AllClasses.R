setClass("gseaResItem",
         representation=list(geneset="character",
           "es"="numeric",
           "nes"="numeric",
           "np"="numeric",
           "fdr"="numeric",
           "fwer"="numeric",
           "geneIndices"="integer",
           "esProfile"="numeric",
           "coreEnrichThr"="numeric"))

setClass("annoGseaResItem",
         representation=list("gsGenes"="character",
           "gsGeneValues"="numeric"),
         contains="gseaResItem")

setClass("annoGseaRes", contains="list")
setClass("annoGseaResList", contains="list") #3 a list of annoGseaRes objects

setClass("GeneSetResult",
         representation=list(
             gsCategory="character",
             gsName="character",
             gsEffSize="integer",
             p="numeric",
             fdr="numeric"),
         contains="VIRTUAL")
            
## Fisher's exact test
setClass("FisherResult",
         representation=list(hits="character"),
         contains="GeneSetResult")

setClass("FisherResultList",
         representation=list(
             inputName="character",
             input="character",
             universe="character"),
         contain="list")


## Camera result
setClass("CameraResult",
         representation=list(
             correlation="numeric",
             hits="character",
             score="numeric"),
         contains="GeneSetResult")
setClass("CameraResultList",
         representation=list(inputName="character"),
         contain="list")
            
