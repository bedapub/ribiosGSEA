#' @importFrom grDevices colorRampPalette
#' @importFrom graphics plot
#' @importFrom stats ave cor filter fisher.test kmeans median p.adjust pchisq phyper pnorm pt qt rnorm runif sd t.test update var
#' @importFrom utils read.table
#' @importFrom methods setClass setGeneric setMethod as is new callGeneric
#' @importFrom XML xmlAttrs xmlTreeParse xmlRoot xmlChildren xmlName
#' @importFrom parallel mclapply
#' @importFrom Rcpp evalCpp
#' @importFrom limma zscoreT squeezeVar lmFit contrasts.fit eBayes topTable camera mroast mroast.default romer
#' @importFrom MASS mvrnorm
#' @importFrom lattice xyplot panel.xyplot trellis.par.get
#' @importFrom BioQC wmwTest readGmt GmtList gsName gsDesc gsGenes gsGeneCount gsSize gsNamespace gsNamespace<- appendGmtList
#' @importFrom globaltest gt
#' @importFrom Biobase exprs exprs<-
#' @importFrom ribiosExpression designMatrix contrastMatrix truncateDgeTable limmaTopTable2dgeTable
#' @importFrom ribiosUtils haltifnot uniqueNonNA whoami
#' @importFrom data.table data.table
#' @importFrom dplyr %>% group_by ungroup arrange select
#' @importFrom httr POST GET content http_error
#' @importFrom jsonlite fromJSON
#' @import ribiosIO
#' @import ribiosMath
#' @import ribiosUtils
#' @export geneSetPerm gsGenes<- gsGeneValues<- isGseaCoreEnrich gseaES gseaNES gseaNP gseaFDR gseaFWER gsGeneIndices gseaESprofile gseaCoreEnrichThr gseaCoreEnrichGenes gseaLeadingEdgeGenes annoGseaResItem annoGseaRes gseaScore gseaScores filterBySize gsNames gsDescs GSEA_DATA_DIR GSEA_ANNOTATION_DIR GSEA_GENESET_DIR JAVA_BIN GSEA_JAR DEFAULT_GMT DEFAULT_CHIP gseaData gseaAnno gseaGeneSet dirGseaGeneSet lsGseaGeneSet javaBin gseaJar defaultGmt defaultChip gseaResQvalue gseaResES gseaFingerprint gseaFingerprintMatrix buildGSEAcomm readMPSGmt writeGmt readDefaultGenesets gmthyper gmthyperList hits pValue fdrValue gsSize gsEffSize minPValue minFdrValue sigGeneSet sigGeneSetTable topGeneSetTable topOrSigGeneSetTable gsFisherTestCore gsListFisherTestCore fisherTest fisherTestEdgeResult list2mat ntpTemplates ntpBiTemplates ntp parseDTG parseGSEAdir parseGSEAedb parseGSEAres biosCamera gmtListCamera nGenes nSamples tpGeneSetInd deltaMean tpGeneSetCor bgDgeInd bgDgeLength bgDgePerc bgDgeDeltaMean bgCorLength bgCorInd bgCorPerc bgCorCluster bgCorSigma randomSeed simulators genesets pFunc pValues generateBenchmarkData benchmark ROC AUC ranks avgRank minRank maxRank rankStat TPR FPR nGenes<- nSamples<- tpGeneSetInd<- deltaMean<- tpGeneSetCor<- bgDgeInd<- bgDgeDeltaMean<- bgCorInd<- bgCorCluster<- bgCorSigma<- randomSeed<- pFunc<- as.matrix.TwoGroupExprsSimulator avgCor mutateBgByParams defaultBgDgeDeltaMeanFunc defaultByCorSigmaFunc randomlyMutateBg newTwoGroupExprsSimulator pValues2BenchmarkResult pFuncFisherExact pFuncCamera pFuncCameraRank pFuncBioQCtStat pFuncMroast pFuncRomer pFuncGlobaltest pFuncFisherMethod pFuncTstatZtest pFuncTstatTtest pFuncTstatWMW pFuncTstatChisq pFuncLimmaAggregated fisherMethod newBenchmarker varParPerformance plotROC plotAUC plotRanks tpDiff isGeMSReachable getJsonResponse insertGmtListToGeMS removeFromGeMS getUserSetsFromGeMS fishersMethod getPvalCol getFDRCol kmeansGeneset ronetGeneSetNamespace readRonetGmt prettyRonetGenesetNames
#' @exportClass TwoGroupExprsSimulator BenchmarkDataset Benchmarker BenchmarkResult CameraResult CameraResultList FisherResult gseaResItem annoGseaResItem annoGseaRes annoGseaResList
#' @exportMethod [
NULL
