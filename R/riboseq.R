
gseaDataDir <- function() {
  return("/data64/bi/httpd_8080/apps/gsea")
}
gseaGeneSetDir <- function() {
  return(file.path(gseaDataDir(),
                   "genesets"))
}
listGseaGeneSets <- function() {
  dir(gseaGeneSetDir(), pattern="*.gmt")
}
gseaGeneSetFile <- function(file.name) {
  file <- file.path(gseaGeneSetDir(),
                    file.name)
  if(!file.exists(file)) {
    stop(file.name, " not found at ", file, "!")
  }
  return(file)
}
readGseaGeneSet <- function(file.name) {
  read_gmt_list(gseaGeneSetFile(file.name))
}


geneSetMedianTest <- function (index, statistics, alternative = "mixed", type = "auto", 
    ranks.only = TRUE, nsim = 9999) 
{
    alternative <- match.arg(alternative, c("mixed", "either", 
        "down", "up", "less", "greater", "two.sided"))
    if (alternative == "two.sided") 
        alternative <- "either"
    if (alternative == "less") 
        alternative <- "down"
    if (alternative == "greater") 
        alternative <- "up"
    type <- match.arg(tolower(type), c("auto", "t", "f"))
    allsamesign <- all(statistics >= 0) || all(statistics <= 
        0)
    if (type == "auto") {
        if (allsamesign) 
            type <- "f"
        else type <- "t"
    }
    if (type == "f" & alternative != "mixed") 
        stop("Only alternative=\"mixed\" is possible with F-like statistics.")
    if (alternative == "mixed") 
        statistics <- abs(statistics)
    if (alternative == "down") {
        statistics <- -statistics
        alternative <- "up"
    }
    if (ranks.only) {
        pvalues <- rankSumTestWithCorrelation(index = index, 
            statistics = statistics, df = Inf)
        p.value <- switch(alternative, down = pvalues["less"], 
            up = pvalues["greater"], either = 2 * min(pvalues), 
            mixed = pvalues["greater"])
    }
    else {
        ssel <- statistics[index]
        ssel <- ssel[!is.na(ssel)]
        nsel <- length(ssel)
        if (nsel == 0) 
            return(1)
        stat <- statistics[!is.na(statistics)]
        msel <- median(ssel)
        if (alternative == "either") 
            posstat <- abs
        else posstat <- function(x) x
        msel <- posstat(msel)
        ntail <- 1
        for (i in 1:nsim) if (posstat(median(sample(stat, nsel))) >= 
            msel) 
            ntail <- ntail + 1
        p.value <- ntail/(nsim + 1)
    }
    as.vector(p.value)
}



fisher.method <- function(p) {
  Xsq <- -2*sum(log(p))
  p.val <- pchisq(Xsq, df = 2*length(p),lower.tail=FALSE)
  return(c(Xsq = Xsq, p.value = p.val))
}

pair.fisher.method <- function(p1, p2) {
  stopifnot(length(p1)==length(p2))
  res <- data.frame(t(sapply(seq(along=p1), function(i)
                             fisher.method(c(p1[i], p2[i])))))
}

finite.t.test.pvalue <- function(x,y, ...) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if(length(x)<2 || length(y)<2)
    return(NA)
  return(t.test(x, y, ...)$p.value)
}
gage <- function(rnks, ind.list, alternative="two.sided") {
  sapply(ind.list, function(x) {
    finite.t.test.pvalue(rnks[x], rnks[-x], alternative=alternative)
  })
}
## doGeneSetTest using a simple data structure (a 2-column data frame of gene names and statistics)
doGeneSetTest <- function(rnk, anno,
                          gmt, minGene=NULL, maxGene=NULL, nsim=9999, ...) {
  stopifnot(is.data.frame(rnk) & ncol(rnk)>=2 & is.numeric(rnk[,2]))
  if(is.null(minGene)) minGene <- 0
  if(is.null(maxGene)) maxGene <- Inf
  sorted.anno <- matchColumn(rnk[,1L], anno, 0L)
  sorted.gs <- sorted.anno$GeneSymbol
  all.index <- lapply(gmt, function(x) na.omit(match(x$genes,
                                                     sorted.gs)))

  all.size <- sapply(gmt, function(x) length(x$genes)) 
  all.effsize <- sapply(all.index, function(x) length(x))
  hasGenes <- all.effsize >= minGene & all.effsize <= maxGene
  sub.index <- all.index[hasGenes]
  ## gs <- sapply(sub.index, function(x) geneSetTest(x, rnk[,2L], alternative="either"))
  ## gs.sim <- sapply(sub.index, function(x) geneSetMedianTest(x, rnk[,2L], alternative="either", ranks.only=FALSE, nsim=nsim))
  gs.sim <- gage(rnk[,2], ind.list=sub.index, alternative="two.sided")
  gs.wmw <- BioQC::wmwTest(rnk[,2], ind.list=sub.index, valType="p.two.sided")
  gs.wmw.up <- BioQC::wmwTest(rnk[,2], ind.list=sub.index, valType="p.greater")
  gs.wmw.down <- BioQC::wmwTest(rnk[,2], ind.list=sub.index, valType="p.less")
  gs.comb <- pair.fisher.method(gs.sim, gs.wmw)$p.value
  gs.median <- sapply(sub.index, function(x) median(rnk[x,2]))
  gs.direction <- ifelse(gs.wmw.up<=gs.wmw.down, "up", "down")
  ## wmw bootstrapping (almost identical to gs.wmw)
  ##sub.index.len <- sapply(sub.index, length)
  ##sub.index.ulen <- unique(sub.index.len)
  ##Nboot <- 99999
  ##sub.index.boot <- lapply(seq(sub.index.ulen),
  ##                         function(i)
  ##                         lapply(1:Nboot, function(x) {
  ##                           sample(1:nrow(rnk), sub.index.ulen[i], replace=TRUE)
  ##                         }))
  ##sub.index.boot.wmw <- lapply(sub.index.boot, function(x)
  ##                             wmwTest(rnk[,2], x, alternative="two.sided"))
  ##gs.wmw.boot.ind <- match(sub.index.len, sub.index.ulen)
  ##gs.wmw.sim <- sapply(seq(along=sub.index),
  ##                     function(i) mean(gs.wmw[i] >= sub.index.boot.wmw[[gs.wmw.boot.ind[i]]]))
  res.raw <- data.frame(geneset=names(sub.index),
                        median=gs.median,
                        direction=gs.direction,
                        p.sim=gs.sim,
                        p.wmw=gs.wmw,
                        p.comb=gs.comb,
                        row.names=NULL)
  gs.info <- data.frame(geneset=names(all.index),
                        size=all.size,
                        effective.size=all.effsize, row.names=NULL)
  res <- merge(gs.info, res.raw, by="geneset", all.x=TRUE)
  return(res)
}

doGeneSetTests <- function(rnks, anno,
                           gmt, minGene=NULL, maxGene=NULL, nsim=9999, ...) {
  resl <- lapply(rnks, function(rnk) {
    doGeneSetTest(rnk, anno=anno,
                  gmt=gmt, minGene=minGene, maxGene=maxGene, nsim=nsim,...)
  })
  res <- cbind(contrast=rep(names(rnks), sapply(resl, nrow)),
               do.call(rbind, resl))
  rownames(res) <- NULL
  return(res)
}
