### R code from vignette source 'PharmacoGx.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: setup
###################################################
options(keep.source=TRUE)


###################################################
### code chunk number 2: install-pkg (eval = FALSE)
###################################################
## source('http://bioconductor.org/biocLite.R')
## biocLite('PharmacoGx')


###################################################
### code chunk number 3: loadlib (eval = FALSE)
###################################################
## library(PharmacoGx)


###################################################
### code chunk number 4: download-example (eval = FALSE)
###################################################
## ## Example
## CGP <- downloadPSet("CGP") 


###################################################
### code chunk number 5: download-sig (eval = FALSE)
###################################################
## ## Example
## CGP.sigs <- downloadSignatures("CGP")


###################################################
### code chunk number 6: load-examples (eval = FALSE)
###################################################
##   ##Using the included example datasets
##   library(PharmacoGx)
##   data("CGPsmall")
##   data("CCLEsmall")
##   CGPsmall <- probeGeneMapping(CGPsmall) 
##   CCLEsmall <- probeGeneMapping(CCLEsmall) 


###################################################
### code chunk number 7: inconsistencies (eval = FALSE)
###################################################
##   library(PharmacoGx)
##   data("CGPsmall")
##   data("CCLEsmall")
##   CGPsmall <- probeGeneMapping(CGPsmall) 
##   CCLEsmall <- probeGeneMapping(CCLEsmall) 
##   common <- intersectPSet(list('CCLE'=CCLEsmall,
##                                'CGP'=CGPsmall),
##                           intersectOn=c("cell.lines", "drugs", 'genes'))
##   
##   CGP.auc <- summarizeSensitivityPhenotype(
##                 common$CGP, 
##                 sensitivity.measure='auc_published', 
##                 summaryStat="median")
##   CCLE.auc <- summarizeSensitivityPhenotype(
##                 common$CCLE, 
##                 sensitivity.measure='auc_published', 
##                 summaryStat="median")
##   
##   CGP.ic50 <- summarizeSensitivityPhenotype(
##                 common$CGP, 
##                 sensitivity.measure='ic50_published', 
##                 summaryStat="median")
##   CCLE.ic50 <- summarizeSensitivityPhenotype(
##                 common$CCLE, 
##                 sensitivity.measure='ic50_published', 
##                 summaryStat="median")
##   
##   common$CGP <- summarizeGeneExpression(common$CGP, 
##                                         cellNames(common$CGP),
##                                         verbose=FALSE)
##   common$CCLE <- summarizeGeneExpression(common$CCLE, 
##                                          cellNames(common$CCLE),
##                                          verbose=FALSE)
##   gg <- geneNames(common[[1]])
##   cc <- cellNames(common[[1]])
##   
##   ge.cor <- sapply(cc, function (x, d1, d2) {
##     return (stats::cor(d1[ , x], d2[ , x], method="spearman",
##                 use="pairwise.complete.obs"))
##   }, d1=rnaData(common$CGP), d2=rnaData(common$CCLE))
##   ic50.cor <- sapply(cc, function (x, d1, d2) {
##     return (stats::cor(d1[, x], d2[ , x], method="spearman",
##                 use="pairwise.complete.obs"))
##   }, d1=t(CGP.ic50), d2=t(CCLE.ic50))
##   auc.cor <- sapply(cc, function (x, d1, d2) {
##     return (stats::cor(d1[ , x], d2[ , x], method="spearman",
##                 use="pairwise.complete.obs"))
##   }, d1=t(CGP.auc), d2=t(CCLE.auc))
##   
##   w1 <- stats::wilcox.test(x=ge.cor, y=auc.cor, conf.int=TRUE, exact=FALSE)
##   w2 <- stats::wilcox.test(x=ge.cor, y=ic50.cor, conf.int=TRUE, exact=FALSE)
##   yylim <- c(-1, 1)
##   ss <- sprintf("GE vs. AUC = %.1E\nGE vs. IC50 = %.1E",
##                 w1$p.value, w2$p.value)
##   boxplot(list("GE"=ge.cor, "AUC"=auc.cor, "IC50"=ic50.cor),
##           main="Concordance between cell lines",
##           ylab=expression(R[s]),
##           sub=ss,
##           ylim=yylim,
##           col="lightgrey",
##           pch=20,
##           border="black")
## 


###################################################
### code chunk number 8: fig2
###################################################

  library(PharmacoGx)
  data("CGPsmall")
  data("CCLEsmall")
  CGPsmall <- probeGeneMapping(CGPsmall) 
  CCLEsmall <- probeGeneMapping(CCLEsmall) 
  common <- intersectPSet(list('CCLE'=CCLEsmall,
                               'CGP'=CGPsmall),
                          intersectOn=c("cell.lines", "drugs", 'genes'))
  
  CGP.auc <- summarizeSensitivityPhenotype(
                common$CGP, 
                sensitivity.measure='auc_published', 
                summaryStat="median")
  CCLE.auc <- summarizeSensitivityPhenotype(
                common$CCLE, 
                sensitivity.measure='auc_published', 
                summaryStat="median")
  
  CGP.ic50 <- summarizeSensitivityPhenotype(
                common$CGP, 
                sensitivity.measure='ic50_published', 
                summaryStat="median")
  CCLE.ic50 <- summarizeSensitivityPhenotype(
                common$CCLE, 
                sensitivity.measure='ic50_published', 
                summaryStat="median")
  
  common$CGP <- summarizeGeneExpression(common$CGP, 
                                        cellNames(common$CGP),
                                        verbose=FALSE)
  common$CCLE <- summarizeGeneExpression(common$CCLE, 
                                         cellNames(common$CCLE),
                                         verbose=FALSE)
  gg <- geneNames(common[[1]])
  cc <- cellNames(common[[1]])
  
  ge.cor <- sapply(cc, function (x, d1, d2) {
    return (stats::cor(d1[ , x], d2[ , x], method="spearman",
                use="pairwise.complete.obs"))
  }, d1=rnaData(common$CGP), d2=rnaData(common$CCLE))
  ic50.cor <- sapply(cc, function (x, d1, d2) {
    return (stats::cor(d1[, x], d2[ , x], method="spearman",
                use="pairwise.complete.obs"))
  }, d1=t(CGP.ic50), d2=t(CCLE.ic50))
  auc.cor <- sapply(cc, function (x, d1, d2) {
    return (stats::cor(d1[ , x], d2[ , x], method="spearman",
                use="pairwise.complete.obs"))
  }, d1=t(CGP.auc), d2=t(CCLE.auc))
  
  w1 <- stats::wilcox.test(x=ge.cor, y=auc.cor, conf.int=TRUE, exact=FALSE)
  w2 <- stats::wilcox.test(x=ge.cor, y=ic50.cor, conf.int=TRUE, exact=FALSE)
  yylim <- c(-1, 1)
  ss <- sprintf("GE vs. AUC = %.1E\nGE vs. IC50 = %.1E",
                w1$p.value, w2$p.value)
  boxplot(list("GE"=ge.cor, "AUC"=auc.cor, "IC50"=ic50.cor),
          main="Concordance between cell lines",
          ylab=expression(R[s]),
          sub=ss,
          ylim=yylim,
          col="lightgrey",
          pch=20,
          border="black")

    



###################################################
### code chunk number 9: load-cmap
###################################################
  library(PharmacoGx)
  require(xtable)
  data(CMAPsmall)
  drug.perturbation <- drugPertubrationSig(CMAPsmall)
  data(HDAC_genes)
  
  res <- apply(drug.perturbation[,,c("tstat", "fdr")], 2, function(x, HDAC){ 
	    return(connectivityScore(x=x, 
	                             y=HDAC[,2,drop=FALSE], 
	                             method="gsea", nperm=100))
	}, HDAC=HDAC_genes)
  
  rownames(res) <- c("Connectivity", "P Value")
  res <- t(res)
  res <- res[order(res[,1], decreasing=T),]
  xtable(res, 
    caption='Connectivity Score results for HDAC inhibitor gene signature.')


###################################################
### code chunk number 10: sessionInfo
###################################################
toLatex(sessionInfo())


