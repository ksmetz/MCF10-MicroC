
# INTIALIZE -----------------------
library(GenomicRanges)
library(mariner)
library(matrixStats)
library(rtracklayer)
library(InteractionSet)
library(scales)


# READ IN -----------------------

## Read in loops -------
loops = as_ginteractions(readRDS("./output/loopProcessing/loops.rds"))

## Read in DEGs  ------
genes = GRanges(readRDS("./output/rnaProcessing/genes.rds"))

## Read in gene loops ------
geneLoops = readRDS("./output/loopGeneOverlap/geneLoops.rds")

## Read in enhancers -------
enh = readRDS("./output/H3K27acProcessing/enhancers.rds")

## Read in repressors -------
rep = GRanges(readRDS("./output/H3K27me3Processing/repressors.rds"))


# LOOP EP CLASSIFICATION -----------

## Find expressed genes -----
smallestGroupSize <- 3
cutoff = 100 # cutoff of 10 = 17,185 expressed genes
expressed <- rowSums(as.matrix(mcols(genes)[,21:29]) >= cutoff) >= smallestGroupSize
sum(expressed) # 12,666 out of 60,660

## Promoter anchors ------
ov1 = findOverlaps(query = anchors(loops, "first"),
                   subject = promoters(genes[expressed],
                                       upstream = 10000,
                                       downstream = 5000))
ov2 = findOverlaps(query = anchors(loops, "second"),
                   subject = promoters(genes[expressed],
                                       upstream = 10000,
                                       downstream = 5000))
PPloops = unique(queryHits(ov1)[queryHits(ov1) %in% queryHits(ov2)])

## EP anchors ----
ov3 = findOverlaps(query = anchors(loops, "first"),
                   subject = enh)
ov4 = findOverlaps(query = anchors(loops, "second"),
                   subject = enh)

EPloops = unique(c(queryHits(ov3)[queryHits(ov3) %in% queryHits(ov2)],
                   queryHits(ov4)[queryHits(ov4) %in% queryHits(ov1)]))
EPloops = EPloops[!(EPloops %in% PPloops)]

## Enhancer anchors -----
EEloops = unique(queryHits(ov3)[queryHits(ov3) %in% queryHits(ov4)])
EEloops = EEloops[!(EEloops %in% PPloops)]
EEloops = EEloops[!(EEloops %in% EPloops)]

## Single-feature anchors -------
PXloops = unique(c(queryHits(ov1)[-(queryHits(ov1) %in% queryHits(ov2))],
                   queryHits(ov2)[-(queryHits(ov2) %in% queryHits(ov1))]))
PXloops = PXloops[!(PXloops %in% PPloops)]
PXloops = PXloops[!(PXloops %in% EPloops)]

EXloops = unique(c(queryHits(ov3)[-(queryHits(ov3) %in% queryHits(ov4))],
                   queryHits(ov4)[-(queryHits(ov4) %in% queryHits(ov3))]))
EXloops = EXloops[!(EXloops %in% PPloops)]
EXloops = EXloops[!(EXloops %in% EPloops)]
EXloops = EXloops[!(EXloops %in% EEloops)]
EXloops = EXloops[!(EXloops %in% PXloops)]

## PLOT: EP pie chart -------
dat = rep("none", length(loops))
dat[PPloops] = "PP"
dat[EPloops] = "EP"
dat[EEloops] = "EE"
dat[EXloops] = "E"
dat[PXloops] = "P"

pdf("./output/enhancerPromoterLoops/Fig2A-EPpie.pdf",
    width = 6, height = 6)

pie(table(dat)[c(4, 5, 6, 3, 2, 1)],
    border = NA, 
    init.angle = -90,
    col = c("grey",
            "khaki2", "gold2",
            "orange2",
            "firebrick3", "pink2"))
dev.off()

## PLOT: Length boxplot -------
pdf("./output/enhancerPromoterLoops/Fig2B-EPlengthBoxplot.pdf",
    width = 6, height = 6)

dat = list("none" = loops$span[-c(PPloops, EPloops, EEloops,
                                  EXloops, PXloops)],
           "P" = loops$span[PXloops],
           "PP" = loops$span[PPloops],
           "EP" = loops$span[EPloops],
           "EE" = loops$span[EEloops],
           "E" = loops$span[EXloops])
boxplot(dat, outline = F,
        col = c("grey",
                "khaki2", "gold2",
                "orange2",
                "firebrick3", "pink2"),
        main = "Loop sizes")
t.test(dat$PP, dat$P, alternative = "less")
t.test(dat$EP, dat$P, alternative = "less")
t.test(dat$EP, dat$E, alternative = "less")
t.test(dat$EE, dat$E, alternative = "less")
t.test(dat$P, dat$E, alternative = "less")
t.test(dat$E, dat$none, alternative = "less")
t.test(dat$PP, dat$EP, alternative = "less")
t.test(dat$PP, dat$EE, alternative = "less")
t.test(dat$EP, dat$EE, alternative = "less") ## Not signif

dev.off()


## PLOT: Count boxplot -------
pdf("./output/enhancerPromoterLoops/Fig2C-EPcountBoxplot.pdf",
    width = 6, height = 6)

dat = list("none" = loops$avgCount[-c(PPloops, EPloops, EEloops,
                                      EXloops, PXloops)],
           "P" = loops$avgCount[PXloops],
           "PP" = loops$avgCount[PPloops],
           "EP" = loops$avgCount[EPloops],
           "EE" = loops$avgCount[EEloops],
           "E" = loops$avgCount[EXloops])
boxplot(dat, outline = F,
        col = c("grey",
                "khaki2", "gold2",
                "orange2",
                "firebrick3", "pink2"),
        main = "Loop counts")
t.test(dat$PP, dat$P, alternative = "greater")
t.test(dat$EP, dat$P, alternative = "greater")
t.test(dat$EP, dat$E, alternative = "greater")
t.test(dat$EE, dat$E, alternative = "greater")
t.test(dat$P, dat$E, alternative = "greater") ## Not signif
t.test(dat$E, dat$none, alternative = "greater")
t.test(dat$PP, dat$EP, alternative = "less") ## Swapped
t.test(dat$PP, dat$EE, alternative = "less") ## Swapped
t.test(dat$EP, dat$EE, alternative = "greater") ## Not signif

dev.off()


# DEG EXPLANATIONS ---------------
overlapSet <- function(geneSet, enhSet, pairSet){
  
  # Find a gene set of interest
  fullSet = unique(geneSet$SYMBOL)
  fullSet = fullSet[fullSet != ""]
  
  # Find promoters + enhancer H3K27ac peaks
  ov = findOverlaps(query = enhSet,
                    subject = promoters(geneSet,
                                        upstream = 10000,
                                        downstream = 5000))
  prom = enhSet[queryHits(ov),]
  nonProm = enhSet[-queryHits(ov),]
  
  # Set 1: Genes with differential H3K27ac in promoters
  set1 = unique(geneSet$SYMBOL[subjectHits(ov)])
  
  # Set 2: Genes with differential H3K27ac in distal (looped) enhancers
  ov = findOverlaps(query = geneLoops,
                    subject = nonProm)
  set2 = unique(geneLoops$SYMBOL[queryHits(ov)])
  
  # Find diff loops with enhancers
  set3 = unique(pairSet$SYMBOL)
  
  # Find overlaps
  dat = table(fullSet %in% set1,
              fullSet %in% set2,
              fullSet %in% set3)
  
  return(dat)
}


setStats <- function(geneSet, enhSet, pairSet,
                     geneCol, enhCol, pairCol){
  
  # Find a gene set of interest
  fullSet = unique(geneSet$SYMBOL)
  fullSet = fullSet[fullSet != ""]
  
  # Find promoters + enhancer H3K27ac peaks
  ov = findOverlaps(query = enhSet,
                    subject = promoters(geneSet,
                                        upstream = 10000,
                                        downstream = 5000))
  prom = enhSet[queryHits(ov),]
  nonProm = enhSet[-queryHits(ov),]
  
  # Set 1: Genes with differential H3K27ac in promoters
  # set1 = fullSet[fullSet %in% geneSet$SYMBOL[subjectHits(ov)]]
  dat1 = table(sign(mcols(geneSet)[subjectHits(ov), geneCol]),
               sign(mcols(enhSet)[queryHits(ov), enhCol]))
  
  # Set 2: Genes with differential H3K27ac in distal (looped) enhancers
  ov = findOverlaps(query = geneLoops,
                    subject = nonProm, 
                    maxgap = 10000)
  # set2 = unique(geneLoops[unique(queryHits(ov)),])
  dat2 = table(sign(mcols(geneLoops)[queryHits(ov), geneCol]),
               sign(mcols(nonProm)[subjectHits(ov), enhCol]))
  if(nrow(dat2) == 3){
    dat2 = dat2[-2,]
  }
  
  # Set 3: Find diff loops with enhancers
  set3 = pairSet[pairSet$SYMBOL %in% fullSet]
  dat3 = table(sign(mcols(set3)[, geneCol]),
               sign(mcols(set3)[, pairCol]))
  
  # Return stats
  return(list("prom" = fisher.test(dat1),
              "enh" = fisher.test(dat2),
              "loop" = fisher.test(dat3),
              "pDat" = dat1,
              "eDat" = dat2,
              "lDat" = dat3))
}

## Overlapping DEGs with differential distal enhancers -------
acUp = overlapSet(genes[genes$AC.sig & genes$AC.log2FoldChange > 0,], 
                  enh[which(enh$AC.log2FoldChange > 0),], 
                  geneLoops[which(geneLoops$AClfc > 0 & geneLoops$ACdiff),])
acDn = overlapSet(genes[genes$AC.sig & genes$AC.log2FoldChange < 0,], 
                  enh[which(enh$AC.log2FoldChange < 0),], 
                  geneLoops[which(geneLoops$AClfc < 0 & geneLoops$ACdiff),])

atUp = overlapSet(genes[genes$AT.sig & genes$AT.log2FoldChange > 0,], 
                  enh[which(enh$AT.log2FoldChange > 0),], 
                  geneLoops[which(geneLoops$ATlfc > 0 & geneLoops$ATdiff),])
atDn = overlapSet(genes[genes$AT.sig & genes$AT.log2FoldChange < 0,], 
                  enh[which(enh$AT.log2FoldChange < 0),], 
                  geneLoops[which(geneLoops$ATlfc < 0 & geneLoops$ATdiff),])

tcUp = overlapSet(genes[genes$TC.sig & genes$TC.log2FoldChange > 0,], 
                  enh[which(enh$TC.log2FoldChange > 0),], 
                  geneLoops[which(geneLoops$TClfc > 0 & geneLoops$TCdiff),])
tcDn = overlapSet(genes[genes$TC.sig & genes$TC.log2FoldChange < 0,], 
                  enh[which(enh$TC.log2FoldChange < 0),], 
                  geneLoops[which(geneLoops$TClfc < 0 & geneLoops$TCdiff),])

allUp = acUp + tcUp + atUp
allDn = acDn + tcDn + atDn

## PLOT: OR Barplots -------
acStat = setStats(geneSet = genes[genes$AC.sig,], 
                  enhSet = enh,  # enh[enh$AC.sig,], 
                  pairSet = geneLoops, # geneLoops[geneLoops$ATdiff,],
                  geneCol = "AC.log2FoldChange",
                  enhCol = "AC.log2FoldChange",
                  pairCol = "AClfc")
atStat = setStats(geneSet = genes[genes$AT.sig,], 
                  enhSet = enh, # enh[enh$AT.sig,], 
                  pairSet = geneLoops, #geneLoops[geneLoops$ATdiff,],
                  geneCol = "AT.log2FoldChange",
                  enhCol = "AT.log2FoldChange",
                  pairCol = "ATlfc")
tcStat = setStats(geneSet = genes[genes$TC.sig,], 
                  enhSet = enh, # enh[enh$TC.sig,], 
                  pairSet = geneLoops, # geneLoops[geneLoops$TCdiff,],
                  geneCol = "TC.log2FoldChange",
                  enhCol = "TC.log2FoldChange",
                  pairCol = "TClfc")

pdf("./output/enhancerPromoterLoops/FigS3E-ORbarplots.pdf",
    width = 6, height = 2.5)

par(mfrow=c(1,3))
barplot(t(atStat$pDat), names.arg = c("-", "+"), xlab = "Gene Change", main = "Gene-Promoter", col = c("orange", "grey"), border = NA)
barplot(t(atStat$eDat), names.arg = c("-", "+"), xlab = "Gene Change", main = "Gene-Enhancer", col = c("firebrick", "grey"), border = NA)
barplot(t(atStat$lDat), names.arg = c("-", "+"), xlab = "Gene Change", main = "Gene-Loop", col = c("steelblue", "grey"), border = NA)

par(mfrow=c(1,3))
barplot(t(tcStat$pDat), names.arg = c("-", "+"), xlab = "Gene Change", main = "Gene-Promoter", col = c("orange", "grey"), border = NA)
barplot(t(tcStat$eDat), names.arg = c("-", "+"), xlab = "Gene Change", main = "Gene-Enhancer", col = c("firebrick", "grey"), border = NA)
barplot(t(tcStat$lDat), names.arg = c("-", "+"), xlab = "Gene Change", main = "Gene-Loop", col = c("steelblue", "grey"), border = NA)

par(mfrow=c(1,3))
barplot(t(acStat$pDat), names.arg = c("-", "+"), xlab = "Gene Change", main = "Gene-Promoter", col = c("orange", "grey"), border = NA)
barplot(t(acStat$eDat), names.arg = c("-", "+"), xlab = "Gene Change", main = "Gene-Enhancer", col = c("firebrick", "grey"), border = NA)
barplot(t(acStat$lDat), names.arg = c("-", "+"), xlab = "Gene Change", main = "Gene-Loop", col = c("steelblue", "grey"), border = NA)

dev.off()


## PLOT: Sig OR Barplots -------
acStat = setStats(geneSet = genes[genes$AC.sig,], 
                  enhSet = enh[enh$AC.sig,], 
                  pairSet = geneLoops[geneLoops$ACdiff,],
                  geneCol = "AC.log2FoldChange",
                  enhCol = "AC.log2FoldChange",
                  pairCol = "AClfc")
atStat = setStats(geneSet = genes[genes$AT.sig,], 
                  enhSet = enh[enh$AT.sig,], 
                  pairSet = geneLoops[geneLoops$ATdiff,],
                  geneCol = "AT.log2FoldChange",
                  enhCol = "AT.log2FoldChange",
                  pairCol = "ATlfc")
tcStat = setStats(geneSet = genes[genes$TC.sig,], 
                  enhSet = enh[enh$TC.sig,], 
                  pairSet = geneLoops[geneLoops$TCdiff,],
                  geneCol = "TC.log2FoldChange",
                  enhCol = "TC.log2FoldChange",
                  pairCol = "TClfc")

pdf("./output/enhancerPromoterLoops/FigS3F-sigORbarplots.pdf",
    width = 6, height = 2.5)

par(mfrow=c(1,3))
barplot(t(atStat$pDat), names.arg = c("-", "+"), xlab = "Gene Change", main = "Gene-Promoter", col = c("orange", "grey"), border = NA)
barplot(t(atStat$eDat), names.arg = c("-", "+"), xlab = "Gene Change", main = "Gene-Enhancer", col = c("firebrick", "grey"), border = NA)
barplot(t(atStat$lDat), names.arg = c("-", "+"), xlab = "Gene Change", main = "Gene-Loop", col = c("steelblue", "grey"), border = NA)

par(mfrow=c(1,3))
barplot(t(tcStat$pDat), names.arg = c("-", "+"), xlab = "Gene Change", main = "Gene-Promoter", col = c("orange", "grey"), border = NA)
barplot(t(tcStat$eDat), names.arg = c("-", "+"), xlab = "Gene Change", main = "Gene-Enhancer", col = c("firebrick", "grey"), border = NA)
barplot(t(tcStat$lDat), names.arg = c("-", "+"), xlab = "Gene Change", main = "Gene-Loop", col = c("steelblue", "grey"), border = NA)

par(mfrow=c(1,3))
barplot(t(acStat$pDat), names.arg = c("-", "+"), xlab = "Gene Change", main = "Gene-Promoter", col = c("orange", "grey"), border = NA)
barplot(t(acStat$eDat), names.arg = c("-", "+"), xlab = "Gene Change", main = "Gene-Enhancer", col = c("firebrick", "grey"), border = NA)
barplot(t(acStat$lDat), names.arg = c("-", "+"), xlab = "Gene Change", main = "Gene-Loop", col = c("steelblue", "grey"), border = NA)

dev.off()


## PLOT: DEG explanation pies -------
piePlotter <- function(input, title){
  pie(c(input[c(1, 2, 4, 3)], sum(input[5:8])), 
      border = NA,
      labels = c("None",
                 "Diff prom",
                 "Diff enh + prom",
                 "Diff enh",
                 "Diff loop"), 
      col = c("grey",
              "gold",
              "orange1",
              "firebrick2",
              "steelblue2"),
      main = paste0(title, "\nn=",
                    sum(input)))
}

pdf("./output/enhancerPromoterLoops/Fig2DF-DEGpies.pdf",
    width = 6, height = 6)

par(mfrow=c(1,1))
piePlotter(allUp, "Up-regulated genes")
piePlotter(allDn, "Down-regulated genes")

dev.off()

## PLOT: DEG explanation pies by comparison -------
pdf("./output/enhancerPromoterLoops/FigS3G-DEGpiesByComp.pdf",
    width = 6, height = 6)

piePlotter(atUp, "AvT Up-regulated genes")
piePlotter(atDn, "AvT Down-regulated genes")

piePlotter(tcUp, "TvC Up-regulated genes")
piePlotter(tcDn, "TvC Down-regulated genes")

piePlotter(acUp, "AvC Up-regulated genes")
piePlotter(acDn, "AvC Down-regulated genes")

dev.off()


## Feature boxplots ----------
lfcBoxplot <- function(loopInput, title=""){
  
  par(mfrow=c(1,1))
  ov = findOverlaps(query = loopInput,
                    subject = enh)
  loopEnh = unique(enh[subjectHits(ov),])
  
  ov = findOverlaps(query = loopEnh,
                    subject = promoters(genes[genes$SYMBOL %in% loopInput$SYMBOL],
                                        upstream = 10000,
                                        downstream = 5000))
  enhK27 = unique(loopEnh[-(queryHits(ov)),])
  proK27 = unique(loopEnh[queryHits(ov),])
  
  ov = findOverlaps(query = loopInput,
                    subject = rep)
  k27me3 = unique(rep[subjectHits(ov),])
  
  dat = list("repr" = k27me3$AC.log2FoldChange,
             "enh" = enhK27$AC.log2FoldChange,
             "pro" = proK27$AC.log2FoldChange,
             "gene" = loopInput$AC.log2FoldChange,
             "loop" = loopInput$AClfc)
  dat1 = list("repr" = k27me3$AC.log2FoldChange[!k27me3$AC.sig],
              "enh" = enhK27$AC.log2FoldChange[!enhK27$AC.sig],
              "pro" = proK27$AC.log2FoldChange[!proK27$AC.sig],
              "gene" = loopInput$AC.log2FoldChange[!loopInput$AC.sig],
              "loop" = loopInput$AClfc[!loopInput$ACdiff])
  dat2 = list("repr" = k27me3$AC.log2FoldChange[k27me3$AC.sig],
              "enh" = enhK27$AC.log2FoldChange[enhK27$AC.sig],
              "pro" = proK27$AC.log2FoldChange[proK27$AC.sig],
              "gene" = loopInput$AC.log2FoldChange[loopInput$AC.sig],
              "loop" = loopInput$AClfc[loopInput$ACdiff])
  boxplot(dat,
          col = alpha("grey", 0.5),
          outline = F,
          las = 1,
          main = title)
  stripchart(dat1,
             col = alpha("grey", 0.5),
             vertical = T, 
             method = "jitter",
             jitter = 0.25,
             pch = 19,
             add = T)
  stripchart(dat2,
             col = alpha(c("grey30", "firebrick", "orange", "gold", "steelblue"),
                         0.5),
             vertical = T, 
             method = "jitter",
             jitter = 0.25,
             pch = 19,
             add = T)
  
  boxplot(dat,
          col = NA,
          outline = F,
          las = 1,
          add = T)
  abline(h=0)
  
  return(list("repr" = t.test(dat$repr),
              "reprSig" = t.test(dat2$repr),
              "enh"  = t.test(dat$enh),
              "enhSig"  = t.test(dat2$enh),
              "pro"  = t.test(dat$pro),
              "proSig"  = t.test(dat2$pro),
              "gene" = t.test(dat$gene),
              "geneSig" = t.test(dat2$gene),
              "loop" = t.test(dat$loop),
              "loopSig" = t.test(dat2$loop)))
  
}



## PLOT: DEG feature boxplots ---------
pdf("./output/enhancerPromoterLoops/Fig2EG-DEGfeatureBoxplots.pdf",
    width = 8, height = 6)

upStatic = geneLoops[geneLoops$AC.sig & 
                       geneLoops$AC.log2FoldChange > 0,]

upregStats = lfcBoxplot(loopInput = upStatic, title = "Upregulated genes")


downStatic = geneLoops[geneLoops$AC.sig &
                         geneLoops$AC.log2FoldChange < 0,]

dnregStats = lfcBoxplot(loopInput = downStatic, title = "Downregulated genes")

dev.off()



# OUTPUT ---------------
EPstatus = rep("none", length(loops))
EPstatus[PPloops] = "PP"
EPstatus[EPloops] = "EP"
EPstatus[EEloops] = "EE"
EPstatus[EXloops] = "E"
EPstatus[PXloops] = "P"

saveRDS(EPstatus, file = "./output/enhancerPromoterLoops/loop-EPstatus.rds")
