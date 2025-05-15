
# INTIALIZE -----------------------
library(GenomicRanges)
library(scales)
library(GenomicFeatures)
library(gprofiler2)
library(org.Hs.eg.db)
library(colorspace)
library(mariner)


# READ IN -----------------------

## Micro-C loops -------
loops = as_ginteractions(readRDS("./output/loopProcessing/loops.rds"))

## Loops with differential genes -------
DEGloops = readRDS("./output/loopGeneOverlap/DEGloops.rds")
DEGdiffloops = readRDS("./output/loopGeneOverlap/DEGdiffLoops.rds")

## Genes -------
genes = GRanges(readRDS("./output/rnaProcessing/genes.rds"))

## Enhancers (H3K27ac) -------
enh = GRanges(readRDS(file = "./output/H3K27acProcessing/enhancers.rds"))

## Repressors (H3K27me3) -------
repr = GRanges(readRDS("./output/H3K27me3Processing/repressors.rds"))


# ID DEG-DIFF LOOP PAIRS -----------------------
ATloopDEGs = DEGdiffloops[DEGdiffloops$ATdiff == T & DEGdiffloops$AT.sig == T,]
TCloopDEGs = DEGdiffloops[DEGdiffloops$TCdiff == T & DEGdiffloops$TC.sig == T,]
ACloopDEGs = DEGdiffloops[DEGdiffloops$ACdiff == T & DEGdiffloops$AC.sig == T,]

gainedLoopGeneLFC = c(ATloopDEGs$AT.log2FoldChange[ATloopDEGs$ATlfc > 0], 
                      TCloopDEGs$TC.log2FoldChange[TCloopDEGs$TClfc > 0],
                      ACloopDEGs$AC.log2FoldChange[ACloopDEGs$AClfc > 0])
lostLoopGeneLFC = c(ATloopDEGs$AT.log2FoldChange[ATloopDEGs$ATlfc < 0], 
                    TCloopDEGs$TC.log2FoldChange[TCloopDEGs$TClfc < 0],
                    ACloopDEGs$AC.log2FoldChange[ACloopDEGs$AClfc < 0])

tmpSt = DEGloops[DEGloops$cluster == "static",]
staticLoopGeneLFC = c(tmpSt$AC.log2FoldChange[tmpSt$AC.sig], 
                      tmpSt$AT.log2FoldChange[tmpSt$AT.sig], 
                      tmpSt$TC.log2FoldChange[tmpSt$TC.sig])
set.seed(476)
staticLoopGeneLFC = staticLoopGeneLFC[sample(1:length(staticLoopGeneLFC), 
                                             size = max(length(gainedLoopGeneLFC),
                                                        length(lostLoopGeneLFC)))]

## PLOT: Looped-DEG LFC boxplots -------
pdf("./output/loopDEGanalysis/Fig4A-loopedDEGboxplots.pdf",
    width = 6, height = 6)
boxplot(list("gained" = gainedLoopGeneLFC,
             "lost" = lostLoopGeneLFC,
             "static" = staticLoopGeneLFC),
        col = c("#a7d7eb", "#a8edd0", "lightgrey"),
        ylim = c(-10, 10),
        outline = F,
        main = "DEG LFC at Loop Anchors")
stripchart(list("gained" = gainedLoopGeneLFC,
                "lost" = lostLoopGeneLFC,
                "static" = staticLoopGeneLFC),
           col = alpha(c("#327894", "#2bbd80", "darkgrey"),
                       0.25),
           vertical = T, 
           method = "jitter",
           jitter = 0.25,
           pch = 19,
           add = T)
boxplot(list("gained" = gainedLoopGeneLFC,
             "lost" = lostLoopGeneLFC,
             "static" = staticLoopGeneLFC),
        col = NA,
        outline = F, 
        add = T)
abline(h=0)

t.test(gainedLoopGeneLFC)
t.test(lostLoopGeneLFC)
t.test(staticLoopGeneLFC)

dev.off()

## PLOT: Looped-DEG direction OR barplots -------
pdf("./output/loopDEGanalysis/Fig4B-loopedDEGdirectionOR.pdf",
    width = 6, height = 6)

dat1 = table(gene = sign(ATloopDEGs$AT.log2FoldChange), loop = sign(ATloopDEGs$ATlfc))
dat2 = table(gene = sign(TCloopDEGs$TC.log2FoldChange), loop = sign(TCloopDEGs$TClfc))
dat3 = table(gene = sign(ACloopDEGs$AC.log2FoldChange), loop = sign(ACloopDEGs$AClfc))
fisher.test(dat1)
fisher.test(dat2)
fisher.test(dat3)
fisher.test(dat1 + dat2 + dat3)

all = dat1+dat2+dat3
rownames(all) = c("+", "-")
colnames(all) = c("+", "-")
barplot(all[c(2,1),], xlab = "Loop Change", 
        border = NA, main = "Differential genes at differential loops",
        col = c("#E6E6E6", "#327894"))

dev.off()


## PLOT: Looped-DEG direction enrichment -------
proms = promoters(genes, upstream = 2000, downstream = 500)
smallestGroupSize <- 3
cutoff = 10
expressed <- rowSums(as.matrix(mcols(genes)[,21:29]) >= cutoff) >= smallestGroupSize

gap = 10000

pdf("./output/loopDEGanalysis/FigS5BC-loopedDEGdirEnrich.pdf",
    width = 4, height = 4)

### Upregulated genes in gained loops -----
focus = genes[genes$AC.log2FoldChange > 0 & genes$AC.sig,] ## Upregulated genes
obs = ACloopDEGs[ACloopDEGs$AC.log2FoldChange > 0 & ACloopDEGs$AClfc > 0,] ## Upreg genes in gained loops
pool = proms[expressed,] ## Expressed genes
pool = pool[!pool$AC.sig] ## Expressed genes not differential at A v C

exp = lapply(1:1000, function(n){
  ## Random sample of genes (equal in number to # of upregulated genes)
  rand = pool[sample(1:length(pool), size = length(focus)),]
  
  ## How many genes overlap with gained loops
  ov = findOverlaps(query = loops,
                    subject = rand, 
                    maxgap = gap)
  looped = loops[queryHits(ov),]
  return(sum(looped$ACdiff & looped$AClfc > 0))
})

plot(density(unlist(exp)), xlim=c(0, 150), 
     main = "Upregulated genes in gained loops")
abline(v = length(obs), col = "red")

length(focus)
median(unlist(exp), na.rm=T)
length(obs)

### Downregulated genes in gained loops -----
focus = genes[genes$AC.log2FoldChange < 0 & genes$AC.sig,]
obs = ACloopDEGs[ACloopDEGs$AC.log2FoldChange < 0 & ACloopDEGs$AClfc < 0,]

exp = lapply(1:1000, function(n){
  rand = pool[sample(1:length(pool), size = length(focus)),]
  
  ov = findOverlaps(query = loops,
                    subject = rand, 
                    maxgap = gap)
  looped = loops[queryHits(ov),]
  return(sum(looped$ACdiff & looped$AClfc < 0))
})

plot(density(unlist(exp)), xlim=c(0, 150), 
     main = "Downregulated genes in lost loops")
abline(v = length(obs), col = "red")

length(focus)
median(unlist(exp), na.rm=T)
length(obs)

dev.off()


## Total unique genes -------
loopedUp = c(ACloopDEGs[ACloopDEGs$AC.log2FoldChange > 0 & ACloopDEGs$AClfc > 0,],
             ATloopDEGs[ATloopDEGs$AT.log2FoldChange > 0 & ATloopDEGs$ATlfc > 0,],
             TCloopDEGs[TCloopDEGs$TC.log2FoldChange > 0 & TCloopDEGs$TClfc > 0,])
length(unique(loopedUp$GENEID))

loopedDn = c(ACloopDEGs[ACloopDEGs$AC.log2FoldChange < 0 & ACloopDEGs$AClfc < 0,],
             ATloopDEGs[ATloopDEGs$AT.log2FoldChange < 0 & ATloopDEGs$ATlfc < 0,],
             TCloopDEGs[TCloopDEGs$TC.log2FoldChange < 0 & TCloopDEGs$TClfc < 0,])
length(unique(loopedDn$GENEID))


# GO ENRICHMENT ----------------

## Identify gene sets --------
allGenes = genes$GENEID
ATgenes = unique(ATloopDEGs$SYMBOL)
TCgenes = unique(TCloopDEGs$SYMBOL)
ACgenes = unique(ACloopDEGs$SYMBOL)

## Run GO enrichment --------
goAT = gost(ATgenes, organism="hsapiens", ordered_query=F, significant=T, 
            user_threshold=0.05, correction_method="fdr", sources=c("GO", "KEGG","REAC", "WP"),
            custom_bg=allGenes) # running GO
goTC = gost(TCgenes, organism="hsapiens", ordered_query=F, significant=T, 
            user_threshold=0.05, correction_method="fdr", sources=c("GO", "KEGG","REAC", "WP"),
            custom_bg=allGenes) # running GO
goAC = gost(ACgenes, organism="hsapiens", ordered_query=F, significant=T, 
            user_threshold=0.05, correction_method="fdr", sources=c("GO", "KEGG","REAC", "WP"),
            custom_bg=allGenes) # running GO

## Pull top GO hits --------
godat = function(goList, source, n = 30){ # plotting the result
  goList = lapply(goList, function(goterm){
    goterm = goterm$result
    goterm = goterm[goterm$term_size<800, ]
    goterm$geneRatio = goterm$intersection_size/goterm$term_size
    # goterm = goterm[!(goterm$term_id %in% goterm$parents),]
    goterm = goterm[order(goterm$p_value),]
    goterm = goterm[goterm$source %in% source,]
    return(goterm)
  })
  
  topTerms = lapply(goList, function(goterm){
    goterm$term_name[1:n]
  })
  topTerms = unique(unlist(topTerms))
  
  topIDs = lapply(goList, function(goterm){
    goterm$term_id[1:n]
  })
  topIDs = unique(unlist(topIDs))
  
  goterm2plot = lapply(goList, function(goterm){
    goterm2plot = goterm[,c("term_name","p_value","intersection_size","term_id", "geneRatio")]
    goterm2plot = goterm2plot[match(topTerms, goterm2plot$term_name),]
    goterm2plot$log10P = -log10(goterm2plot$p_value)
    goterm2plot$term_name = topTerms
    goterm2plot$term_id = topIDs
    goterm2plot$term_name = factor(goterm2plot$term_name, levels=rev(goterm2plot$term_name))
    goterm2plot[is.na(goterm2plot)] = 0
    return(goterm2plot)
  })
  
  names(goterm2plot) = names(goList)
  return(goterm2plot)
}

datAll = godat(list("AvT" = goAT, 
                    "TvC" = goTC, 
                    "AvC" = goAC), source = "GO:BP", n = 10)

## Manually classify GO terms ------------
sets = data.frame("ID" = datAll$AvT$term_id,
                  "term" = datAll$AvT$term_name,
                  "cat" = "other")
sets$cat[sets$ID %in% paste0("GO:", 
                             c("0007156",
                               "0098742",
                               "0016339",
                               "0034329",
                               "0030155",
                               "0034330",
                               "0022407",
                               "0007268",
                               "0098916"))] = "adhesion"
sets$cat[sets$ID %in% paste0("GO:", 
                             c("0070848",
                               "0071363",
                               "0090287",
                               "0008285",
                               "0071560",
                               "0071559",
                               "0006044"))] = "proliferation"
sets$cat[sets$ID %in% paste0("GO:", 
                             c("0007416",
                               "0048729",
                               "0050808",
                               "0045596",
                               "0001501",
                               "0002009",
                               "0048598",
                               "0035107",
                               "0035108",
                               "0007507",
                               "0032989",
                               "0030855",
                               "0043588",
                               "0008544",
                               "0048589",
                               "0001568",
                               "0060173",
                               "0048736"))] = "development"
cols = c("#1b9e77", "#7570b3", "#d95f02")
sets$col = cols[as.factor(sets$cat)]

## PLOT: Looped-DEG GO enrichment -------
goplot = function(dat, title, rFactor = 50){
  plot(x=0, y=0, type='n',
       yaxt = 'n', xaxt = 'n',
       ylab = NA, xlab = NA,
       bty = 'n',
       xlim = c(-16, length(dat)),
       ylim = c(1, nrow(dat[[1]])),
       xpd = T, main = title)
  for(i in 1:length(dat)){
    radius = sqrt(dat[[i]]$log10P/pi)
    symbols(x = rep(i, nrow(dat[[i]])),
            y = nrow(dat[[i]]):1,
            circles=radius/6,
            inches = F,
            bg = "black",
            fg = NA, add=T)
  }
  axis(side = 1, lwd = 0,
       at = 1:length(dat),
       labels = names(dat))
  text(x = rep(0.5, times = nrow(dat[[1]])),
       y = nrow(dat[[1]]):1,
       labels = dat[[1]]$term_name, adj = c(1, 0.5), 
       col = sets$col,
       xpd = T)
}


pdf("./output/loopDEGanalysis/Fig4C-loopedDEGgoEnrich.pdf",
    width = 8, height = 8)
par(mfrow=c(1,1))
par(mar=c(3,3,3,3))
goplot(datAll, title = "GO term enrichment")
dev.off()


# VISUALIZE DEG-DIFF LOOP PAIRS ----------------
lollipopPlot = function(input, title, loopFCcol, geneFCcol, GOdat,
                        col1 = "red", col2 = "blue"){
  gLFC = mcols(input)[,geneFCcol]
  gLFCsumm = split(gLFC, input$SYMBOL)
  gLFCsumm = lapply(gLFCsumm, mean)
  
  lLFC = mcols(input)[,loopFCcol]
  lLFCsumm = split(lLFC, input$SYMBOL)
  
  o = order(unlist(gLFCsumm), decreasing = T)
  
  datRange = round(range(c(gLFCsumm, lLFC), na.rm=T), digits = 2)
  plot(x = 1:length(gLFCsumm),
       y = gLFCsumm[o],
       type = 'n',
       xaxt = 'n', xlab = '', ylab = "LFC",
       las = 2,
       ylim = datRange,
       main = title)
  
  abline(h=0, col="grey")
  
  lapply(1:length(gLFCsumm), function(n){
    segments(x0 = n,
             x1 = n,
             y0 = gLFCsumm[o][[n]],
             y1 = lLFCsumm[o][[n]],
             col = "grey")
  })
  
  lapply(1:length(lLFCsumm), function(n){
    dat = lLFCsumm[o][[n]]
    points(x = rep(n, times = length(dat)),
           y = dat,
           pch = 19,
           col = alpha("grey", 0.5))
  })
  lapply(1:length(lLFCsumm), function(n){
    ldat = lLFCsumm[o][[n]]
    gdat = gLFCsumm[o][[n]]
    if(all(!is.na(ldat))){
      points(x = rep(n, times=sum(sign(ldat) == sign(gdat))),
             y = ldat[which(sign(ldat) == sign(gdat))],
             col = alpha("black", 0.5), pch = 19)
    }
  })
  
  points(x = 1:length(gLFCsumm),
         y = gLFCsumm[o],
         pch = 19, 
         col = c(lighten(col2, 0.5, space = "HLS"), 
                 lighten(col1, 0.5, space = "HLS"))[as.factor(unlist(lapply(gLFCsumm[o], sign)))])
  points(x = which(unlist(Map(function(G, L){
    any(sign(unlist(G)) == 1 & sign(unlist(L)) == 1)
  }, G = gLFCsumm[o], L = lLFCsumm[o]))),
  y = gLFCsumm[o][which(unlist(Map(function(G, L){
    any(sign(unlist(G)) == 1 & sign(unlist(L)) == 1)
  }, G = gLFCsumm[o], L = lLFCsumm[o])))],
  col = col1, pch = 19)
  points(x = which(unlist(Map(function(G, L){
    any(sign(unlist(G)) == -1 & sign(unlist(L)) == -1)
  }, G = gLFCsumm[o], L = lLFCsumm[o]))),
  y = gLFCsumm[o][which(unlist(Map(function(G, L){
    any(sign(unlist(G)) == -1 & sign(unlist(L)) == -1)
  }, G = gLFCsumm[o], L = lLFCsumm[o])))],
  col = col2, pch = 19)
  
  GOlist = select(org.Hs.eg.db,
                  keys = GOdat$term_id, 
                  keytype = "GOALL",
                  columns = "SYMBOL")
  highlighted = names(gLFCsumm)[o] %in% GOlist$SYMBOL
  text(x = which(!highlighted),
       y = rep(datRange[1]-diff(datRange)*.05, 
               times = sum(!highlighted)),
       labels = names(gLFCsumm)[o][!highlighted],
       srt = 90,
       adj = c(1, 0.5),
       col = "grey",
       cex = 0.5,
       xpd = T)
  
  GOlist = split(GOlist, GOlist$GO)
  lapply(GOlist, function(g){
    terms = names(gLFCsumm)[o] %in% g$SYMBOL
    text(x = which(terms),
         y = rep(datRange[1]-diff(datRange)*.05, 
                 times = length(terms)),
         labels = names(gLFCsumm)[o][terms],
         srt = 90,
         adj = c(1, 0.5),
         col = sets$col[match(unique(g$GOALL), sets$ID)],
         cex = 0.75,
         xpd = T)
  })
  
}


## PLOT: Looped-DEG pairwise -------
pdf("./output/loopDEGanalysis/Fig5A-C_lollipopPlots.pdf",
    width = 12, height = 12)

par(mfrow=c(3,1))
par(mar=c(4,4,4,4))

lollipopPlot(input = ATloopDEGs, 
             GOdat = datAll$AvT,
             title = "DEGs in loop anchors (A v T)",
             geneFCcol = "AT.log2FoldChange",
             loopFCcol = "ATlfc",
             col1 = "orange2",
             col2 = "steelblue3")

lollipopPlot(input = TCloopDEGs, 
             GOdat = datAll$TvC,
             title = "DEGs in loop anchors (T v C)",
             geneFCcol = "TC.log2FoldChange",
             loopFCcol = "TClfc",
             col1 = "firebrick3",
             col2 = "orange2")

lollipopPlot(input = ACloopDEGs, 
             GOdat = datAll$AvC,
             title = "DEGs in loop anchors (A v C)",
             geneFCcol = "AC.log2FoldChange",
             loopFCcol = "AClfc",
             col1 = "firebrick3",
             col2 = "steelblue3")

dev.off()


## Number of corresponding changes per comparison -------
lollipopStats = function(input, loopFCcol, geneFCcol){
  gLFC = mcols(input)[,geneFCcol]
  gLFCsumm = split(gLFC, input$SYMBOL)
  gLFCsumm = lapply(gLFCsumm, mean)
  
  lLFC = mcols(input)[,loopFCcol]
  lLFCsumm = split(lLFC, input$SYMBOL)
  
  o = order(unlist(gLFCsumm), decreasing = T)
  
  tmp = (unlist(Map(function(G, L){
    any(sign(unlist(G)) == sign(unlist(L)))
  }, G = gLFCsumm[o], L = lLFCsumm[o])))
  
  return(tmp)
  
}

## Upregulated in A vs T/C
AvT = lollipopStats(input = ATloopDEGs[ATloopDEGs$AT.log2FoldChange < 0,], 
                    loopFCcol = "ATlfc",geneFCcol = "AT.log2FoldChange")
sum(AvT)/length(AvT)
AvC = lollipopStats(input = ACloopDEGs[ACloopDEGs$AC.log2FoldChange < 0,], 
                    loopFCcol = "AClfc",geneFCcol = "AC.log2FoldChange")
sum(AvC)/length(AvC)

## Upregulated in T vs A/C
TvA = lollipopStats(input = ATloopDEGs[ATloopDEGs$AT.log2FoldChange > 0,], 
                    loopFCcol = "ATlfc",geneFCcol = "AT.log2FoldChange")
sum(TvA)/length(TvA)
TvC = lollipopStats(input = TCloopDEGs[TCloopDEGs$TC.log2FoldChange < 0,], 
                    loopFCcol = "TClfc",geneFCcol = "TC.log2FoldChange")
sum(TvC)/length(TvC)

## Upregulated in C vs A/T
CvA = lollipopStats(input = ACloopDEGs[ACloopDEGs$AC.log2FoldChange > 0,], 
                    loopFCcol = "AClfc",geneFCcol = "AC.log2FoldChange")
sum(CvA)/length(CvA)
CvT = lollipopStats(input = TCloopDEGs[TCloopDEGs$TC.log2FoldChange > 0,], 
                    loopFCcol = "TClfc",geneFCcol = "TC.log2FoldChange")
sum(CvT)/length(CvT)


## Same vs Opp LFC boxplots ----------
lfcBoxplot <- function(input1, input2, title=""){
  par(mfrow=c(1,1))
  ov = findOverlaps(query = input1,
                    subject = enh)
  loopEnh = unique(enh[subjectHits(ov),])
  
  ov = findOverlaps(query = loopEnh,
                    subject = promoters(genes[genes$SYMBOL %in% input1$SYMBOL],
                                        upstream = 10000,
                                        downstream = 5000))
  enhK27 = unique(loopEnh[-(queryHits(ov)),])
  proK27 = unique(loopEnh[queryHits(ov),])
  
  ov = findOverlaps(query = input1,
                    subject = repr)
  k27me3 = unique(repr[subjectHits(ov),])
  
  dat1 = list("repr" = k27me3$AC.log2FoldChange,
              "enh" = enhK27$AC.log2FoldChange,
              "pro" = proK27$AC.log2FoldChange,
              "gene" = input1$AC.log2FoldChange,
              "loop" = input1$AClfc)
  
  
  ov = findOverlaps(query = input2,
                    subject = enh)
  loopEnh = unique(enh[subjectHits(ov),])
  
  ov = findOverlaps(query = loopEnh,
                    subject = promoters(genes[genes$SYMBOL %in% input2$SYMBOL],
                                        upstream = 10000,
                                        downstream = 5000))
  enhK27 = unique(loopEnh[-(queryHits(ov)),])
  proK27 = unique(loopEnh[queryHits(ov),])
  
  ov = findOverlaps(query = input2,
                    subject = repr)
  k27me3 = unique(repr[subjectHits(ov),])
  
  dat2 = list("repr" = k27me3$AC.log2FoldChange,
              "enh" = enhK27$AC.log2FoldChange,
              "pro" = proK27$AC.log2FoldChange,
              "gene" = input2$AC.log2FoldChange,
              "loop" = input2$AClfc)
  
  dat = list(dat1$repr, dat2$repr,
             dat1$enh, dat2$enh,
             dat1$pro, dat2$pro,
             dat1$gene, dat2$gene,
             dat1$loop, dat2$loop)
  
  boxplot(dat,
          col = alpha(rep(c("grey", "firebrick", "orange", "gold", "steelblue"), each = 2),
                      rep(c(0.75, 0.25), times = 5)),
          outline = F,
          # ylim = c(-1, 10),
          main = title)
  stripchart(dat,
             col = alpha(rep(c("grey", "firebrick", "orange", "gold", "steelblue"), each = 2),
                         rep(c(0.5, 0.25), times = 5)),
             vertical = T, 
             method = "jitter",
             jitter = 0.25,
             pch = 19,
             add = T)
  
  boxplot(dat,
          col = NA,
          outline = F,
          add = T)
  abline(h=0)
  
  return(list("repr" = t.test(dat1$repr, dat2$repr),
              "enh"  = t.test(dat1$enh, dat2$enh),
              "pro"  = t.test(dat1$pro, dat2$pro),
              "gene" = t.test(dat1$gene, dat2$gene),
              "loop" = t.test(dat1$loop, dat2$loop)))
}

### PLOT: Looped-DEG feature boxplots -------
pdf("./output/loopDEGanalysis/Fig5DE-loopedDEGfeatureBoxplots.pdf",
    width = 12, height = 6)

bothGained = geneLoops[geneLoops$AC.sig & geneLoops$ACdiff &
                         geneLoops$AC.log2FoldChange > 0 &
                         geneLoops$AClfc > 0,]
oppGained = geneLoops[geneLoops$AC.sig & geneLoops$ACdiff &
                        geneLoops$AC.log2FoldChange > 0 &
                        geneLoops$AClfc < 0,]
lfcBoxplot(input1 = bothGained, input2 = oppGained)

bothLost = geneLoops[geneLoops$AC.sig & geneLoops$ACdiff &
                       geneLoops$AC.log2FoldChange < 0 &
                       geneLoops$AClfc < 0,]
oppLost = geneLoops[geneLoops$AC.sig & geneLoops$ACdiff &
                      geneLoops$AC.log2FoldChange < 0 &
                      geneLoops$AClfc > 0,]
lfcBoxplot(input1 = bothLost, input2 = oppLost)

dev.off()

## Other combinations 
lfcBoxplotAT <- function(input1, input2, title=""){
  par(mfrow=c(1,1))
  ov = findOverlaps(query = input1,
                    subject = enh)
  loopEnh = unique(enh[subjectHits(ov),])
  
  ov = findOverlaps(query = loopEnh,
                    subject = promoters(genes[genes$SYMBOL %in% input1$SYMBOL],
                                        upstream = 10000,
                                        downstream = 5000))
  enhK27 = unique(loopEnh[-(queryHits(ov)),])
  proK27 = unique(loopEnh[queryHits(ov),])
  
  ov = findOverlaps(query = input1,
                    subject = repr)
  k27me3 = unique(repr[subjectHits(ov),])
  
  dat1 = list("repr" = k27me3$AT.log2FoldChange,
              "enh" = enhK27$AT.log2FoldChange,
              "pro" = proK27$AT.log2FoldChange,
              "gene" = input1$AT.log2FoldChange,
              "loop" = input1$ATlfc)
  
  
  ov = findOverlaps(query = input2,
                    subject = enh)
  loopEnh = unique(enh[subjectHits(ov),])
  
  ov = findOverlaps(query = loopEnh,
                    subject = promoters(genes[genes$SYMBOL %in% input2$SYMBOL],
                                        upstream = 10000,
                                        downstream = 5000))
  enhK27 = unique(loopEnh[-(queryHits(ov)),])
  proK27 = unique(loopEnh[queryHits(ov),])
  
  ov = findOverlaps(query = input2,
                    subject = repr)
  k27me3 = unique(repr[subjectHits(ov),])
  
  dat2 = list("repr" = k27me3$AT.log2FoldChange,
              "enh" = enhK27$AT.log2FoldChange,
              "pro" = proK27$AT.log2FoldChange,
              "gene" = input2$AT.log2FoldChange,
              "loop" = input2$ATlfc)
  
  dat = list(dat1$repr, dat2$repr,
             dat1$enh, dat2$enh,
             dat1$pro, dat2$pro,
             dat1$gene, dat2$gene,
             dat1$loop, dat2$loop)
  
  boxplot(dat,
          col = alpha(rep(c("grey", "firebrick", "orange", "gold", "steelblue"), each = 2),
                      rep(c(0.75, 0.25), times = 5)),
          outline = F,
          # ylim = c(-1, 10),
          main = title)
  stripchart(dat,
             col = alpha(rep(c("grey", "firebrick", "orange", "gold", "steelblue"), each = 2),
                         rep(c(0.5, 0.25), times = 5)),
             vertical = T, 
             method = "jitter",
             jitter = 0.25,
             pch = 19,
             add = T)
  
  boxplot(dat,
          col = NA,
          outline = F,
          add = T)
  abline(h=0)
  
  return(list("repr" = t.test(dat1$repr, dat2$repr),
              "enh"  = t.test(dat1$enh, dat2$enh),
              "pro"  = t.test(dat1$pro, dat2$pro),
              "gene" = t.test(dat1$gene, dat2$gene),
              "loop" = t.test(dat1$loop, dat2$loop)))
}

### PLOT: Looped-DEG feature boxplots (AvT) -------
pdf("./output/loopDEGanalysis/FigS5D-loopedDEGfeatureBoxplotsAvT.pdf",
    width = 12, height = 6)

bothGained = geneLoops[geneLoops$AT.sig & geneLoops$ATdiff &
                         geneLoops$AT.log2FoldChange > 0 &
                         geneLoops$ATlfc > 0,]
oppGained = geneLoops[geneLoops$AT.sig & geneLoops$ATdiff &
                        geneLoops$AT.log2FoldChange > 0 &
                        geneLoops$ATlfc < 0,]
lfcBoxplotAT(input1 = bothGained, input2 = oppGained)

bothLost = geneLoops[geneLoops$AT.sig & geneLoops$ATdiff &
                       geneLoops$AT.log2FoldChange < 0 &
                       geneLoops$ATlfc < 0,]
oppLost = geneLoops[geneLoops$AT.sig & geneLoops$ATdiff &
                      geneLoops$AT.log2FoldChange < 0 &
                      geneLoops$ATlfc > 0,]
lfcBoxplotAT(input1 = bothLost, input2 = oppLost)

dev.off()


lfcBoxplotTC <- function(input1, input2, title=""){
  par(mfrow=c(1,1))
  ov = findOverlaps(query = input1,
                    subject = enh)
  loopEnh = unique(enh[subjectHits(ov),])
  
  ov = findOverlaps(query = loopEnh,
                    subject = promoters(genes[genes$SYMBOL %in% input1$SYMBOL],
                                        upstream = 10000,
                                        downstream = 5000))
  enhK27 = unique(loopEnh[-(queryHits(ov)),])
  proK27 = unique(loopEnh[queryHits(ov),])
  
  ov = findOverlaps(query = input1,
                    subject = repr)
  k27me3 = unique(repr[subjectHits(ov),])
  
  dat1 = list("repr" = k27me3$TC.log2FoldChange,
              "enh" = enhK27$TC.log2FoldChange,
              "pro" = proK27$TC.log2FoldChange,
              "gene" = input1$TC.log2FoldChange,
              "loop" = input1$TClfc)
  
  
  ov = findOverlaps(query = input2,
                    subject = enh)
  loopEnh = unique(enh[subjectHits(ov),])
  
  ov = findOverlaps(query = loopEnh,
                    subject = promoters(genes[genes$SYMBOL %in% input2$SYMBOL],
                                        upstream = 10000,
                                        downstream = 5000))
  enhK27 = unique(loopEnh[-(queryHits(ov)),])
  proK27 = unique(loopEnh[queryHits(ov),])
  
  ov = findOverlaps(query = input2,
                    subject = repr)
  k27me3 = unique(repr[subjectHits(ov),])
  
  dat2 = list("repr" = k27me3$TC.log2FoldChange,
              "enh" = enhK27$TC.log2FoldChange,
              "pro" = proK27$TC.log2FoldChange,
              "gene" = input2$TC.log2FoldChange,
              "loop" = input2$TClfc)
  
  dat = list(dat1$repr, dat2$repr,
             dat1$enh, dat2$enh,
             dat1$pro, dat2$pro,
             dat1$gene, dat2$gene,
             dat1$loop, dat2$loop)
  
  boxplot(dat,
          col = alpha(rep(c("grey", "firebrick", "orange", "gold", "steelblue"), each = 2),
                      rep(c(0.75, 0.25), times = 5)),
          outline = F,
          # ylim = c(-1, 10),
          main = title)
  stripchart(dat,
             col = alpha(rep(c("grey", "firebrick", "orange", "gold", "steelblue"), each = 2),
                         rep(c(0.5, 0.25), times = 5)),
             vertical = T, 
             method = "jitter",
             jitter = 0.25,
             pch = 19,
             add = T)
  
  boxplot(dat,
          col = NA,
          outline = F,
          add = T)
  abline(h=0)
  
  return(list("repr" = t.test(dat1$repr, dat2$repr),
              "enh"  = t.test(dat1$enh, dat2$enh),
              "pro"  = t.test(dat1$pro, dat2$pro),
              "gene" = t.test(dat1$gene, dat2$gene),
              "loop" = t.test(dat1$loop, dat2$loop)))
}


### PLOT: Looped-DEG feature boxplots (TvC) -------
pdf("./output/loopDEGanalysis/FigS5E-loopedDEGfeatureBoxplotsTvC.pdf",
    width = 12, height = 6)

bothGained = geneLoops[geneLoops$TC.sig & geneLoops$TCdiff &
                         geneLoops$TC.log2FoldChange > 0 &
                         geneLoops$TClfc > 0,]
oppGained = geneLoops[geneLoops$TC.sig & geneLoops$TCdiff &
                        geneLoops$TC.log2FoldChange > 0 &
                        geneLoops$TClfc < 0,]
lfcBoxplotTC(input1 = bothGained, input2 = oppGained)

bothLost = geneLoops[geneLoops$TC.sig & geneLoops$TCdiff &
                       geneLoops$TC.log2FoldChange < 0 &
                       geneLoops$TClfc < 0,]
oppLost = geneLoops[geneLoops$TC.sig & geneLoops$TCdiff &
                      geneLoops$TC.log2FoldChange < 0 &
                      geneLoops$TClfc > 0,]
lfcBoxplotTC(input1 = bothLost, input2 = oppLost)

dev.off()
