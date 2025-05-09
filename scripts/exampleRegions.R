
# INITIALIZE --------------------------------------------

library(plotgardener)
library(mariner)
library(rtracklayer)


# READ IN --------------------------------------------

## .hic files --------
hicFiles = list.files("/Volumes/MCF10/CantataData/proc/celltype/hicFiles/",
                      full.names = T)
names(hicFiles) = c("A", "C", "T")
hicFiles = hicFiles[c(1, 3, 2)]

## loops -------
loops = readRDS("./output/loopProcessing/loops.rds")

## Gene loops ----
geneLoops = readRDS("./output/loopGeneOverlap/geneLoops.rds")

## TADs -------
tads = readRDS("./output/tadProcessing/tads.RDS")

## IS boundaries -------
bounds = readRDS("./output/tadProcessing/tadBounds.RDS")

## IS tracks -------
IS = readRDS("./output/tadProcessing/IStracks.RDS")

## DEG info --------
genes = readRDS("./output/rnaProcessing/genes.rds")

## CTCF signal --------
ctcfSig = list.files("./input/ctcf/signal/", pattern = "*.bw", full.names = T)
names(ctcfSig) = gsub("_CTCF_pooled_FE.bw", "", gsub("GSE98551_", "", basename(ctcfSig)))

## Diff CTCF peaks -----
dir = "~/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/data/external/miscSteinLab/ctcf/csaw/"
csawFiles = list.files("./input/ctcf/csaw/", pattern = "*txt", full.names = T)
diffCTCF = lapply(csawFiles, read.table, header = T)
diffCTCF = lapply(diffCTCF, GRanges)
names(diffCTCF) = gsub(".txt", "", gsub("CTCF_csaw_DB_results_", "",
                                    basename(csawFiles)))
diffCTCF = unlist(as(diffCTCF, "GRangesList"))
diffCTCF = reduce(diffCTCF)

## RNA signal --------
rnaSig = list.files("./input/rna/signal/", full.names = T)
names(rnaSig) = gsub("_combined.bw", "", basename(rnaSig))

## H3K27ac signal --------
k27acSig = list.files("./input/h3k27ac/signal/", full.names = T)
names(k27acSig) = gsub("_H3K27ac_pooled.Aligned.sortedByCoord.out_macs2_FE.bw", 
                       "", basename(k27acSig))

## ATAC signal -------
atacSig = list.files("./input/atac/signal/", full.names = T)
names(atacSig) = gsub(".mRp.clN.bigWig", "", basename(atacSig))

## Enhancers -----
enh = readRDS("./output/H3K27acProcessing/enhancers.rds")

## CALDER comps -----------
comps = readRDS("./output/compartments/compartments.RDS")


# FUNCTIONS --------------------------------------------

## Function for plotting loops of all clusters
loopPlotter <- function(plot, loopset){
  loopset = as_ginteractions(loopset)
  annoPixels(plot, data = loopset[loopset$cluster == "static",], type = "box", col = "black", half = "top")
  annoPixels(plot, data = loopset[loopset$cluster %in% c("up.early", "up.late", "up.mid"),], type = "box", col = "firebrick", half = "top", lwd = 2)
  annoPixels(plot, data = loopset[loopset$cluster %in% c("down.early", "down.late", "down.mid"),], type = "box", col = "firebrick", half = "top", lwd = 2)
}


# MICRO-C FEATURE PLOT ---------
## Plot region ----
plotWindow = GRanges(
  data.frame(seqnames = "chr1",
             start = 158000000,
             stop = 160000000))

## Plot settings --------
plotW = 8
plotB = 0.05

hicH = 2
sigH = 0.25
tadH = 0.5

pageB = 1
pageW = plotW + pageB*2
pageH = pageB*2 + (hicH + sigH*2 + tadH + plotB*2)*3

IScol = "#1f4297"
tadCol = "#1f4297"


## PLOT: Structure example --------
pdf("./output/exampleRegions/Fig1B_exampleRegion.pdf",
    width = pageW, height = pageH)

### Set coordinates -------
chr = as.character(seqnames(plotWindow))
wStart = start(plotWindow)
wEnd = end(plotWindow)

### Universal params -------
p = pgParams(chrom = chr,
             chromstart = wStart,
             chromend = wEnd,
             width = plotW,
             x = pageB,
             norm = "KR",
             resolution = 10000)

### Make page -------
pageCreate(width = pageW, height = pageH,
           showGuides = F)

### Plot Hi-C -------
dat = readHic(file = hicFiles["A"], params = p)
Z = quantile(dat$counts, 0.95)
hicA = plotHicRectangle(data = hicFiles["A"],
                        params = p,
                        y = pageB, 
                        height = hicH,
                        zrange = c(0, Z))

dat = readHic(file = hicFiles["T"], params = p)
Z = quantile(dat$counts, 0.95)
hicT = plotHicRectangle(data = hicFiles["T"],
                        params = p,
                        y = pageB + (hicH + sigH*2 + tadH + plotB*2), 
                        height = hicH,
                        zrange = c(0, Z))

dat = readHic(file = hicFiles["C"], params = p)
Z = quantile(dat$counts, 0.95)
hicC = plotHicRectangle(data = hicFiles["C"],
                        params = p,
                        y = pageB + (hicH + sigH*2 + tadH + plotB*2)*2, 
                        height = hicH,
                        zrange = c(0, Z))

#### Anno: legends -------
annoHeatmapLegend(plot = hicA,
                  x = pageB + plotW + pageB*.1,
                  y = pageB,
                  height = hicH*.5,
                  width = pageB*0.15)
annoHeatmapLegend(plot = hicT,
                  x = pageB + plotW + pageB*.1,
                  y = pageB + (hicH + sigH*2 + tadH + plotB*2),
                  height = hicH*.5,
                  width = pageB*0.15)
annoHeatmapLegend(plot = hicC,
                  x = pageB + plotW + pageB*.1,
                  y = pageB + (hicH + sigH*2 + tadH + plotB*2)*2,
                  height = hicH*.5,
                  width = pageB*0.15)

#### Anno: loops -------
loopPlotter(plot = hicA, loopset = loops)
loopPlotter(plot = hicT, loopset = loops)
loopPlotter(plot = hicC, loopset = loops)

### Plot CALDER comps --------
plotRanges(data = comps$MCF10A, 
           collapse = T,
           params = p, 
           x = pageB,
           y = pageB + hicH + plotB,
           height = sigH, 
           fill = colorby("V9", 
                          palette = colorRampPalette(c("blue", "white", "red"))))

plotRanges(data = comps$MCF10AT1, 
           collapse = T,
           params = p, 
           x = pageB,
           y = pageB + (hicH + sigH*2 + tadH  + plotB*2) + hicH + plotB,
           height = sigH, 
           fill = colorby("V9", 
                          palette = colorRampPalette(c("blue", "white", "red"))))

plotRanges(data = comps$MCF10CA1a, 
           collapse = T,
           params = p, 
           x = pageB,
           y = pageB + (hicH + sigH*2 + tadH + plotB*2)*2 + hicH + plotB,
           height = sigH, 
           fill = colorby("V9", 
                          palette = colorRampPalette(c("blue", "white", "red"))))

#### Anno: IS bounds -------
plotRanges(data = bounds, 
           collapse = T,
           params = p, 
           y = pageB + hicH + plotB*2 + sigH,
           height = sigH, 
           fill = "darkgrey")
plotRanges(data = bounds[bounds$sig == T,], 
           collapse = T,
           params = p, 
           y = pageB + hicH + plotB*2 + sigH,
           height = sigH, 
           fill = "firebrick")

plotRanges(data = bounds, 
           collapse = T,
           params = p, 
           y = pageB + (hicH + sigH*2 + tadH + plotB*2) + hicH + plotB*2 + sigH,
           height = sigH, 
           fill = "darkgrey")
plotRanges(data = bounds[bounds$sig == T,], 
           collapse = T,
           params = p, 
           y = pageB + (hicH + sigH*2 + tadH + plotB*2) + hicH + plotB*2 + sigH,
           height = sigH, 
           fill = "firebrick")

plotRanges(data = bounds, 
           collapse = T,
           params = p, 
           y = pageB + (hicH + sigH*2 + tadH + plotB*2)*2 + hicH + plotB*2 + sigH,
           height = sigH, 
           fill = "darkgrey")
plotRanges(data = bounds[bounds$sig == T,], 
           collapse = T,
           params = p, 
           y = pageB + (hicH + sigH*2 + tadH + plotB*2)*2 + hicH + plotB*2 + sigH,
           height = sigH, 
           fill = "firebrick")

### Plot IS -------
sigRange = c(-2, 0.5)
plotSignal(data = IS$MCF10_A[is.finite(IS$MCF10_A$score),],
           params = p,
           y = pageB + hicH + plotB*2 + sigH,
           height = sigH,
           range = sigRange,
           linecolor = IScol)
plotSignal(data = IS$MCF10_T[is.finite(IS$MCF10_T$score),],
           params = p,
           y = pageB + (hicH + sigH*2 + tadH + plotB*2) + hicH + plotB*2 + sigH,
           height = sigH,
           range = sigRange,
           linecolor = IScol)
plotSignal(data = IS$MCF10_C[is.finite(IS$MCF10_C$score),],
           params = p,
           y = pageB + (hicH + sigH*2 + tadH + plotB*2)*2 + hicH + plotB*2 + sigH,
           height = sigH,
           range = sigRange,
           linecolor = IScol)

### Plot TADs -------
tads$width = pairdist(tads)
plotPairsArches(data = tads[tads$A],
                params = p,
                y = pageB + hicH + plotB*2 + sigH*2,
                height = tadH,
                flip = T,
                alpha = 1,
                fill = tadCol,
                archHeight = "width")

plotPairsArches(data = tads[tads$T],
                params = p,
                y = pageB + (hicH + sigH*2 + tadH + plotB*2) + hicH + plotB*2 + sigH*2,
                height = tadH,
                flip = T,
                alpha = 1,
                fill = tadCol,
                archHeight = "width")

plotPairsArches(data = tads[tads$C],
                params = p,
                y = pageB + (hicH + sigH*2 + tadH + plotB*2)*2 + hicH + plotB*2 + sigH*2,
                height = tadH,
                flip = T,
                alpha = 1,
                fill = tadCol,
                archHeight = "width")

### Anno: Coords --------
annoGenomeLabel(plot = hicC,
                x = pageB,
                y = pageB + (hicH + sigH*2 + tadH + plotB*2)*3,
                at = seq(ceiling(wStart/2e5)*2e5, 
                         floor(wEnd/2e5)*2e5, 
                         by=2e5), 
                scale = "Mb")

### Anno: Labels --------
plotText(label = "MCF10A",
         fontface = "bold",
         x = pageB + .1,
         y = pageB + .1,
         just = c("top", "left"))
plotText(label = "MCF10AT1",
         fontface = "bold",
         x = pageB + .1,
         y = pageB + (hicH + sigH*2 + tadH + plotB*2) + .1,
         just = c("top", "left"))
plotText(label = "MCF10CA1a",
         fontface = "bold",
         x = pageB + .1,
         y = pageB + (hicH + sigH*2 + tadH + plotB*2)*2 + .1,
         just = c("top", "left"))

plotText(label = "Comp",
         x = pageB - .1,
         y = pageB + hicH + plotB + sigH*0.5,
         fontcolor = IScol,
         just = c("center", "right"))

plotText(label = "IS",
         x = pageB - .1,
         y = pageB + hicH + plotB*2 + sigH*1.5,
         fontcolor = IScol,
         just = c("center", "right"))

plotText(label = "TADs",
         x = pageB - .1,
         y = pageB + hicH + plotB*2 + sigH*2 + tadH*0.5,
         fontcolor = tadCol,
         just = c("center", "right"))

dev.off()



# FULL PLOT --------------------------------------------

## Find regions of interest --------
goi = c("WNT5A", "COL12A1",
        "SPRY1", "SCNN1G")
loi = geneLoops[geneLoops$SYMBOL %in% goi,]
loi = loi[c(1, 10, 20, 23),]

## Find number of unique plot windows
spans = GRanges(
  data.frame(seqnames = seqnames(anchors(loi, "first")),
             start = start(anchors(loi, "first")),
             stop = end(anchors(loi, "second"))))
plotWindows = spans + width(spans) - 1
plotWindows = reduce(plotWindows)


## Plot settings --------
plotW = 8
plotB = 0.05

hicH = 2
sigH = 0.25
tadH = 1
geneH = 1

pageB = 1
pageW = plotW + pageB*2
pageH = hicH*3 + (sigH*3)*6 + tadH + geneH + pageB*2 + plotB*6
Tbuff = (hicH + sigH*2 + plotB*3)
Cbuff = (hicH + sigH*2 + plotB*3)*2
bBuff = (hicH + sigH*2 + plotB*3)*3

IScol = "#0C2C84"
tadCol = "#0C2C84"
ctcfCol = "#225EA8"
k27col = "#1D91C0"
atacCol = "#2cc79e"
peakCol = "darkgrey"
diffCol = "firebrick"
rnaCol = "#FFBD16"

geneCols = data.frame("gene" = genes$SYMBOL)
geneCols$color = "grey50"
geneCols$color[genes$AT.sig == T] = "black"

## PLOT: Gene examples --------
pdf("./output/exampleRegions/Fig2HI-4DE_exampleRegions.pdf",
    width = pageW, height = pageH)

for(n in 1:length(plotWindows)){ 
  print(paste0("Starting page ", n, " out of ", length(plotWindows)))
  
  ### Set coordinates -------
  chr = as.character(seqnames(plotWindows[n,]))
  wStart = start(plotWindows[n,])
  wEnd = end(plotWindows[n,])
  
  ### Universal params -------
  p = pgParams(chrom = chr,
               chromstart = wStart,
               chromend = wEnd,
               width = plotW,
               x = pageB,
               norm = "SCALE",
               resolution = 5000)
  
  ### Make page -------
  pageCreate(width = pageW, height = pageH,
             showGuides = F)
  
  ### Plot Hi-C -------
  dat = readHic(file = hicFiles["A"], params = p)
  Z = quantile(dat$counts, 0.95)
  hicA = plotHicRectangle(data = hicFiles["A"],
                          params = p,
                          y = pageB, 
                          height = hicH,
                          zrange = c(0, Z))
  
  dat = readHic(file = hicFiles["T"], params = p)
  Z = quantile(dat$counts, 0.95)
  hicT = plotHicRectangle(data = hicFiles["T"],
                          params = p,
                          y = pageB + Tbuff, 
                          height = hicH,
                          zrange = c(0, Z))
  
  dat = readHic(file = hicFiles["C"], params = p)
  Z = quantile(dat$counts, 0.95)
  hicC = plotHicRectangle(data = hicFiles["C"],
                          params = p,
                          y = pageB + Cbuff, 
                          height = hicH,
                          zrange = c(0, Z))
  
  ### Anno: legends -------
  annoHeatmapLegend(plot = hicA,
                    x = pageB + plotW + pageB*.1,
                    y = pageB,
                    height = hicH*.5,
                    width = pageB*0.15)
  annoHeatmapLegend(plot = hicT,
                    x = pageB + plotW + pageB*.1,
                    y = pageB + Tbuff,
                    height = hicH*.5,
                    width = pageB*0.15)
  annoHeatmapLegend(plot = hicC,
                    x = pageB + plotW + pageB*.1,
                    y = pageB + Cbuff,
                    height = hicH*.5,
                    width = pageB*0.15)
  
  ### Anno: loops -------
  loopPlotter(plot = hicA, loopset = loops)
  loopPlotter(plot = hicT, loopset = loops)
  loopPlotter(plot = hicC, loopset = loops)
  
  ### Plot CALDER comps --------
  plotRanges(data = comps$MCF10A, 
             collapse = T,
             params = p, 
             x = pageB,
             y = pageB + hicH + plotB,
             height = sigH, 
             fill = colorby("V9", 
                            palette = colorRampPalette(c("blue", "white", "red"))))
  
  plotRanges(data = comps$MCF10AT1, 
             collapse = T,
             params = p, 
             x = pageB,
             y = pageB + Tbuff + hicH + plotB,
             height = sigH, 
             fill = colorby("V9", 
                            palette = colorRampPalette(c("blue", "white", "red"))))
  
  plotRanges(data = comps$MCF10CA1a, 
             collapse = T,
             params = p, 
             x = pageB,
             y = pageB + Cbuff + hicH + plotB,
             height = sigH, 
             fill = colorby("V9", 
                            palette = colorRampPalette(c("blue", "white", "red"))))
  
  ### Anno: IS bounds -------
  plotRanges(data = bounds[bounds$sig == F,], 
             collapse = T,
             params = p, 
             y = pageB + hicH + plotB*2 + sigH,
             height = sigH, 
             fill = alpha(peakCol, 0.5))
  plotRanges(data = bounds[bounds$sig == T,], 
             collapse = T,
             params = p, 
             y = pageB + hicH + plotB*2 + sigH,
             height = sigH, 
             fill = alpha(diffCol, 0.3))
  
  plotRanges(data = bounds[bounds$sig == F,], 
             collapse = T,
             params = p, 
             y = pageB + Tbuff + hicH + plotB*2 + sigH,
             height = sigH, 
             fill = alpha(peakCol, 0.5))
  plotRanges(data = bounds[bounds$sig == T,], 
             collapse = T,
             params = p, 
             y = pageB + Tbuff + hicH + plotB*2 + sigH,
             height = sigH, 
             fill = alpha(diffCol, 0.3))
  
  plotRanges(data = bounds[bounds$sig == F,], 
             collapse = T,
             params = p, 
             y = pageB + Cbuff + hicH + plotB*2 + sigH,
             height = sigH, 
             fill = alpha(peakCol, 0.5))
  plotRanges(data = bounds[bounds$sig == T,], 
             collapse = T,
             params = p, 
             y = pageB + Cbuff + hicH + plotB*2 + sigH,
             height = sigH, 
             fill = alpha(diffCol, 0.3))
  
  ### Plot IS -------
  sigRange = c(-2, 0.5)
  plotSignal(data = IS$MCF10_A[is.finite(IS$MCF10_A$score),],
             params = p,
             y = pageB + hicH + plotB*2 + sigH,
             height = sigH,
             range = sigRange,
             linecolor = IScol)
  plotSignal(data = IS$MCF10_T[is.finite(IS$MCF10_T$score),],
             params = p,
             y = pageB + Tbuff + hicH + plotB*2 + sigH,
             height = sigH,
             range = sigRange,
             linecolor = IScol)
  plotSignal(data = IS$MCF10_C[is.finite(IS$MCF10_C$score),],
             params = p,
             y = pageB + Cbuff + hicH + plotB*2 + sigH,
             height = sigH,
             range = sigRange,
             linecolor = IScol)
  
  tempY = pageB + bBuff - plotB
  
  ### Plot TADs -------
  tads$width = pairdist(tads)
  plotPairsArches(data = tads,
                  params = p,
                  y = tempY,
                  height = tadH,
                  flip = T,
                  alpha = 1,
                  fill = tadCol,
                  archHeight = "width")
  
  tempY = tempY + tadH
  
  ### Anno: Diff CTCF peaks ----------
  plotRanges(data = resize(GRanges(diffCTCF), width = 5000, fix = "center"),
             collapse = T,
             params = p,
             y = tempY,
             height = sigH*3,
             fill = alpha(diffCol, 0.3))
  
  ### Plot CTCF -------
  tmp = lapply(ctcfSig, function(c){readBigwig(c, params = p)})
  max = max(unlist(lapply(tmp, function(c){
    quantile(c$score, 1)
  })))
  
  ctcfA = plotSignal(data = ctcfSig[[1]], 
                     params = p,
                     height = sigH,
                     y = tempY, 
                     range = c(0, max),
                     linecolor = ctcfCol,
                     fill = ctcfCol)
  ctcfT = plotSignal(data = ctcfSig[[2]], 
                     params = p,
                     height = sigH,
                     y = tempY + (sigH), # + plotB 
                     range = c(0, max),
                     linecolor = ctcfCol,
                     fill = ctcfCol)
  ctcfC = plotSignal(data = ctcfSig[[3]], 
                     params = p,
                     height = sigH,
                     y = tempY + (sigH)*2, # + plotB
                     range = c(0, max),
                     linecolor = ctcfCol,
                     fill = ctcfCol)
  
  tempY = tempY + (sigH)*3 # + plotB
  
  ### Anno: Enhancers -------
  plotRanges(data = resize(enh[!enh$AT.sig & !enh$TC.sig & !enh$AC.sig], 
                           width = 5000, fix = "center"),
             collapse = T,
             params = p,
             y = tempY,
             height = sigH*6,
             fill = alpha(peakCol, 0.2))
  plotRanges(data = resize(enh[enh$AT.sig | enh$TC.sig | enh$AC.sig], 
                           width = 5000, fix = "center"),
             collapse = T,
             params = p,
             y = tempY,
             height = sigH*6,
             fill = alpha(diffCol, 0.3))
  
  ### Plot H3K27ac -------
  tmp = lapply(k27acSig, function(c){readBigwig(c, params = p)})
  max = max(unlist(lapply(tmp, function(c){
    quantile(c$score, 1)
  })))
  
  k27A = plotSignal(data = k27acSig[[1]], 
                    params = p,
                    height = sigH,
                    y = tempY, 
                    range = c(0, max),
                    linecolor = k27col,
                    fill = k27col)
  k27T = plotSignal(data = k27acSig[[2]], 
                    params = p,
                    height = sigH,
                    y = tempY + (sigH), # + plotB  
                    range = c(0, max),
                    linecolor = k27col,
                    fill = k27col)
  k27C = plotSignal(data = k27acSig[[3]], 
                    params = p,
                    height = sigH,
                    y = tempY + (sigH)*2, # + plotB  
                    range = c(0, max),
                    linecolor = k27col,
                    fill = k27col)
  
  tempY = tempY + (sigH)*3 # + plotB
  
  ### Plot ATAC -------
  tmp = lapply(atacSig, function(c){readBigwig(c, params = p)})
  max = max(unlist(lapply(tmp, function(c){
    quantile(c$score, 1)
  })))
  
  atacA = plotSignal(data = atacSig[[1]], 
                     params = p,
                     height = sigH,
                     y = tempY, 
                     range = c(0, max),
                     linecolor = atacCol,
                     fill = atacCol)
  atacT = plotSignal(data = atacSig[[2]], 
                     params = p,
                     height = sigH,
                     y = tempY + (sigH), # + plotB  
                     range = c(0, max),
                     linecolor = atacCol,
                     fill = atacCol)
  atacC = plotSignal(data = atacSig[[3]], 
                     params = p,
                     height = sigH,
                     y = tempY + (sigH)*2, # + plotB  
                     range = c(0, max),
                     linecolor = atacCol,
                     fill = atacCol)
  
  tempY = tempY + (sigH)*3 # + plotB 
  
  ### Plot RNA -------
  tmp = lapply(rnaSig, function(c){readBigwig(c, params = p)})
  max = max(unlist(lapply(tmp, function(c){
    quantile(c$score, 0.99)
  })))
  
  rnaA = plotSignal(data = rnaSig[[1]], 
                    params = p,
                    height = sigH,
                    y = tempY, 
                    range = c(0, max),
                    linecolor = rnaCol,
                    fill = rnaCol)
  rnaT = plotSignal(data = rnaSig[[2]], 
                    params = p,
                    height = sigH,
                    y = tempY + (sigH), # + plotB  
                    range = c(0, max),
                    linecolor = rnaCol,
                    fill = rnaCol)
  rnaC = plotSignal(data = rnaSig[[3]], 
                    params = p,
                    height = sigH,
                    y = tempY + (sigH)*2, # + plotB  
                    range = c(0, max),
                    linecolor = rnaCol,
                    fill = rnaCol)
  tempY = tempY + (sigH)*3 # + plotB 
  
  ### Plot genes --------
  plotG = plotGenes(params = p,
                    y = tempY, 
                    height = geneH,
                    geneHighlights = geneCols)
  
  #### Anno: Coords --------
  annoGenomeLabel(plot = plotG,
                  x = pageB,
                  y = tempY + geneH,
                  at = seq(ceiling(wStart/2e5)*2e5, 
                           floor(wEnd/2e5)*2e5, 
                           by=2e5), 
                  scale = "Mb")
  tryCatch(annoGenomeLabel(plot = plotG,
                           x = pageB,
                           y = tempY + geneH,
                           at = seq(ceiling(wStart/1e6)*1e6, 
                                    floor(wEnd/1e6)*1e6, 
                                    by=1e6), 
                           tcl = 1,
                           fontsize = 0), 
           error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  ### Anno: Labels --------
  plotText(label = "MCF10A",
           fontface = "bold",
           x = pageB + .1,
           y = pageB + .1,
           just = c("top", "left"))
  plotText(label = "MCF10AT1",
           fontface = "bold",
           x = pageB + .1,
           y = pageB + Tbuff + .1,
           just = c("top", "left"))
  plotText(label = "MCF10CA1a",
           fontface = "bold",
           x = pageB + .1,
           y = pageB + Cbuff + .1,
           just = c("top", "left"))
  
  plotText(label = "Compartment",
           x = pageB - .1,
           y = pageB + hicH + sigH*.5 + plotB,
           fontcolor = IScol,
           just = c("center", "right"))
  plotText(label = "IS",
           x = pageB - .1,
           y = pageB + hicH + sigH*1.5 + plotB*2,
           fontcolor = IScol,
           just = c("center", "right"))
  tempY = pageB + bBuff
  
  plotText(label = "TADs",
           x = pageB - .1,
           y = tempY + tadH*0.5,
           fontcolor = tadCol,
           just = c("center", "right"))
  tempY = tempY + tadH
  
  plotText(label = "CTCF",
           x = pageB - .1,
           y = tempY + (sigH) + sigH*0.5,
           fontcolor = ctcfCol,
           just = c("center", "right"))
  tempY = tempY + (sigH)*3
  
  plotText(label = "H3K27ac",
           x = pageB - .1,
           y = tempY + (sigH) + sigH*0.5,
           fontcolor = k27col,
           just = c("center", "right"))
  tempY = tempY + (sigH)*3
  
  plotText(label = "ATAC",
           x = pageB - .1,
           y = tempY + (sigH) + sigH*0.5,
           fontcolor = atacCol,
           just = c("center", "right"))
  tempY = tempY + (sigH)*3
  
  plotText(label = "RNA",
           x = pageB - .1,
           y = tempY + (sigH) + sigH*0.5,
           fontcolor = rnaCol,
           just = c("center", "right"))
  tempY = tempY + (sigH)*3
  
  print(paste0("Finished page ", n, " out of ", length(plotWindows)))
}

dev.off()
