
# INTIALIZE -----------------------
library(GenomicRanges)
library(mariner)
library(matrixStats)
library(rtracklayer)


# READ IN -----------------------

## Read in loops -------
loops = as_ginteractions(readRDS("./output/loopProcessing/loops.rds"))

## Read in CTCF peaks ------
peakFiles = list.files("./input/ctcf/peaks",
                       pattern = "*.narrowPeak",
                       full.names = T)
names(peakFiles) = gsub("GSE98551_", 
                        "",
                        gsub("_CTCF_pooled_peaks_passIDR.05.narrowPeak",
                             "",
                             basename(peakFiles)))

peaks = lapply(peakFiles, import, format = "narrowPeak")
allPeaks = unlist(as(peaks, "GRangesList"))
ctcf = reduce(allPeaks)


## Read in diff CTCF peaks ------
ACctcf = read.table("./input/ctcf/csaw/CTCF_csaw_DB_results_AvsC.txt",
                    sep = "\t", header = T)
ATctcf = read.table("./input/ctcf/csaw/CTCF_csaw_DB_results_AvsT.txt",
                    sep = "\t", header = T)
TCctcf = read.table("./input/ctcf/csaw/CTCF_csaw_DB_results_TvsC.txt",
                    sep = "\t", header = T)


# CTCF OVERLAP -----------------------

## Anchor overlap: 77.6% of anchors overlap CTCF peaks
ov = countOverlaps(query = regions(loops),
                   subject = ctcf, 
                   maxgap = 10000)
sum(ov >= 1)/length(ov)

## Loop overlap: 72.4% of loops have double overlap, 95.0% have single overlap
Lov = countOverlaps(query = loops,
                    subject = ctcf, 
                    maxgap = 10000,
                    use.region = "first")
Rov = countOverlaps(query = loops,
                    subject = ctcf, 
                    maxgap = 10000,
                    use.region = "second")

both = which(Rov > 0 & Lov > 0)
either = which((Rov > 0 & Lov == 0) | (Rov == 0 & Lov > 0))
either = which(Rov > 0 | Lov > 0)
neither = which(Rov == 0 & Lov == 0)

length(both)/length(loops)
length(either)/length(loops)
length(neither)/length(loops)


## PLOT: CTCF overlap by loop strength -----------------------
pdf("./output/loopCTCFoverlap/FigS2E-CTCFbyLoopStrength.pdf",
    width = 6,
    height = 6)
cuts = cut(loops$maxCount, 
           unique(quantile(loops$maxCount, seq(0, 1, 0.05))), 
           include.lowest = T)
barplot(table(cuts)/table(cuts), space = 0, las = 2, names.arg = NA,
        main = "CTCF Occupancy by Loop Strength")
barplot(table(cuts[either])/table(cuts), col = "firebrick1", space = 0, border = NA, add = T, axes = F, names.arg = NA)
barplot(table(cuts[both])/table(cuts), col = "firebrick4", space = 0, border = NA, add = T, axes = F, names.arg = NA)
barplot(table(cuts)/table(cuts), space = 0, col = NA, add = T, axes = F, names.arg = NA)
dev.off()


## CTCF-dependent trends ----------
diffLoops = loops[loops$cluster != "static",]
Lov = countOverlaps(query = diffLoops,
                    subject = ctcf, 
                    maxgap = 10000,
                    use.region = "first")
Rov = countOverlaps(query = diffLoops,
                    subject = ctcf, 
                    maxgap = 10000,
                    use.region = "second")

both = which(Rov > 0 & Lov > 0)
either = which((Rov > 0 & Lov == 0) | (Rov == 0 & Lov > 0))
neither = which(Rov == 0 & Lov == 0)


## PLOT: Loop metrics by CTCF overlap -----------------------
pdf("./output/loopCTCFoverlap/FigS2F-loopMetricsByCTCF.pdf",
    width = 12,
    height = 6)
par(mfrow=c(1,3))
dat = list("neither" = diffLoops$maxCount[neither],
           "either" = diffLoops$maxCount[either],
           "both" = diffLoops$maxCount[both])
boxplot(dat, main = "Max loop counts", 
        col = c("grey", "firebrick1", "firebrick4"))

dat = list("neither" = abs(diffLoops$maxLFC)[neither],
           "either" = abs(diffLoops$maxLFC)[either],
           "both" = abs(diffLoops$maxLFC)[both])
boxplot(dat, main = "Max abs LFC", 
        col = c("grey", "firebrick1", "firebrick4"))

dat = list("neither" = diffLoops$span[neither],
           "either" = diffLoops$span[either],
           "both" = diffLoops$span[both])
boxplot(dat, main = "Loop length", outline=F, 
        col = c("grey", "firebrick1", "firebrick4"))
dev.off()


# DIFF LOOP CTCF ----------

## Identify CTCF sites that change in the same or opposite direction as differential loops
ctcfClass <- function(csawOutput, loopSet, LFCcol){
  ov = findOverlaps(query = loopSet,
                    subject = GRanges(csawOutput[csawOutput$logFC.up > csawOutput$logFC.down,]),
                    maxgap = 10000)
  hits1 = loopSet[queryHits(ov),]
  ov = findOverlaps(query = loopSet,
                    subject = GRanges(csawOutput[csawOutput$logFC.up < csawOutput$logFC.down,]),
                    maxgap = 10000)
  hits2 = loopSet[queryHits(ov),]
  
  supported = unique(c(hits1$name[mcols(hits1)[,LFCcol] > 0],
                       hits2$name[mcols(hits2)[,LFCcol] < 0]))
  opposite = unique(c(hits1$name[mcols(hits1)[,LFCcol] < 0],
                      hits2$name[mcols(hits2)[,LFCcol] > 0]))
  
  return(list(supported, opposite))
}

CTCFclasses = Map(ctcfClass,
                  csawOutput = list(ATctcf, ACctcf, TCctcf),
                  loopSet = list(diffLoops[diffLoops$ATdiff],
                                 diffLoops[diffLoops$ACdiff],
                                 diffLoops[diffLoops$TCdiff]),
                  LFCcol = list("ATlfc", "AClfc", "TClfc"))

supported = unique(c(CTCFclasses[[1]][[1]],
                     CTCFclasses[[2]][[1]],
                     CTCFclasses[[3]][[1]]))
opposite = unique(c(CTCFclasses[[1]][[2]],
                    CTCFclasses[[2]][[2]],
                    CTCFclasses[[3]][[2]]))
set = diffLoops$name

## PLOT: CTCF status at differential loops -----------
dat = c(length(supported), # 277
        sum(!opposite %in% supported), # 137
        sum(!(set %in% c(supported, opposite))))
percents = signif(dat/sum(dat)*100, digits = 2)

pdf("./output/loopCTCFoverlap/FigS2I-CTCFdiffLoopPie.pdf", width = 6, height = 6)

par(mfrow=c(1,1))
pie(c(length(supported), # 277
      sum(!opposite %in% supported), # 137
      sum(!(set %in% c(supported, opposite)))), # 1171
    main = "CTCF at changing loops",
    border = NA,
    col = c("firebrick", "darkgrey", "lightgrey"),
    labels = c(paste0("same (", percents[[1]], "%)"), 
               paste0("opposite (", percents[[2]], "%)"),
               paste0("static CTCF (", percents[[3]], "%)")))

dev.off()
