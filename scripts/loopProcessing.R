
# INITIALIZE  -----------------------------------------------------------------------------

library(mariner)
library(DESeq2)
library(InteractionSet)
library(scales)



# READ IN     -----------------------------------------------------------------------------

## Loops
# # Comprehensive
# celldir = "./input/loops/celltype/comprehensive//"
# megadir = "./input/loops/mega/comprehensive/mega_loops.txt"

# Conservative
celldir = "./input/loops/celltype/conservative/"
megadir = "./input/loops/mega/conservative/mega_loops.txt"

loopFiles = list.files(celldir,  
                       pattern = "*loops.txt", 
                       full.names = T)
loopFiles = c(loopFiles, megadir)
names(loopFiles) = c("MCF10_A", "MCF10_C", "MCF10_T", "mega")
loopFiles = loopFiles[c("MCF10_A", "MCF10_T", "MCF10_C", "mega")]

loopSets = lapply(loopFiles, read.table, header = T)
loopsGR = lapply(loopSets, as_ginteractions)


## Tech rep hic files (for counts)
trHicPaths = list.files("/Volumes/MCF10/CantataData/proc/techrep/hicFiles/",
                        full.names = T)
trHicPaths = trHicPaths[c(1:8, 17:24, 9:16)]



# RUN     -----------------------------------------------------------------------------

## MERGING LOOPS ---------
## Bin to 10kb
loops10k = lapply(loopsGR, assignToBins, binSize = 10000)

## Merge loops
mergedLoops <- 
  mergePairs(x = loops10k, 
             radius = 20000,
             method = "max",
             pos = "center")

## Order + rebin accordingly
mergedLoops = swapAnchors(mergedLoops)
mergedLoops = sort(mergedLoops)
mergedLoops = assignToBins(mergedLoops, binSize = 10000, 
                           pos1="center", pos2="center")


## Add loop call table
dat = mariner:::.makeSrcMatrix(mergedLoops)
mcols(mergedLoops) = dat

## Add in metadata columns
mergedLoops = aggMetadata(mergedLoops, 
                          columns = "APScoreAvg", 
                          funs = "mean")


## PULL COUNTS ---------
# Read in counts from hic files; runs in ~6m
counts = pullHicPixels(x = mergedLoops, 
                       files = trHicPaths, 
                       binSize = 10000, 
                       norm = "NONE")


## DIFFERENTIAL ANALYSIS ---------
# Pull count matrix for DESeq
countMatrix=counts(counts)
countMatrix=suppressWarnings(apply(countMatrix, 2, as.numeric))
rownames(countMatrix) = paste0("loop_", 1:nrow(countMatrix))

# Make DNA copy matrix
loops = as.data.frame(mergedLoops)
dnacopy <- matrix(1, nrow=nrow(countMatrix), ncol=ncol(countMatrix))

dnacopy[loops$seqnames1 == "chr1" & loops$start1 >= 143200000,  1:8] <- 1.5
dnacopy[loops$seqnames1 == "chr8" & loops$start1 >= 54600000,  1:8] <- 1.4
dnacopy[loops$seqnames1 == "chr9" & loops$end1 <= 39000000,  1:8] <- 0.75
dnacopy[loops$seqnames1 == "chr20", 9:16] <- 1.4
dnacopy[loops$seqnames1 == "chr3" & loops$start1 >= 94000000,  17:24] <- 1.4
dnacopy[loops$seqnames1 == "chr9" & loops$end1 <= 39000000,  17:24] <- 1.3
dnacopy[loops$seqnames1 == "chr9" & loops$start1 >= 65000000,  17:24] <- 1.5
dnacopy[loops$seqnames1 == "chr10", 17:24] <- 1.25


# Create coldata
colData <- data.frame(cell = factor(rep(c("A", "T", "C"), each = 8)),
                      brep = factor(rep(rep(c(1, 2), each = 4), times = 3)),
                      trep = factor(rep(c(1:4), times = 6)))
rownames(colData) = colnames(countMatrix)

# Create DESeq Data Set from matrix
dds <- DESeqDataSetFromMatrix(countMatrix, 
                              colData=colData, 
                              design =~ trep + brep + cell)

# Account for differences in DNA copy #
dds <- estimateSizeFactors(dds, normMatrix=dnacopy)

# Run DESeq
dds <- DESeq(dds, 
             test="LRT", 
             full= ~ trep + brep + cell, 
             reduced = ~ trep + brep)

# Transform the dds (rlog, vst or ntd)
dds.trans <- vst(dds)

# Run results and build list
contrasts = list("resAT" = c("cell","T","A"), 
                 "resAC" = c("cell","C","A"), 
                 "resTC" = c("cell","C","T"))

resList = lapply(contrasts, function(comp){
  results(dds, contrast=comp, independentFiltering = F)
})

# Make list of all significant loops
L = log(1.5, base=2)
p = 0.025

sigList = lapply(resList, function(res){
  rownames(res)[which(res$padj <= p &
                        abs(res$log2FoldChange) >= L)]
})

allSig = unique(unlist(sigList))
length(allSig)


# Add additional loop info to DF
# Loop IDs
loops$name = rownames(resList$resAT)

# Loop spans
loops$span = loops$end2 - loops$start1

# Average + max raw counts
loops$avgCount = rowMeans(countMatrix)
loops$maxCount = rowMaxs(countMatrix)

# DESeq info
loops$ATlfc = resList$resAT$log2FoldChange
loops$AClfc = resList$resAC$log2FoldChange
loops$TClfc = resList$resTC$log2FoldChange
loops$padj = resList$resAT$padj
loops$ATdiff = loops$name %in% sigList$resAT
loops$ACdiff = loops$name %in% sigList$resAC
loops$TCdiff = loops$name %in% sigList$resTC

# Summary DESeq info
# Collect all LFC values
lfc = data.frame(row.names=rownames(resList[[1]]))

for (i in 1:length(resList)){
  newCol = resList[[i]][,"log2FoldChange"]
  lfc = cbind(lfc, newCol)
}

# Set NA's to 0
lfc[is.na(lfc)] = 0

# Find the comparison of max LFC, and the value
maxLFCcomp = max.col(abs(lfc), ties.method = "first")
maxLFC = lfc[cbind(1:nrow(lfc), maxLFCcomp)]

key = c("AT", "AC", "TC")
maxLFCcomp = key[maxLFCcomp]

# Add them to the loop DF
loops = cbind(loops,
              countMatrix,
              maxLFC,
              maxLFCcomp)

# Subset differential loops
diffLoops = loops[loops$name %in% allSig,]


## PLOT: P-value hist -------
dim(loops)
dim(diffLoops)

pdf(file="./output/loopProcessing/pvalue-hist.pdf", 
    width=6, height=6)
hist(loops$padj, 
     breaks=seq(0,1,by=0.01),
     col = rep(c("red", "grey"), 
               times = c(sum(seq(0,1,by=0.01) <= p),
                         sum(seq(0,1,by=0.01) > p))))
dev.off()


## PLOT: Loop length dist -------
breakSize = 50000

pdf(file="./output/loopProcessing/loopLengths.pdf", 
    width=6, height=6)
hist(loops$span[loops$span <= 2e6],
     breaks = seq(0, 2e6, by = breakSize),
     main = "",
     xlab = "Loop lengths")
hist(diffLoops$span[diffLoops$span <= 2e6],
     breaks = seq(0, 2e6, by = breakSize),
     col = "grey25",
     add = T)
legend("topright",
       legend = c("all loops", "differential loops"),
       fill = c("grey", "grey25"))
dev.off()



## PLOT: Loop PCA -------
sample.pal = c("steelblue", "orange", "firebrick")

pdf(file="./output/loopProcessing/PCAplot.pdf", 
    width=6, height=6)

# Plot PCA of transformed counts
pc = prcomp(t(assay(dds.trans)))
plot(pc$x[,1], pc$x[,2], 
     pch=19, cex=2.5, 
     col=rep(sample.pal, each=8), las=1,
     xlab=paste0("PC1: ", round((pc$sdev[1]^2/sum(pc$sdev^2))*100, digits=0),"% variance"),
     ylab=paste0("PC2: ", round((pc$sdev[2]^2/sum(pc$sdev^2))*100, digits=0),"% variance"), tck=.01,
     main="Loop Count PCA")
text(pc$x[,1], pc$x[,2], 
     label=rep(rep(c(1,2), each=4), times=3), 
     col="white", font=2)

dev.off()


## PLOT: Diff loop overlaps  -----
library(UpSetR)

pdf("./output/loopProcessing/diffLoopOverlaps.pdf",
    width = 8, height = 8)
dat = diffLoops[,c("MCF10_A", "MCF10_T", "MCF10_C")]
dat = as.data.frame(do.call(cbind, lapply(dat, as.numeric)))
upset(dat, 
      sets.x.label = "Total diff loop #", 
      mainbar.y.label = "Diff loops called")
dev.off()




## Read Info -----------
pdf(file="./output/loopProcessing/readInfo.pdf", 
    width=8, height=8)
par(mfrow=c(2,2))

### PLOT 1: Count boxplots for each sample (gene avg) ----
boxplot(countMatrix, 
        outline=F, 
        xaxt='n', main="Counts per loop", 
        col=rep(sample.pal, each=8))
axis(side=1, las=2, at=1:24,
     labels=gsub("_contact_map.hic", "", colnames(countMatrix)))

### PLOT 2: Total summed counts for each sample ----
bp = barplot(colSums(countMatrix),
             xaxt='n', main="Summed counts per sample", 
             col=rep(sample.pal, each=8))
axis(side=1,las=2,at=bp,
     labels=gsub("_contact_map.hic", "", colnames(countMatrix)))

### PLOT 3: Histogram of total counts per loop ----
hist(rowSums(countMatrix), 
     main="Summed counts per loop", breaks=100)
rug(rowSums(countMatrix))

### PLOT 4: Histogram of total counts per gene for a single sample ----
tmp=cbind(rowSums(countMatrix[,1:8]),
          rowSums(countMatrix[,8:16]),
          rowSums(countMatrix[,17:24]))
bp = barplot(colSums(tmp), 
             main="Summed counts per loop per sample", 
             col=sample.pal)
axis(side=1,las=1,at=bp,
     labels=c("A", "T", "C"))

### PAGE 2: Plot counts for each chromosome (to check for processing errors) ----
par(mfrow=c(4,6))
par(mar=rep(2, times=4))
for (chr in c(paste0("chr", c(1:22, "X")))){
  tmp = loops[loops$seqnames1 == chr, 27:50]
  tmp=suppressWarnings(apply(tmp, 2, as.numeric))
  barplot(colSums(tmp)/colSums(countMatrix), 
          col=rep(sample.pal, each=8), main=chr)
}
par(mar=c(5,4,4,2))
par(mfrow=c(1,1))

dev.off()


## PLOT: Check diff loops by chromosome  -----
pdf("./output/loopProcessing/loopsByChrom.pdf",
    width = 10, height = 6)
par(mfrow=c(1,2))
barplot(table(loops$seqnames1),
        main = "Loops called by chromosome",
        las = 2)
barplot(table(diffLoops$maxLFCcomp, diffLoops$seqnames1),
        col = c("white", "darkgrey", "grey40"),
        main = "Differential loops by chromosome",
        las = 2)
legend("topright",
       legend = c("A vs C", "A vs T", "T vs C"),
       fill = c("white", "darkgrey", "grey40"),
       title = "Max LFC Comparison",
       bty = 'n')
dev.off()

par(mfrow=c(1,1))

## PLOT: MA plots  -----

pdf("./output/loopProcessing/loopMA-byChrom.pdf",
    width = 12, height = 12)
par(mfrow=c(5,5))
par(mar=c(2,2,3,3))
for(c in unique(loops$seqnames1)){
  dat = loops[loops$seqnames1 == c,]
  plot(x = dat$avgCount, y = dat$ATlfc, pch = 19, col=alpha('black', .1),
       main = c, ylim=c(-4,4))
  abline(h=0)
  points(x = dat$avgCount[dat$ATdiff == T],
         y = dat$ATlfc[dat$ATdiff == T],
         pch = 19, col = alpha("red", 0.1))
}
plot(x = loops$avgCount, y = loops$ATlfc, pch = 19, col=alpha('black', .1),
     main = "A vs T", ylim=c(-4,4))
abline(h=0)
points(x = loops$avgCount[loops$ATdiff == T],
       y = loops$ATlfc[loops$ATdiff == T],
       pch = 19, col = alpha("red", 0.1))

par(mfrow=c(5,5))
for(c in unique(loops$seqnames1)){
  dat = loops[loops$seqnames1 == c,]
  plot(x = dat$avgCount, y = dat$AClfc, pch = 19, col=alpha('black', .1),
       main = c, ylim=c(-4,4))
  abline(h=0)
  points(x = dat$avgCount[dat$ACdiff == T],
         y = dat$AClfc[dat$ACdiff == T],
         pch = 19, col = alpha("red", 0.1))
}
plot(x = loops$avgCount, y = loops$AClfc, pch = 19, col=alpha('black', .1),
     main = "A vs C", ylim=c(-4,4))
abline(h=0)
points(x = loops$avgCount[loops$ACdiff == T],
       y = loops$AClfc[loops$ACdiff == T],
       pch = 19, col = alpha("red", 0.1))
dev.off()



## OUTPUT  -----
# save(loops, countMatrix, dds, dds.trans, diffLoops, allSig,
#      file = "./output/loopProcessing/MCF10_processedLoops-comprehensive.rda")

save(loops, countMatrix, dds, dds.trans, diffLoops, allSig,
     file = "./output/loopProcessing/MCF10_processedLoops-conservative.rda")
