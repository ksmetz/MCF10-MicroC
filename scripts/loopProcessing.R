
# INITIALIZE  ---------------------------------------------

library(mariner)
library(DESeq2)
library(InteractionSet)
library(scales)
library(RColorBrewer)
library(colorspace)
library(pheatmap)



# READ IN     ---------------------------------------------

## Loop calls -------
celldir = "./input/loops/celltype/"
megadir = "./input/loops/mega/mega_loops.txt"

loopFiles = list.files(celldir,  
                       pattern = "*loops.txt", 
                       full.names = T)
loopFiles = c(loopFiles, megadir)
names(loopFiles) = c("MCF10_A", "MCF10_C", "MCF10_T", "mega")
loopFiles = loopFiles[c("MCF10_A", "MCF10_T", "MCF10_C", "mega")]

loopSets = lapply(loopFiles, read.table, header = T)
loopsGR = lapply(loopSets, as_ginteractions)


## Tech rep hic files -------
trHicPaths = list.files("/Volumes/MCF10/CantataData/proc/techrep/hicFiles/",
                        full.names = T)
trHicPaths = trHicPaths[c(1:8, 17:24, 9:16)]


## Denylist -------
deny = readRDS("./output/denyList/denylist.rds")



# MERGING LOOPS     ------------------------------------------------------

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

## Remove loops in denylist regions (29,627 --> 29,205)
ov = findOverlaps(query = mergedLoops,
                  subject = deny)
mergedLoops = mergedLoops[-(unique(queryHits(ov))),]


# PULL COUNTS     ------------------------------------------------------
# Read in counts from hic files; runs in ~6m
counts = pullHicPixels(x = mergedLoops, 
                       files = trHicPaths, 
                       binSize = 10000, 
                       norm = "NONE")


# DIFFERENTIAL ANALYSIS     ------------------------------------------------------
# Pull count matrix for DESeq
countMatrix=counts(counts)
countMatrix=suppressWarnings(apply(countMatrix, 2, as.numeric))
rownames(countMatrix) = paste0("loop_", 1:nrow(countMatrix))

# Make DNA copy matrix
loops = as.data.frame(mergedLoops)
dnacopy <- matrix(1, nrow=nrow(countMatrix), ncol=ncol(countMatrix))

# Note: Based on readouts from NeoLoopFinder CNV
dnacopy[loops$seqnames1 == "chr1" & loops$start1 >= 143200000,  1:8] <- 1.39
dnacopy[loops$seqnames1 == "chr1" & loops$start1 >= 143200000,  9:16] <- 0.96
dnacopy[loops$seqnames1 == "chr1" & loops$start1 >= 143200000,  17:24] <- 0.93

dnacopy[loops$seqnames1 == "chr1" & loops$start1 >= 201000000,  1:8] <- 1.62
dnacopy[loops$seqnames1 == "chr1" & loops$start1 >= 201000000,  9:16] <- 0.87
dnacopy[loops$seqnames1 == "chr1" & loops$start1 >= 201000000,  17:24] <- 0.86

dnacopy[loops$seqnames1 == "chr3" & loops$start1 >= 69000000,  1:8] <- 0.98
dnacopy[loops$seqnames1 == "chr3" & loops$start1 >= 69000000,  9:16] <- 0.93
dnacopy[loops$seqnames1 == "chr3" & loops$start1 >= 69000000,  17:24] <- 1.44

dnacopy[loops$seqnames1 == "chr8" & loops$start1 >= 99000000,  1:8] <- 1.56
dnacopy[loops$seqnames1 == "chr8" & loops$start1 >= 99000000,  9:16] <- 0.92
dnacopy[loops$seqnames1 == "chr8" & loops$start1 >= 99000000,  17:24] <- 0.94

dnacopy[loops$seqnames1 == "chr9" & loops$end1 <= 39000000,  1:8] <- 0.82
dnacopy[loops$seqnames1 == "chr9" & loops$end1 <= 39000000,  9:16] <- 1.11
dnacopy[loops$seqnames1 == "chr9" & loops$end1 <= 39000000,  17:24] <- 1.49

dnacopy[loops$seqnames1 == "chr9" & loops$start1 >= 65000000,  1:8] <- 0.95
dnacopy[loops$seqnames1 == "chr9" & loops$start1 >= 65000000,  9:16] <- 1.02
dnacopy[loops$seqnames1 == "chr9" & loops$start1 >= 65000000,  17:24] <- 1.35

dnacopy[loops$seqnames1 == "chr10", 1:8] <- 1.04
dnacopy[loops$seqnames1 == "chr10", 9:16] <- 1.04
dnacopy[loops$seqnames1 == "chr10", 17:24] <- 1.25

dnacopy[loops$seqnames1 == "chr13" & loops$start >= 51000000, 1:8] <- 1.34
dnacopy[loops$seqnames1 == "chr13" & loops$start >= 51000000, 9:16] <- 0.99
dnacopy[loops$seqnames1 == "chr13" & loops$start >= 51000000, 17:24] <- 1.10

dnacopy[loops$seqnames1 == "chr19" & loops$start >= 34000000, 1:8] <- 1.06
dnacopy[loops$seqnames1 == "chr19" & loops$start >= 34000000, 9:16] <- 1.37
dnacopy[loops$seqnames1 == "chr19" & loops$start >= 34000000, 17:24] <- 1.00

dnacopy[loops$seqnames1 == "chr20", 1:8] <- 0.98
dnacopy[loops$seqnames1 == "chr20", 9:16] <- 1.46
dnacopy[loops$seqnames1 == "chr20", 17:24] <- 0.92


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


## Add additional loop info to DF -----
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

## Summary DESeq info -----
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


## PLOT: Loop PCA -------
sample.pal = c("steelblue", "orange", "firebrick")

pdf(file="./output/loopProcessing/FigS1D-PCAplot.pdf", 
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


## PLOT: Diff loop barplot -------
diffLoopByComp = matrix(c(sum(loops$ATlfc[loops$ATdiff == T] > 0),
                          sum(loops$ATlfc[loops$ATdiff == T] < 0),
                          sum(loops$TClfc[loops$TCdiff == T] > 0),
                          sum(loops$TClfc[loops$TCdiff == T] < 0),
                          sum(loops$AClfc[loops$ACdiff == T] > 0),
                          sum(loops$AClfc[loops$ACdiff == T] < 0)),
                        ncol = 3, byrow = F)
colnames(diffLoopByComp) = c("A vs T", "T vs C", "A vs C")

pdf("./output/loopProcessing/FigS2A-diffLoopByComp.pdf",
    width = 6, height = 6)
barplot(diffLoopByComp,
        col = c("#6099B0", "#A0D1BC"),
        border = NA,
        main = "Number of differential loops")
legend("topleft",
       legend = c("Weakened", "Strengthened"),
       text.col = c("#A0D1BC", "#6099B0"),
       bty = 'n')
dev.off()


## PLOT: Loop number by map -------
n = rowSums(loops[,c("MCF10_A", "MCF10_T", "MCF10_C")])
dat = c(rev(table(n[loops$MCF10_A])), 
        rev(table(n[loops$MCF10_T])), 
        rev(table(n[loops$MCF10_C])))
dat = matrix(dat, ncol = 3, byrow = F)

pdf(file="./output/loopProcessing/Fig1E-loopNumbers.pdf", 
    width=8, height=8)
barplot(dat,
        border = NA,
        names.arg = c("MCF10A", "MCF10AT1", "MCF10CA1a"),
        las = 1,
        yaxt = 'n',
        col = rev(brewer.pal(n = 5, "YlGnBu"))[2:4],
        main = paste0(nrow(loops), " Unique\nChromatin Loops"))
axis(side = 2, at = seq(0, 30e3, 5e3), las = 2, 
     xpd = T, lwd = 0, line = -0.8)
abline(h = seq(0, 30e3, 5e3), col = "grey")
dev.off()



# CLUSTERING     ------------------------------------------------------

# Cluster settings
k = 4
k.pal = brewer.pal(n=4, "Spectral")

# Function for combining replicate columns
combineReps <- function(matrix, new.colnames)
{
  newMatrix = c()
  for (i in c(1, 9, 17))
  {
    tempMatrix = matrix[,i:(i+7)]
    newCol = apply(tempMatrix,1,mean)
    newMatrix  = cbind(newMatrix,newCol)
  }
  rownames(newMatrix) = rownames(matrix)
  
  if(!is.null(new.colnames)){
    colnames(newMatrix) = new.colnames
  }
  
  return(newMatrix)
}


## Normalize
# Build a matrix of transformed loop counts
countMat = assay(dds.trans)

# Center and scale data
countMat.norm <- (countMat - rowMeans(countMat))/rowSds(countMat + .5)

# Subset for significant loops
countMatSig.norm <- countMat.norm[which(rownames(countMat.norm) %in% allSig),]

# Combine replicates
countMat.combo <- combineReps(countMat, 
                              new.colnames=c("A", "T", "C"))
countMat.norm.combo <- combineReps(countMat.norm, 
                                   new.colnames=c("A", "T", "C"))
countMatSig.norm.combo <- combineReps(countMatSig.norm, 
                                      new.colnames=c("A", "T", "C"))

## Cluster
# Set seed to preserve manual ordering
seed = 32
set.seed(seed)

# Perform clustering
clust = kmeans(countMatSig.norm, centers=k)
cut = clust$cluster

# Order clusters
clusterOrder = c(3, 2, 4, 1) 
names(clusterOrder) = c("up.early", "up.late", 
                        "down.early", "down.late")



## PLOT: Line Plots -----------
# pdf("./figs/MCF10analysis/prelimLoops/MCF10_clusters/clusterLinePlots.pdf", 
#     width=9, height=6)
pdf("./output/loopProcessing/Fig1H-clusterLines.pdf", 
    width=3, height=12)

# Plot each cluster
par(mar=c(3,3,3,3))
par(mfrow=c(2,ceiling(k/2)))
par(mfrow=c(4, 1))

n=0
for (i in clusterOrder) {
  n=n+1
  
  # grab data
  dat = countMatSig.norm.combo[cut == i,]
  
  # make empty plot
  plot(x=1:3, 
       y=colMeans(countMatSig.norm.combo[cut==i,]), 
       type="n",
       main=paste(names(clusterOrder)[n],"\nn=",table(cut)[i],sep=""),
       ylim=c(-1,1), 
       xaxt="n", 
       xlab="Cancer progression", 
       ylab="Relative loop strength", 
       bty = 'n', las = 2)
  
  # add fill (std dev)
  polygon(c(1:3, 3:1), 
          c(colMeans(dat)-colSds(dat), 
            rev(colMeans(dat)+colSds(dat))), 
          col=lighten(k.pal[n], amount=.85), border=NA)
  
  # and transparent grey lines
  #apply(countMatSig.norm.combo[cut==i,],1,lines,col=adjustcolor("grey",alpha.f=0.4), x=c(0, .25, .5, .75, 1,2,3,4))
  abline(h=0, col=alpha("black", .2))
  
  # add median line
  lines(x=1:3,
        y=colMeans(countMatSig.norm.combo[cut==i,]),
        col=k.pal[n],
        lwd=2)
  
  # add axis
  axis(1, at=1:3, 
       labels=c("A", "T", "C"))
}

dev.off()



## PLOT: Heatmap -----------

# Plot heatmap of k clusters 
png(filename="./output/loopProcessing/Fig1H-clusterHeatmap.png",
    width=4, height=8, units="in", res=300)

# Plot heatmap
par(mfrow=c(1,1))
par(mar=c(5,5,5,5))

pheatmap(countMatSig.norm.combo[order(match(cut, clusterOrder)),], 
         cluster_rows = F, show_rownames = F,
         cluster_cols = F,
         annotation_row = data.frame(cluster = match(cut, clusterOrder)[order(match(cut, clusterOrder))],
                                     row.names = rownames(countMatSig.norm.combo[order(match(cut, clusterOrder)),])),
         annotation_colors = list(cluster = k.pal),
         annotation_legend = F,
         labels_col = c("MCF10A", "MCF10AT1", "MCF10CA1a"),
         angle_col = "0")


dev.off()


## TABLE: Loop info -----------
loops = cbind(loops,
              data.frame(
                cluster = cut[match(x=rownames(countMat.norm.combo),
                                    table=names(cut))],
                A_ZSCR = countMat.norm.combo[,1],
                T_ZSCR = countMat.norm.combo[,2],
                C_ZSCR = countMat.norm.combo[,3],
                max_ZSCR = rowMaxs(countMat.norm.combo),
                A_VST = countMat.combo[,1],
                T_VST = countMat.combo[,2],
                C_VST = countMat.combo[,3]))

# Change name of cluster to descriptive name
loops$cluster[which(is.na(loops$cluster) == TRUE)] = "static"
for (i in 1:k){
  loops$cluster[which(loops$cluster == clusterOrder[i])] = names(clusterOrder[i])
}

## 53.8% of diff loops weakened
sum(table(loops$cluster)[1:2])/sum(loops$cluster != "static")

## 46.2% of diff loops strengthened
sum(table(loops$cluster)[4:5])/sum(loops$cluster != "static")

# Save to table
write.table(x = loops,
            file = "./output/loopProcessing/TableS2-loops.txt",
            sep = "\t", quote = F, row.names = F, col.names = T)


# OUTPUT      -----------------------------------------------------------------------------
saveRDS(loops, 
        file="./output/loopProcessing/loops.rds")

save(countMatrix, 
     countMat.combo, countMat.norm.combo, countMatSig.norm.combo, 
     cut, clusterOrder,
     dds, dds.trans, diffLoops, allSig,
     file = "./output/loopProcessing/loopData.rda")
