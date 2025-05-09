
# INITIALIZE  ---------------------------------------------
library(DESeq2)
library(scales)
library(mariner)
library(GenomicRanges)



# READ IN     ---------------------------------------------

## CNVs from NeoLoopFinder
dir = "./input/CNVs/"
pro = lapply(list.files(path = dir, pattern = "*1000.*profile*", 
                        full.names = T), 
             read.table, col.names = c("seqnames", "start", "stop", "val"))


## Loops
loops = as_ginteractions(readRDS("./output/loopProcessing/loops.rds"))

## Loop count matrix
load("./output/loopProcessing/loopData.rda", verbose = T)
seqDepthFactors = estimateSizeFactorsForMatrix(countMatrix)

## EP count matrix
load("./output/abcProcessing/diffContactObjects.rda", verbose = T)
EPseqDepthFactors = estimateSizeFactorsForMatrix(countMatrix)

## EP Contacts
contacts = readRDS("./output/abcProcessing/abcPairs.rds")



# LOOP CNV PLOT    ---------------------------------------------

## Reorder CNV chromosomes
pro = lapply(pro, function(df){
  df$seqnames = factor(df$seqnames, levels = paste0('chr', c(1:22, "X")))
  df = df[order(df$seqnames),]
  return(df)
})


## Find CNV values for loops
ov = findOverlaps(subject = GRanges(pro[[1]]), 
                  query = loops,
                  select = "first")
loopBands = split(loops, ov)

# Function for highlighting region of CNV plot
highlightRegion <- function(df, chrom, start, stop, color, alpha, 
                            ymin = 0, ymax = 3, plot = T){
  region = which(df[,1] == chrom & df[,2] >= start & df[,3] <= stop)
  if(plot == T){
    rect(xleft = min(region), xright = max(region),
         ytop = ymax, ybottom = ymin, 
         col = alpha(color, alpha), border = NA)
  }
  return(df[region,])
}


## PLOT: CNVs for loops across genome --------
pdf("./output/loopCNV/FigS1C-CNVplot.pdf",
    width = 12, height = 3)

par(mfrow=c(1,1))
par(mar=c(4,4,1,1))
plot(x = 1:nrow(pro[[1]]),
     y = pro[[1]]$val,
     type = 'l', col = "steelblue",
     ylim = c(0, 3),
     xlab = "Genomic position",
     ylab = "Neoloop CNV factor")
lines(x = 1:nrow(pro[[3]]),
      y = pro[[3]]$val, col = "orange")
lines(x = 1:nrow(pro[[2]]),
      y = pro[[2]]$val, col = "firebrick")
abline(h = c(0, 1, 2), lty = 2)

cols = pro[[1]]$seqnames
levels(cols) = c(rep(c("black", "white"), times = 11), "black")
points(x = 1:nrow(pro[[1]]),
       y = rep(0, times = nrow(pro[[1]])),
       pch = 19, 
       col = cols)

tmp = highlightRegion(df = pro[[1]], 
                      chrom = "chr1", start = 143200000, stop = 201000000, 
                      color = "steelblue", alpha = 1.5/5)
mean(tmp$val) # 1.669
tmp = highlightRegion(df = pro[[1]], 
                      chrom = "chr1", start = 201000000, stop = 248956422, 
                      color = "steelblue", alpha = 2/5)
mean(tmp$val) # 2.152
tmp = highlightRegion(df = pro[[2]], 
                      chrom = "chr3", start = 69000000, stop = 198295559, 
                      color = "firebrick", alpha = 1.5/5)
mean(tmp$val) # 1.501
tmp = highlightRegion(df = pro[[1]], 
                      chrom = "chr8", start = 99000000, stop = 145138636, 
                      color = "steelblue", alpha = 1.5/5)
mean(tmp$val) # 1.710
tmp = highlightRegion(df = pro[[2]], 
                      chrom = "chr9", start = 1, stop = 39000000, 
                      color = "firebrick", alpha = 1.5/5)
mean(tmp$val) # 1.465
tmp = highlightRegion(df = pro[[1]], 
                      chrom = "chr9", start = 1, stop = 39000000, 
                      color = "steelblue", alpha = 0.8/5)
mean(tmp$val) # 0.796
tmp = highlightRegion(df = pro[[2]], 
                      chrom = "chr9", start = 65000000, stop = 138394717, 
                      color = "firebrick", alpha = 1.4/5)
mean(tmp$val) # 1.366
tmp = highlightRegion(df = pro[[3]], 
                      chrom = "chr20", start = 1, stop = 64444167, 
                      color = "orange", alpha = 1.2/5)
mean(tmp$val) # 1.223
tmp = highlightRegion(df = pro[[3]], 
                      chrom = "chr19", start = 34000000, stop = 58617616, 
                      color = "orange", alpha = 1/5)
mean(tmp$val) # 1.004
tmp = highlightRegion(df = pro[[2]],
                      chrom = "chr10", start = 1, stop = 133000000, 
                      color = "firebrick", alpha = 1.2/5)
mean(tmp$val) # 1.182
tmp = highlightRegion(df = pro[[1]],
                      chrom = "chr13", start = 51000000, stop = 114364328, 
                      color = "steelblue", alpha = 1.4/5)
mean(tmp$val) # 1.012


dev.off()


# LOOP DESeq FACTOR PLOT    ---------------------------------------------

# Find the total loop counts for each band
karyoFactors = lapply(loopBands, function(d){
  tmp = as.data.frame(mcols(d))[,17:40]
  tmp = t(apply(tmp, 1, function(r){r/seqDepthFactors}))
  return(estimateSizeFactorsForMatrix(tmp, type = "poscounts"))
})
karyoFactors = do.call("rbind", karyoFactors)

# Combine bands + loop count summary info
regions = GRanges(pro[[1]][as.numeric(names(loopBands)),])
mcols(regions) = cbind(mcols(regions), karyoFactors)

regionDF = as.data.frame(regions)
regionDF = regionDF[complete.cases(regionDF),]

regionDF$A = rowMeans(regionDF[,7:14])
regionDF$T = rowMeans(regionDF[,15:22])
regionDF$C = rowMeans(regionDF[,23:30])


# Plot DESeq normalization factors across genome
plot(x = 1:nrow(regionDF),
     y = regionDF$A,
     type = 'l', col = "steelblue",
     ylim = c(0, 3),
     xlab = "Genomic position",
     ylab = "Loop count DESeq factor")
lines(x = 1:nrow(regionDF),
      y = regionDF$T, col = "orange")
lines(x = 1:nrow(regionDF),
      y = regionDF$C, col = "firebrick")
abline(h = c(0, 1, 2), lty = 2)

cols = unique(regionDF$seqnames)
cols = factor(regionDF$seqnames, levels = cols)
levels(cols) = c(rep(c("black", "white"), times = 11), "black")
points(x = 1:nrow(regionDF),
       y = rep(0, times = nrow(regionDF)),
       pch = 19, 
       col = cols)

tmp = highlightRegion(df = regionDF, 
                      chrom = "chr1", start = 143200000, stop = 201000000, 
                      color = "steelblue", alpha = 1.4/5)
mean(tmp$A) # 1.388
mean(tmp$T) # 0.960
mean(tmp$C) # 0.932
tmp = highlightRegion(df = regionDF, 
                      chrom = "chr1", start = 201000000, stop = 248956422, 
                      color = "steelblue", alpha = 1.6/5)
mean(tmp$A) # 1.617
mean(tmp$T) # 0.874
mean(tmp$C) # 0.858
tmp = highlightRegion(df = regionDF, 
                      chrom = "chr3", start = 69000000, stop = 198295559, 
                      color = "firebrick", alpha = 1.44/5)
mean(tmp$A) # 0.978
mean(tmp$T) # 0.932
mean(tmp$C) # 1.441
tmp = highlightRegion(df = regionDF, 
                      chrom = "chr8", start = 99000000, stop = 145138636, 
                      color = "steelblue", alpha = 1.5/5)
mean(tmp$A) # 1.560
mean(tmp$T) # 0.919
mean(tmp$C) # 0.938
tmp = highlightRegion(df = regionDF,
                      chrom = "chr9", start = 1, stop = 39000000, 
                      color = "firebrick", alpha = 1.5/5)
tmp = highlightRegion(df = regionDF,
                      chrom = "chr9", start = 1, stop = 39000000, 
                      color = "steelblue", alpha = 0.8/5)
mean(tmp$A) # 0.818
mean(tmp$T) # 1.112
mean(tmp$C) # 1.487
tmp = highlightRegion(df = regionDF,
                      chrom = "chr9", start = 65000000, stop = 138394717, 
                      color = "firebrick", alpha = 1.4/5)
mean(tmp$A) # 0.944
mean(tmp$T) # 1.017
mean(tmp$C) # 1.344
tmp = highlightRegion(df = regionDF,
                      chrom = "chr10", start = 1, stop = 133000000, 
                      color = "firebrick", alpha = 1.2/5)
mean(tmp$A) # 1.043
mean(tmp$T) # 1.043
mean(tmp$C) # 1.246
tmp = highlightRegion(df = regionDF,
                      chrom = "chr13", start = 51000000, stop = 114364328, 
                      color = "steelblue", alpha = 1.4/5)
mean(tmp$A) # 1.339
mean(tmp$T) # 0.989
mean(tmp$C) # 1.105
tmp = highlightRegion(df = regionDF,
                      chrom = "chr19", start = 34000000, stop = 58617616, 
                      color = "orange", alpha = 1.4/5)
mean(tmp$A) # 1.056
mean(tmp$T) # 1.370
mean(tmp$C) # 0.996
tmp = highlightRegion(df = regionDF,
                      chrom = "chr20", start = 1, stop = 64444167, 
                      color = "orange", alpha = 1.4/5)
mean(tmp$A) # 0.984 
mean(tmp$T) # 1.456
mean(tmp$C) # 0.920



# CONTACT DESeq FACTOR PLOT    ---------------------------------------------

abc2gi <- function(abc){
  gi = as_ginteractions(
    data.frame(seqnames1 = abc$chr,
               start1 = abc$start,
               end1 = abc$end,
               seqnames2 = abc$chr,
               start2 = abc$TargetGeneTSS,
               end2 = abc$TargetGeneTSS,
               abc[, -c(1:3)]))
}
contacts = abc2gi(contacts)

## Find CNV values for contacts
ov = findOverlaps(subject = GRanges(pro[[1]]), 
                  query = contacts,
                  select = "first")
contactBands = split(contacts, ov)

# Find the total loop counts for each band
karyoFactors = lapply(contactBands, function(d){
  tmp = as.data.frame(mcols(d))[,164:187]
  tmp = t(apply(tmp, 1, function(r){r/EPseqDepthFactors}))
  return(estimateSizeFactorsForMatrix(tmp, type = "poscounts"))
})
karyoFactors = do.call("rbind", karyoFactors)

# Combine bands + loop count summary info
regions = GRanges(pro[[1]][as.numeric(names(contactBands)),])
mcols(regions) = cbind(mcols(regions), karyoFactors)

regionDF = as.data.frame(regions)
regionDF = regionDF[complete.cases(regionDF),]

regionDF$A = rowMeans(regionDF[,7:14])
regionDF$T = rowMeans(regionDF[,15:22])
regionDF$C = rowMeans(regionDF[,23:30])


# Plot DESeq normalization factors across genome
plot(x = 1:nrow(regionDF),
     y = regionDF$A,
     type = 'l', col = "steelblue",
     ylim = c(0, 3),
     xlab = "Genomic position",
     ylab = "Contact count DESeq factor")
lines(x = 1:nrow(regionDF),
      y = regionDF$T, col = "orange")
lines(x = 1:nrow(regionDF),
      y = regionDF$C, col = "firebrick")
abline(h = c(0, 1, 2), lty = 2)

cols = unique(regionDF$seqnames)
cols = factor(regionDF$seqnames, levels = cols)
levels(cols) = c(rep(c("black", "white"), times = 11), "black")
points(x = 1:nrow(regionDF),
       y = rep(0, times = nrow(regionDF)),
       pch = 19, 
       col = cols)

tmp = highlightRegion(df = regionDF, 
                      chrom = "chr1", start = 143200000, stop = 201000000, 
                      color = "steelblue", alpha = 1.4/5)
mean(tmp$A) # 1.325
mean(tmp$T) # 0.959
mean(tmp$C) # 0.950
tmp = highlightRegion(df = regionDF, 
                      chrom = "chr1", start = 201000000, stop = 248956422, 
                      color = "steelblue", alpha = 1.6/5)
mean(tmp$A) # 1.565
mean(tmp$T) # 0.867
mean(tmp$C) # 0.883
tmp = highlightRegion(df = regionDF, 
                      chrom = "chr3", start = 69000000, stop = 198295559, 
                      color = "firebrick", alpha = 1.44/5)
mean(tmp$A) # 0.945
mean(tmp$T) # 0.912
mean(tmp$C) # 1.428
tmp = highlightRegion(df = regionDF, 
                      chrom = "chr8", start = 99000000, stop = 145138636, 
                      color = "steelblue", alpha = 1.5/5)
mean(tmp$A) # 1.477
mean(tmp$T) # 0.891
mean(tmp$C) # 0.925
tmp = highlightRegion(df = regionDF,
                      chrom = "chr9", start = 1, stop = 39000000, 
                      color = "firebrick", alpha = 1.5/5)
tmp = highlightRegion(df = regionDF,
                      chrom = "chr9", start = 1, stop = 39000000, 
                      color = "steelblue", alpha = 0.8/5)
mean(tmp$A) # 0.813
mean(tmp$T) # 1.063
mean(tmp$C) # 1.434
tmp = highlightRegion(df = regionDF,
                      chrom = "chr9", start = 65000000, stop = 138394717, 
                      color = "firebrick", alpha = 1.4/5)
mean(tmp$A) # 0.929
mean(tmp$T) # 0.988
mean(tmp$C) # 1.320
tmp = highlightRegion(df = regionDF,
                      chrom = "chr10", start = 1, stop = 133000000, 
                      color = "firebrick", alpha = 1.2/5)
mean(tmp$A) # 0.969
mean(tmp$T) # 1.020
mean(tmp$C) # 1.229
tmp = highlightRegion(df = regionDF,
                      chrom = "chr13", start = 51000000, stop = 114364328, 
                      color = "steelblue", alpha = 1.4/5)
mean(tmp$A) # 1.189
mean(tmp$T) # 0.999
mean(tmp$C) # 1.074
tmp = highlightRegion(df = regionDF,
                      chrom = "chr19", start = 34000000, stop = 58617616, 
                      color = "orange", alpha = 1.4/5)
mean(tmp$A) # 1.018
mean(tmp$T) # 1.219
mean(tmp$C) # 0.970
tmp = highlightRegion(df = regionDF,
                      chrom = "chr20", start = 1, stop = 64444167, 
                      color = "orange", alpha = 1.4/5)
mean(tmp$A) # 0.960 
mean(tmp$T) # 1.331
mean(tmp$C) # 0.963
