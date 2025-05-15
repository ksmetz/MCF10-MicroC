
# INITIALIZE --------------------------
library(mariner)
library(dplyr)
library(GenomicRanges)
library(InteractionSet)


# READ IN --------------------------

## TADs from SpectralTAD ------
tads = readRDS("./output/tadCalling/tads.RDS")

## IS from FAN-C ------
# Insulation scores (tech rep level)
techIS = lapply(
  list.files("./input/tads/techrepIS/", full.names = T),
  read.table,
  col.names = c("chr",
                "start",
                "end",
                "name",
                "score",
                "strand"))
names(techIS) = gsub("_contact_map_10kb-IS_100kb.bed", "",
                 list.files("./input/tads/techrepIS/"))
techIS = techIS[c(1:8, 17:24, 9:16)]

# Insulation scores (cell type level)
cellIS = lapply(
  list.files("./input/tads/celltypeIS/", full.names = T),
  read.table,
  col.names = c("chr",
                "start",
                "end",
                "name",
                "score",
                "strand"))
names(cellIS) = gsub("_merged_map_10kb-IS_100kb.bed", "",
                     list.files("./input/tads/celltypeIS/"))
cellIS = cellIS[c(1, 3, 2)]


## IS Boundaries from FAN-C --------
# Insulation score boundaries (cell type)
bounds = lapply(
  list.files("./input/tads/boundaries/", full.names = T),
  read.table,
  col.names = c("chr",
                "start",
                "end",
                "name",
                "score",
                "strand"))
names(bounds) = gsub("_merged_map_10kb-ISboundaries.bed", "",
                     list.files("./input/tads/boundaries/"))
bounds = bounds[c(1, 3, 2)]

## Denylist -------
deny = readRDS("./output/denyList/denylist.rds")



# PROCESS TADs + BOUNDS --------------------------

## Filter TADs ---------
#Combine TAD levels + remove duplicate TADs
tads = lapply(tads, bind_rows)
tads = lapply(tads, function(t){
  as.data.frame(t[!duplicated(t[,1:3]),])
})

# Remove small TADs
tads = lapply(tads, function(t){
  t = t[(t$end - t$start) >= 200000,]
})


# Convert to data.frame
tads = lapply(tads, function(t){
  data.frame(chr1 = t$chr,
             start1 = t$start,
             end1 = t$start + 10000,
             chr2 = t$chr,
             start2 = t$end - 10000,
             end2 = t$end,
             score = t$Sil_Score,
             level = t$Level)
})

# Bin to 10kb
tads = lapply(tads, assignToBins, binSize = 10000)

## Merge TADs ---------
# Merge TADs (21,682)
mergedTads <- 
  mergePairs(x = tads, 
             radius = 50000,
             method = "manhattan",
             pos = "center")

# Order + rebin accordingly
mergedTads = swapAnchors(mergedTads)
mergedTads = sort(mergedTads)
mergedTads = assignToBins(mergedTads, binSize = 10000, 
                          pos1="center", pos2="center")

# Add TAD call table 
dat = mariner:::.makeSrcMatrix(mergedTads)
mcols(mergedTads) = dat

# Add in metadata columns 
mergedTads = aggMetadata(mergedTads, 
                         columns = "score", 
                         funs = \(x){mean(x, na.rm=T)})
mergedTads = aggMetadata(mergedTads, 
                         columns = "level", 
                         funs = c("min", "max", "list"))
# 21,682 TADs
# 27,407 regions


## Process boundaries ------------

# Remove boundaries with low scores
bounds = lapply(bounds, function(df){
  df[which(df$score >= 0.5),]
})

# Remove duplicate boundaries
bounds = lapply(bounds, function(df){
  df[!duplicated(df),]
})


## Merge boundaries ---------
# Combine boundaries across cell types 
mergedBounds = do.call(rbind, bounds)
mergedBounds = mergedBounds[!duplicated(mergedBounds[,1:3]),]
mergedBounds = GRanges(mergedBounds)
# 20,294 bounds


## Filter boundaries ---------
# Remove boundaries with NA or Inf values
ISmat = cbind(techIS$MCF10_A_1_1[,1:3],
              lapply(techIS, function(df){df$score}))
keep = which(complete.cases(ISmat) & rowSums(!sapply(ISmat[,-(1:3)], is.finite)) == 0)
validBins = GRanges(ISmat[keep,])

ov = findOverlaps(query = mergedBounds, # 20,294
                  subject = validBins)  # 276,996 

mergedBounds = mergedBounds[unique(queryHits(ov)),] # 19,814


# Remove FAN-C boundaries overlapping denylist 
ov = findOverlaps(query = mergedBounds,
                  subject = deny) 
mergedBounds = mergedBounds[-(unique(queryHits(ov))),] # 19,433


## TAD-boundary overlaps --------
# Find SpectralTAD TADs overlapping with FAN-C IS boundaries 
boundaryOverlaps = findOverlaps(query = mergedTads,
                                subject = mergedBounds,
                                maxgap = 50000)
boundaryOverlaps1 = findOverlaps(query = mergedTads,
                                 subject = mergedBounds,
                                 maxgap = 50000,
                                 use.region = "first")
boundaryOverlaps2 = findOverlaps(query = mergedTads,
                                 subject = mergedBounds,
                                 maxgap = 50000,
                                 use.region = "second")

hits = unique(queryHits(boundaryOverlaps))
length(hits)/length(mergedTads) # 89.1%

hits = unique(subjectHits(boundaryOverlaps))
length(hits)/length(mergedBounds) # 76.5%

hits1 = unique(queryHits(boundaryOverlaps1))
hits2 = unique(queryHits(boundaryOverlaps2))
bothHits = hits1[hits1 %in% hits2]
length(bothHits)/length(mergedTads) # 61.1%

# Subset SpectralTAD TADs for those with IS bounds on both ends
mergedTads = mergedTads[bothHits,] # 21,682 --> 13,231
mergedTads = sort(mergedTads)
mergedTads = reduceRegions(mergedTads) # 13,231
# 13,253 TADs
# 17,150 regions

# Find boundaries of TADs (a subset of all IS boundaries)
tadBounds = regions(mergedTads) # 17,118
start(tadBounds) = start(tadBounds) + 1

# Comparing # of TAD bounds with full set of boundaries
length(tadBounds) # 17,118
length(mergedBounds) # 19,433
length(tadBounds)/length(mergedBounds) # 88.1%


## PLOT: TAD number by map -----
n = rowSums(as.data.frame(mcols(mergedTads)[,c("A", "T", "C")]))
dat = c(rev(table(n[mergedTads$A])), 
        rev(table(n[mergedTads$T])), 
        rev(table(n[mergedTads$C])))
dat = matrix(dat, ncol = 3, byrow = F)

library(RColorBrewer)

pdf(file="./output/tadProcessing/Fig1D-TADnumbers.pdf", 
    width=8, height=8)
barplot(dat,
        border = NA,
        names.arg = c("MCF10A", "MCF10AT1", "MCF10CA1a"),
        las = 1,
        yaxt = 'n',
        col = rev(brewer.pal(n = 5, "YlGnBu"))[2:4],
        main = "TADs called per map")
axis(side = 2, at = seq(0, 7.5e3, 2.5e3), las = 2, 
     xpd = T, lwd = 0, line = -0.8)
abline(h = seq(0, 12.5e3, 2.5e3), col = "grey")
dev.off()


# DIFFERENTIAL BOUNDS --------------------------

# Pulling IS for boundary bins 
ov = findOverlaps(query = resize(tadBounds, width = 50000, fix = "center"),
                  subject = validBins)

tadBoundIS = ISmat[keep, ][subjectHits(ov),] # 85,352 rows (boundaries)
tadBoundIS = split(tadBoundIS, queryHits(ov)) # 17,097 data.frames
tadBoundIS = lapply(tadBoundIS, function(df){ # 17,097 lists of 24 IS values
  dat = colMeans(df[, -(1:3)])
  return(dat)
}) 
tadBoundIS = do.call(rbind, tadBoundIS)

# Add boundary coordinates
tadBoundIS = cbind(as.data.frame(tadBounds)[unique(queryHits(ov)), c(1:3)],
                   tadBoundIS)

# Function for finding differential regions (by t-test)
findDiffBins <- function(mat){
  test = apply(mat, 1, function(r){
    val = t.test(x = r[1:8], y = r[9:16], alternative = "two.sided")
    return(val$p.value)
  })
  adjTest = p.adjust(test, method = "fdr")
  diff = apply(mat, 1, function(r){
    mean(r[9:16] - r[1:8])
  })
  return(list("pvalue" = test,
              "padj" = adjTest,
              "diff" = diff))
}

# Find differential bins among full IS matrix
ac = findDiffBins(mat = tadBoundIS[, c(4:11, 20:27)])
at = findDiffBins(mat = tadBoundIS[, c(4:11, 12:19)])
tc = findDiffBins(mat = tadBoundIS[, c(12:19, 20:27)])

# P-value histograms
hist(ac$pvalue)
hist(at$pvalue)
hist(tc$pvalue)

hist(ac$padj)
hist(at$padj)
hist(tc$padj)

# More changes between A/C > T/C > A/T
p = 0.01
sum(at$padj < p) # 567
sum(tc$padj < p) # 1693
sum(ac$padj < p) # 2314

# Differences in IS score
plot(density(ac$diff))
lines(density(ac$diff[ac$padj < p]), col = "red")
plot(density(at$diff))
lines(density(at$diff[at$padj < p]), col = "red")
plot(density(tc$diff))
lines(density(tc$diff[tc$padj < p]), col = "red")

range(abs(ac$diff[ac$padj < p]))
range(abs(at$diff[at$padj < p]))
range(abs(tc$diff[tc$padj < p]))

# Add differential info to TAD boundary IS matrix
tadBoundDiffInfo = cbind(tadBoundIS,
                         ATpadj = at$padj,
                         ATdiff = at$diff,
                         ATsig = at$padj < p,
                         TCpadj = tc$padj,
                         TCdiff = tc$diff,
                         TCsig = tc$padj < p,
                         ACpadj = ac$padj,
                         ACdiff = ac$diff,
                         ACsig = ac$padj < p)
tadBoundDiffInfo$sig = apply(tadBoundDiffInfo[, c("ATsig", "TCsig", "ACsig")], 1, any)
sum(tadBoundDiffInfo$sig)/nrow(tadBoundDiffInfo) # 19.8%

## 70-89% of boundaries get weaker with cancer progression
table(sign(tadBoundDiffInfo$ATdiff[tadBoundDiffInfo$ATsig == T]))/sum(tadBoundDiffInfo$ATsig == T)
table(sign(tadBoundDiffInfo$TCdiff[tadBoundDiffInfo$TCsig == T]))/sum(tadBoundDiffInfo$TCsig == T)
table(sign(tadBoundDiffInfo$ACdiff[tadBoundDiffInfo$ACsig == T]))/sum(tadBoundDiffInfo$ACsig == T)

## PLOT: # of diff boundaries ------
diffISbyComp = matrix(c(sum(tadBoundDiffInfo$ATdiff[tadBoundDiffInfo$ATsig == T] < 0),
                        sum(tadBoundDiffInfo$ATdiff[tadBoundDiffInfo$ATsig == T] > 0),
                        sum(tadBoundDiffInfo$TCdiff[tadBoundDiffInfo$TCsig == T] < 0),
                        sum(tadBoundDiffInfo$TCdiff[tadBoundDiffInfo$TCsig == T] > 0),
                        sum(tadBoundDiffInfo$ACdiff[tadBoundDiffInfo$ACsig == T] < 0),
                        sum(tadBoundDiffInfo$ACdiff[tadBoundDiffInfo$ACsig == T] > 0)),
                      ncol = 3, byrow = F)
colnames(diffISbyComp) = c("A vs T", "T vs C", "A vs C")

pdf("./output/tadProcessing/FigS2B-diffISbyComp.pdf",
    width = 6, height = 7)
barplot(diffISbyComp,
        col = c("#6099B0", "#A0D1BC"),
        border = NA,
        main = "Number of differential boundaries")
legend("topleft",
       legend = c("Weakened", "Strengthened"),
       text.col = c("#A0D1BC", "#6099B0"),
       bty = 'n')
dev.off()


## PLOT: Direction bias ------
sampler <- function(i, diffCol, valCol){
  # Sample all boundaries
  selection = sample(1:nrow(tadBoundDiffInfo), 
                     sum(tadBoundDiffInfo[, diffCol]))
  ISdat = tadBoundDiffInfo[selection,]
  
  # Find % weakened
  percWeakened = sum(ISdat[, valCol] > 0)/nrow(ISdat)
  return(percWeakened)
}

# Random sample (expected)
n = 1000
randAT = unlist(lapply(1:n, sampler, diffCol = "ATsig", valCol = "ATdiff"))
randTC = unlist(lapply(1:n, sampler, diffCol = "TCsig", valCol = "TCdiff"))
randAC = unlist(lapply(1:n, sampler, diffCol = "ACsig", valCol = "ACdiff"))

# Observed
obs = diffISbyComp[2,]/colSums(diffISbyComp)

# P-values
p = c(sum(randAT <= obs[1])/length(randAT),
      sum(randTC <= obs[2])/length(randTC),
      sum(randAC <= obs[3])/length(randAC))

# Plot
plotBias <- function(rand, title, ob, pval, diffCol){
  plot(density(rand), main = title, 
       xlim = c(0.5, 0.9), ylim = c(0, 50))
  abline(v = ob, col = "red")
  legend("topleft",
         legend = c(paste0("Random sample (", n, " perm.)"),
                    paste0("Observed (n = ", sum(tadBoundDiffInfo[, diffCol]), ")"),
                    paste0("P-value: ", 
                           signif(pval, 3))),
         col = c("black", "red", "black"),
         text.col = c("black", "red", "black"),
         bty = 'n',
         lty = c(1, 1, NA))
}


pdf("./output/tadProcessing/FigS2D-diffISdirBias.pdf",
    width = 6, height = 6)

par(mar=c(3,3,4,1))

plotBias(rand = randAT, title = "Percent of weakened boundaries (A vs T)", 
         ob = obs[1], pval = p[1], diffCol = "ATsig")
plotBias(rand = randTC, title = "Percent of weakened boundaries (T vs C)", 
         ob = obs[2], pval = p[2], diffCol = "TCsig")
plotBias(rand = randAC, title = "Percent of weakened boundaries (A vs C)", 
         ob = obs[3], pval = p[3], diffCol = "ACsig")

dev.off()


# CLUSTERING IS --------------------------

# Use cell-type level IS scores for clustering
cellISmat = cbind(cellIS$MCF10_A[,1:3],
                  lapply(cellIS, function(df){df$score}))
keep = which(complete.cases(ISmat) & rowSums(!sapply(ISmat[,-(1:3)], is.finite)) == 0)
validBins = GRanges(ISmat[keep,])

# Pulling IS for boundary bins
ov = findOverlaps(query = resize(tadBounds, width = 50000, fix = "center"),
                  subject = validBins)

cellBoundIS = cellISmat[keep, ][subjectHits(ov),] # 85,352 rows (boundaries)
cellBoundIS = split(cellBoundIS, queryHits(ov)) # 17,097 data.frames
cellBoundIS = lapply(cellBoundIS, function(df){ # 17,129 lists of 24 IS values
  dat = colMeans(df[, -(1:3)])
  return(dat)
}) 
cellBoundIS = do.call(rbind, cellBoundIS)

# Center and scale data
norm <- t(scale(t(cellBoundIS)))
normSig <- norm[tadBoundDiffInfo$sig,]

# Perform clustering
k = 4
set.seed(32)
clust = kmeans(normSig, centers=k)
cut = clust$cluster

cluster.order = c(3, 4, 2, 1)
names(cluster.order) = c("up.early", "up.late",
                         "down.early", "down.late")

## PLOT: Line plots -----
# Set graphical params
library(RColorBrewer)
library(scales)

pdf("./output/tadProcessing/Fig1G-clusterLines.pdf",
    width = 3, height = 12)

par(mar=c(3,3,3,3))
par(mfrow=c(4, 1))
k.pal = brewer.pal(n=4, "Spectral")

# Plot each cluster
n=0
for (i in cluster.order) {
  n=n+1
  
  # grab data
  dat = matrix(normSig[cut == i,], ncol=3)
  
  # make empty plot
  plot(x=1:3, y=colMeans(dat), type="n", ylim=c(-1.2,1.2), xaxt="n", 
       main=paste(names(cluster.order)[n],"\nn=",table(cut)[i],sep=""),
       xlab="Cancer Progression", ylab="Relative IS", 
       las = 2, bty = 'n', )
  abline(h=0, col = "grey")
  
  # add fill (std dev)
  polygon(c(1:3, 3:1), c(colMeans(dat)-colSds(dat), rev(colMeans(dat)+colSds(dat))), col=alpha(k.pal[n], .25), border=NA)
  
  # add median line
  lines(x=c(1:3), y=colMeans(dat),col=k.pal[n],lwd=2)
  
  # add axis
  axis(1, at=c(1:3), labels=c("A", "T", "C"))
}

dev.off()

## PLOT: Heatmap --------
library(pheatmap)
png(filename="./output/tadProcessing/Fig1G-clusterHeatmap.png",
    width=4, height=8, units="in", res=300)

par(mfrow=c(1,1))
par(mar=c(5,5,4,5))

pheatmap(normSig[order(match(cut, cluster.order)),], 
         cluster_rows = F, show_rownames = F,
         cluster_cols = F,
         annotation_row = data.frame(cluster = match(cut, cluster.order)[order(match(cut, cluster.order))],
                                     row.names = rownames(normSig[order(match(cut, cluster.order)),])),
         annotation_colors = list(cluster = k.pal),
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         annotation_legend = F,
         angle_col = "0")

dev.off()

## Diff TADs --------------
## Find TADs with 1 or 2 diff boundaries
diffTadBounds = GRanges(tadBoundDiffInfo[tadBoundDiffInfo$sig,])
boundaryOverlaps = findOverlaps(query = mergedTads,
                                subject = diffTadBounds)

## 38.4% (5,084) of TADs have at least 1 diff boundary 
hits = unique(queryHits(boundaryOverlaps)) 
length(hits)/length(mergedTads) 
length(hits)


## TABLE: Boundary info -----------
tadBounds = cbind(tadBoundDiffInfo,
                  data.frame(
                    cluster = cut[match(x=rownames(norm),
                                        table=names(cut))],
                    A_IS = cellBoundIS[,1],
                    T_IS = cellBoundIS[,2],
                    C_IS = cellBoundIS[,3],
                    A_ZSCR = norm[,1],
                    T_ZSCR = norm[,2],
                    C_ZSCR = norm[,3],
                    max_ZSCR = rowMaxs(norm)))

# Change name of cluster to descriptive name
tadBounds$cluster[which(is.na(tadBounds$cluster) == TRUE)] = "static"
for (i in 1:k){
  tadBounds$cluster[which(tadBounds$cluster == cluster.order[i])] = names(cluster.order[i])
}

# Save to table
write.table(x = tadBounds,
            file = "./output/tadProcessing/TableS3-TADbounds.txt",
            sep = "\t", quote = F, row.names = F, col.names = T)

## 28.8% of bounds lose IS (strengthened boundary)
sum(table(tadBounds$cluster)[1:2])/sum(tadBounds$cluster != "static")

## 71.2 of bounds lose IS (weakened boundary)
sum(table(tadBounds$cluster)[4:5])/sum(tadBounds$cluster != "static")

## TABLE: TAD info -----------
write.table(x = as.data.frame(mergedTads)[,1:13],
            file = "./output/tadProcessing/TableS4-TADs.txt",
            sep = "\t", quote = F, row.names = F)


# OUTPUT --------------------------
saveRDS(mergedTads, "./output/tadProcessing/tads.RDS")
saveRDS(tadBounds, "./output/tadProcessing/tadBounds.RDS")
saveRDS(cellIS, "./output/tadProcessing/IStracks.RDS")


