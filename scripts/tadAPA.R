
# INITIALIZE --------------------------
library(mariner)
library(InteractionSet)
library(plotgardener)


# READ IN --------------------------

# tads
tads = readRDS("./output/tadProcessing/tads.RDS")
bounds = readRDS("./output/tadProcessing/tadBounds.RDS")

# hic files
hicFiles <- list.files(path = "/Volumes/MCF10/CantataData/proc/celltype/hicFiles/",
                       full.names = T)
names(hicFiles) = c("A",
                    "C",
                    "T")
hicFiles = hicFiles[c(1, 3, 2)]



# RUN --------------------------

# TAD APAs: Boundaries -----------
bound2mat <- function(x, window) {
  d <- x + window
  out <- GInteractions(d, d)
  mcols(out) = mcols(x)
  return(out)
}

## Extract 11x11 count matrices from all 3 hic files -----
diffBounds = bounds[bounds$sig == T,]

compCols = c("ATsig", "TCsig", "ACsig")
valCols = c("ATdiff", "TCdiff", "ACdiff")
fileComps = list(hicFiles[c(1,2)], 
                 hicFiles[c(2,3)],
                 hicFiles[c(1,3)])
n = 100

matrices <- Map(function(diffCol, valCol, files){
  
  sortedBounds = diffBounds[diffBounds[, diffCol] == T,]
  sortedBounds = sortedBounds[order(abs(sortedBounds[, valCol]), decreasing = T),]
  sortedBounds = GRanges(sortedBounds)
  
  set = sortedBounds
  set$change = sign(mcols(set)[, valCol])
  set = c(set[mcols(set)[, "change"] > 0][1:n],
          set[mcols(set)[, "change"] < 0][1:n])
  set = bound2mat(set, 250000)
  
  mat = pullHicMatrices(set,
                        binSize=10000,
                        files=files,
                        norm = "SCALE",
                        half='both')
}, diffCol = compCols, valCol = valCols, files = fileComps)

## Subset into gained/lost loops -------
res <- lapply(matrices, rowData)
gained <- Map(function(m, r, c){
  m[which(r[, c] > 0)]
}, m = matrices, r = res, c = valCols)
lost <- Map(function(m, r, c){
  m[which(r[, c] < 0)]
}, m = matrices, r = res, c = valCols)


## Aggregate gained/lost for control & treated -----
gainedA <- lapply(gained, function(m){
  apply(counts(m)[,,,1], c(1,2), median, na.rm=TRUE)
}) |> Reduce(f = '+')
gainedB <- lapply(gained, function(m){
  apply(counts(m)[,,,2], c(1,2), median, na.rm=TRUE)
}) |> Reduce(f = '+')

lostA <- lapply(lost, function(m){
  apply(counts(m)[,,,1], c(1,2), median, na.rm=TRUE)
}) |> Reduce(f = '+')
lostB <- lapply(lost, function(m){
  apply(counts(m)[,,,2], c(1,2), median, na.rm=TRUE)
}) |> Reduce(f = '+')


## Find common scale for gained/lost -------
gainedScale <- c(0, max(gainedA/3, gainedB/3))
lostScale <- c(0, max(lostA/3, lostB/3))

gainedScale = c(0, 40)
lostScale = gainedScale


## Visualize ---------

plotH = 2
plotB = 0.1
pageB = 0.5
pageW = plotH*2 + pageB*2
pageH = plotH + pageB*2

png("./output/tadAPA/Fig1G-gainedAPA.png",
    width = pageW, height = pageH, units = 'in', res = 300)

pageCreate(height = pageH, width = pageW,
           showGuides = F)

p = pgParams(width = plotH,
             height = plotH,
             y = pageB)

p1 = plotMatrix(gainedA/3, 
                params = p,
                x = pageB,
                zrange=gainedScale)
p2 = plotMatrix(gainedB/3, 
                params = p,
                x = pageB + plotH + plotB,
                zrange=gainedScale)
plotText("Weakened boundaries", 
         fontface = "bold",
         fontsize = 14,
         fontcolor = "#58869D",
         x = pageW/2,
         y = pageB - 0.3)
plotText(paste0("n = ", n*3), 
         fontsize = 12,
         x = pageW/2,
         y = pageB - 0.125)

dev.off()


png("./output/tadAPA/Fig1G-lostAPA.png",
    width = pageW, height = pageH, units = 'in', res = 300)

pageCreate(height = pageH, width = pageW,
           showGuides = F)

p3 = plotMatrix(lostA/3, 
                params = p,
                x = pageB,
                zrange=lostScale)
p4 = plotMatrix(lostB/3, 
                params = p,
                x = pageB + plotH + plotB,
                zrange=lostScale)
plotText("Strengthened boundaries", 
         fontface = "bold",
         fontsize = 14,
         fontcolor = "#86BEA2",
         x = pageW/2,
         y = pageB - 0.3)
plotText(paste0("n = ", n*3), 
         fontsize = 12,
         x = pageW/2,
         y = pageB - 0.125)

dev.off()


