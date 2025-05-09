
# INITIALIZE --------------------------
library(mariner)
library(plotgardener)


# READ IN --------------------------

# loops
loops = as_ginteractions(readRDS("./output/loopProcessing/loops.rds"))

# hic files
hicFiles <- list.files(path = "/Volumes/MCF10/CantataData/proc/celltype/hicFiles/",
                       full.names = T)
names(hicFiles) = c("A",
                    "C",
                    "T")
hicFiles = hicFiles[c(1, 3, 2)]



# LOOP APAs    -----------------------------------------------------

## Extract 11x11 count matrices from all 3 hic files -----
sortedLoops = loops[order(loops$padj, decreasing = F),]
sortedLoops = sortedLoops[sortedLoops$span >= 130000,]

compCols = c("ATdiff", "TCdiff", "ACdiff")
lfcCols = c("ATlfc", "TClfc", "AClfc")
fileComps = list(hicFiles[c(1,2)], 
                 hicFiles[c(2,3)],
                 hicFiles[c(1,3)])
n = 100

matrices <- Map(function(diffCol, lfcCol, files){
  
  set = sortedLoops[mcols(sortedLoops)[,diffCol] == T,]
  set = c(set[mcols(set)[,lfcCol] > 0, ][1:n],
          set[mcols(set)[,lfcCol] < 0, ][1:n])
  
  mat = set |>
    pixelsToMatrices(buffer=5) |>
    removeShortPairs() |>
    pullHicMatrices(binSize=10000,
                    files=files,
                    norm = "SCALE",
                    half='upper')
}, diffCol = compCols, lfcCol = lfcCols, files = fileComps)

## Subset into gained/lost loops -------
res <- lapply(matrices, rowData)
gained <- Map(function(m, r, c){
  m[which(r[, c] > 0)]
}, m = matrices, r = res, c = lfcCols)
lost <- Map(function(m, r, c){
  m[which(r[, c] < 0)]
}, m = matrices, r = res, c = lfcCols)

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


## PLOT: APA plots ---------
plotH = 2
plotB = 0.1
pageB = 0.5
pageW = plotH*2 + pageB*2
pageH = plotH + pageB*2

png("./output/loopAPA/Fig1H-gainedAPA.png",
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
plotText("Strengthened loops", 
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


png("./output/loopAPA/Fig1H-lostAPA.png",
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
plotText("Weakened loops", 
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

