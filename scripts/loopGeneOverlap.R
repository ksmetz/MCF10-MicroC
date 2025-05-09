
# INITIALIZE -----------------------------------------
library(mariner)
library(scales)
library(GenomicRanges)


# READ IN --------------------------------------------

## Loops --------
loops = as_ginteractions(readRDS("./output/loopProcessing/loops.rds"))

## PE RNA-seq gene info ---------
genes = GRanges(readRDS("./output/rnaProcessing/genes.rds"))


# OVERLAP LOOPS + GENES --------------------------------------------

## Find expressed gene promoters ------
smallestGroupSize <- 3
cutoff = 10
expressed <- rowSums(as.matrix(mcols(genes)[,21:29]) >= cutoff) >= smallestGroupSize

## Overlap loops + gene promoters ---------
proms = promoters(genes, upstream = 2000, downstream = 500)
gap = 10000

ov = findOverlaps(query = loops,
                  subject = proms[expressed,], 
                  maxgap = gap)

## Identify loop classes ---------
geneLoops = loops[queryHits(ov),]
mcols(geneLoops) = cbind(mcols(geneLoops), mcols(proms[expressed,][subjectHits(ov),]))

geneDiffLoops = geneLoops[geneLoops$cluster != "static",]

diffGeneLoops = geneLoops[geneLoops$AC.sig == T | geneLoops$AT.sig == T | geneLoops$TC.sig == T,]

diffGeneDiffLoops = diffGeneLoops[diffGeneLoops$cluster != "static",]

nonGeneLoops = loops[-unique(queryHits(ov)),]

nonLoopedGenes = genes[expressed,][-unique(subjectHits(ov)),]


# DEG ENRICHMENT ---------
## % differential genes in loops ---------
## 49% of all genes in loops are DEG vs 65% of genes in diff loops
genesInLoops = (unique(geneLoops$GENEID))
diffGenesInLoops = (unique(diffGeneLoops$GENEID))
length(diffGenesInLoops)/length(genesInLoops)

genesInDiffLoops = (unique(geneDiffLoops$GENEID))
diffGenesInDiffLoops = (unique(diffGeneDiffLoops$GENEID))
length(diffGenesInDiffLoops)/length(genesInDiffLoops)

## % differential genes in random sample ---------
## This is more than you'd expect by chance!
differential = (genes$AT.sig | genes$AC.sig | genes$TC.sig)
diffGenes = genes$GENEID[differential]
i = 1000
ratios = lapply(1:i, function(i){
  rand = genesInLoops[sample(1:length(genesInLoops), 
                             length(genesInDiffLoops))]
  randDiff = rand[(rand %in% diffGenes)]
  ratio = length(randDiff)/length(rand)
  return(ratio)
}) |> unlist()


## PLOT: DEG Enrichment ---------
pdf("./output/loopGeneOverlap/FigS5A-DEGdiffLoopEnrichment.pdf",
    width = 6, height = 6)

plot(density(ratios*100), xlim = c(30, 70),
     main = "Percent of genes differentially expressed")
abline(v = length(diffGenesInDiffLoops)/length(genesInDiffLoops)*100, 
       col = "red")
abline(v = length(diffGenesInLoops)/length(genesInLoops)*100, 
       col = "black", lty = 2)
abline(v = sum(differential)/sum(expressed)*100, 
       col = "grey", lty = 2)
legend("topleft",
       legend = c(paste0("Random sample (n = ", i, ")"),
                  paste0("Genes in static loops (",
                         signif(length(diffGenesInLoops)/length(genesInLoops)*100, digits = 3),
                         "%)"),
                  paste0("All genes (", 
                         signif(sum(differential)/sum(expressed)*100, digits = 3),
                         "%)"),
                  paste0("Genes in diff loops (", 
                         signif(length(diffGenesInDiffLoops)/length(genesInDiffLoops)*100, digits = 3),
                         "%)")),
       col = c("black", "black", "grey", "red"),
       lty = c(1, 2, 2, 1),
       bty = 'n')

dev.off()


# OUTPUT ------------------
saveRDS(geneLoops, "./output/loopGeneOverlap/geneLoops.rds")
saveRDS(diffGeneLoops, "./output/loopGeneOverlap/DEGloops.rds")
saveRDS(diffGeneDiffLoops, "./output/loopGeneOverlap/DEGdiffLoops.rds")