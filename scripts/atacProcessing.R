
# INTIALIZE -----------------------
library(GenomicRanges)
library(bedtoolsr)


# READ IN -----------------------

## Read in Stein ATAC peaks (narrowPeak) ----------
atacPeaks = list.files("./input/atac/peaks/", full.names = T)

## Read in ATAC alignments -------
atacBAMs = list.files("./input/atac/alignments/", full.names = T,
                      pattern = ".*bam$")

## Read in H3K27ac alignments --------
k27BAMs = list.files("./input/h3k27ac/alignments/", full.names = T,
                     pattern = ".*bam$")


# CONSENSUS PEAKS -----------------------
atacPeaks = lapply(atacPeaks, 
                   read.table, 
                   col.names = c("chr", "start", "stop", 
                                 "name", "score", "strand",
                                 "signalValue", "pValue", "qvalue", "peak"))
atacPeaks = lapply(atacPeaks, GRanges)

atacPeaks = do.call("c", atacPeaks) # 381893; 189-4811 bp
atacPeaks = reduce(atacPeaks) # 199299; 189-5323 bp

# EXTRACT ATAC COUNTS --------
atacCounts = bt.multicov(bams = paste0(atacBAMs, collapse = " "),
                         bed = as.data.frame(atacPeaks)) # 3:21-3:31
colnames(atacCounts) = c("chr",
                         "start",
                         "end",
                         "width",
                         "strand",
                         "MCF10A",
                         "MCF10AT1",
                         "MCF10CA1a")


# PREPARE ABC INPUT -----------------------
## Rank by ATAC counts per cell type per bp -------
ranks = lapply(atacCounts[, 6:8]/atacCounts[, 4], order, decreasing = T)

ranks = mapply("c", ranks[[1]], ranks[[2]], ranks[[3]], SIMPLIFY = F)

ranks = unlist(ranks)

ranks = unique(ranks)

topCounts = GRanges(atacCounts[ranks[1:150000], ])

## Extend regions to at least 500bp ---------
topCounts[width(topCounts) < 500,] = resize(topCounts[width(topCounts) < 500,],
                                            width = 500,
                                            fix = "center")

topCounts = as.data.frame(topCounts)
colnames(topCounts) = c("chr", "start", "end", "width", "strand",
                        "MCF10A", "MCF10AT1", "MCF10CA1a")

## TABLE: Input for ABC --------
write.table(topCounts,
            "./output/atacPeakProcessing/putativeEnhancers.txt", 
            quote = F, row.names = F, col.names = F, sep = "\t")




# EXTRACT H3K27ac COUNTS -----------------------
## Pull H3K27ac counts for ATAC-Seq regions --------
k27Counts = bt.multicov(bams = paste0(k27BAMs, collapse = " "),
                        bed = as.data.frame(atacPeaks)) # 3:21-3:31
colnames(k27Counts) = c("chr",
                        "start",
                        "end",
                        "width",
                        "strand",
                        "MCF10A_R1",
                        "MCF10A_R2",
                        "MCF10AT1_R1",
                        "MCF10AT1_R2",
                        "MCF10CA1a_R1",
                        "MCF10CA1a_R2")

# OUTPUTS ----------
saveRDS(atacPeaks, file = "./output/atacProcessing/atacPeaks.rds")

saveRDS(atacCounts,
        file = "./output/atacProcessing/ATACcounts.rds")

saveRDS(k27Counts,
        file = "./output/atacProcessing/K27acCounts.rds")