
# INTIALIZE -----------------------
library(GenomicRanges)
library(mariner)


# READ IN -----------------------

## ATAC peaks ------
enh = GRanges(readRDS(file = "./output/H3K27acProcessing/putativeEnhancers.rds"))

## Micro-C loops --------
loops = as_ginteractions(readRDS(file = "./output/loopProcessing/loops.rds"))

# DIFFERENTIAL LOOP ANCHORS -----------------------

## Gained loops -----
gainedLoops = loops[(loops$ATdiff & loops$ATlfc > 0) | (loops$ACdiff & loops$AClfc > 0) | (loops$TCdiff & loops$TClfc > 0),] 
ov = findOverlaps(query = gainedLoops,
                  subject = enh)
write.table(x = enh[unique(subjectHits(ov)),], 
            quote = F, col.names = F, row.names = F, sep = "\t",
            file = "./output/motifPrep/gainedLoops.bed")

## Lost loops -------
lostLoops = loops[(loops$ATdiff & loops$ATlfc < 0) | (loops$ACdiff & loops$AClfc < 0) | (loops$TCdiff & loops$TClfc < 0),] 
ov = findOverlaps(query = lostLoops,
                  subject = enh)
write.table(x = enh[unique(subjectHits(ov)),], 
            quote = F, col.names = F, row.names = F, sep = "\t",
            file = "./output/motifPrep/lostLoops.bed")

## Static loops -------
staticLoops = loops[loops$cluster == "static",] 
ov = findOverlaps(query = staticLoops,
                  subject = enh)
write.table(x = enh[unique(subjectHits(ov)),], 
            quote = F, col.names = F, row.names = F, sep = "\t",
            file = "./output/motifPrep/bkgd.bed")

