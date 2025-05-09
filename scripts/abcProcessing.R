
# INTIALIZE -----------------------
library(DESeq2)

# READ IN -----------------------

## ABC output -----------
abc = lapply(list.files("./input/abc/",
                        recursive = T, 
                        pattern = "EnhancerPredictionsAllPutative.tsv$", 
                        full.names = T),
             read.table, header = T)
names(abc) = c("A", "T", "C")

## Denylist --------
deny = readRDS("./output/denyList/denylist.rds")

## Read in DEGs  ------
genes = GRanges(readRDS("./output/rnaProcessing/genes.rds"))

## Read in enhancers -------
enh = readRDS("./output/H3K27acProcessing/putativeEnhancers.rds")

## .hic files -------
hicFiles = list.files("/Volumes/MCF10/CantataData/proc/techrep/hicFiles/",
                        full.names = T)
hicFiles = hicFiles[c(1:8, 17:24, 9:16)]


# FILTER ABC OUTPUT -----------------------
## Merge into data.frame ---------
abcMerge = cbind(abc$A[, c(1:5, 11:13, 15:16, 19:22, 25)],
                 abc$A[, -c(1:5, 11:13, 15:16, 19:22, 25)],
                 abc$T[, -c(1:5, 11:13, 15:16, 19:22, 25)],
                 abc$C[, -c(1:5, 11:13, 15:16, 19:22, 25)])
colnames(abcMerge) = c(colnames(abc$A[, c(1:5, 11:13, 15:16, 19:22, 25)]),
                       paste0("A_", colnames(abc$A[, -c(1:5, 11:13, 15:16, 19:22, 25)])),
                       paste0("T_", colnames(abc$T[, -c(1:5, 11:13, 15:16, 19:22, 25)])),
                       paste0("C_", colnames(abc$C[, -c(1:5, 11:13, 15:16, 19:22, 25)])))

## Convert to genomicRanges ------
all2gr <- function(all){
  gi = GRanges(
    data.frame(seqnames = all$chr,
               start = rowMins(cbind(all$start, all$TargetGeneTSS)),
               end = rowMaxs(cbind(all$start, all$TargetGeneTSS)),
               all[, -c(1:3)]))
}

## Remove denylist regions -------
ov = findOverlaps(query = all2gr(abcMerge),
                  subject = deny)
abcMerge = abcMerge[-(unique(queryHits(ov))),] # 11.5M --> 6.8M

## Subsetting for valid EP pairs ---------
thresh = 0.025
abcMerge$valid = ((abcMerge$A_ABC.Score > thresh | abcMerge$T_ABC.Score > thresh | abcMerge$C_ABC.Score > thresh) &
               (abcMerge$A_ABC.Score != 0 & abcMerge$T_ABC.Score != 0 & abcMerge$C_ABC.Score != 0) & 
               (abcMerge$A_ABC.Score != 1 & abcMerge$T_ABC.Score != 1 & abcMerge$C_ABC.Score != 1))

pairs = abcMerge[abcMerge$valid,]

dim(abcMerge) # 6.786M
dim(pairs) # 148k


# ADD E/P DATA -----------------------
## Add in target gene info to ABC data --------
geneDat = mcols(genes)[match(pairs$TargetGeneEnsembl_ID, genes$GENEID),]

pairs = pairs[!is.na(rownames(geneDat)),] # 148449 --> 148161
geneDat = geneDat[!is.na(rownames(geneDat)),]

pairs = cbind(pairs, geneDat)

colnames(pairs)[68:104] = paste0("RNA_", colnames(pairs)[68:104])

## Add enhancer info to ABC -----------
ov = findOverlaps(query = GRanges(pairs),
                  subject = enh)
pairs = cbind(pairs[queryHits(ov),], mcols(enh[subjectHits(ov),]))
colnames(pairs)[105:130] = paste0("ENH_", colnames(pairs)[105:130])

## Add promoter H3K27ac info to ABC -------
## Identify promoter peaks
ov = findOverlaps(query = promoters(genes,
                                    upstream = 2000,
                                    downstream = 500),
                  subject = enh)
proms = enh[subjectHits(ov),]
proms$gene = genes$SYMBOL[queryHits(ov)]

## Aggregate data for all promoter regions per gene
tmp = data.frame(proms)
tmp1 = aggregate(tmp[,6:11], by=list(tmp$gene), FUN = sum)
tmp2 = aggregate(tmp[,c(13:18, 20:24, 26:30)], 
                 by=list(tmp$gene), FUN = mean, na.rm=T)
tmp3 = aggregate(tmp[,c(19, 25, 31)],
                 by=list(tmp$gene), FUN = any)

## Find gene promoters on only 1 chromosome 
chrs = aggregate(tmp$seqnames, by=list(tmp$gene), FUN = list)
chrs = unlist(lapply(lapply(chrs$x, unique), length))

geneProms = cbind(tmp1[chrs == 1,],
                  tmp2[chrs == 1,],
                  tmp3[chrs == 1,])
rownames(geneProms) = geneProms$Group.1

## Format promoter info
geneProms = geneProms[,-c(8, 25)]
geneProms = geneProms[,c(2:13, 24, 
                         14:18, 25,
                         19:23, 26)]
colnames(geneProms) = paste0("PRO_", colnames(geneProms))

## Match promoters to ABC target genes, add to ABC pairs
pairs = cbind(pairs, 
              geneProms[match(pairs$TargetGene, rownames(geneProms)),])

# ADD EP CONTACT DATA -----------------------
## Convert ABC data.frame to GenomicInteractions -----
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
abcgi = abc2gi(pairs)

## Resize to 10kb ---------
abcgi = swapAnchors(abcgi)
abcgi = sort(abcgi)
abcgi = assignToBins(abcgi, binSize = 10000, 
                     pos1="center", pos2="center")

## Remove denylist regions ------ (150442 --> 150340)
ov = findOverlaps(query = abcgi,
                  subject = deny)
abcgi = abcgi[-(unique(queryHits(ov))),]


## Extract counts ---------
# Read in counts from hic files
countMatrix = pullHicPixels(x = abcgi, 
                            files = hicFiles, 
                            binSize = 10000, 
                            norm = "NONE")

## Differential analysis ---------
# Pull count matrix for DESeq
fullCountMatrix=counts(countMatrix)
fullCountMatrix=suppressWarnings(apply(fullCountMatrix, 2, as.numeric))

rownames(fullCountMatrix) = paste0("contact_", 1:nrow(fullCountMatrix))
fullabcgi = abcgi

# Subset
identical(nrow(fullCountMatrix), length(fullabcgi))
set = (rowMeans(fullCountMatrix) < 70) & (fullabcgi$distance >= 10000)
length(set)
abcgi = fullabcgi[set,]
countMatrix = fullCountMatrix[set,]

# Make DNA copy matrix
dnacopy <- matrix(1, nrow=nrow(countMatrix), ncol=ncol(countMatrix))
tmp = as.data.frame(abcgi)

# Note: Based on readouts from NeoLoopFinder CNV
dnacopy[tmp$seqnames1 == "chr1" & tmp$start1 >= 143200000,  1:8] <- 1.32
dnacopy[tmp$seqnames1 == "chr1" & tmp$start1 >= 143200000,  9:16] <- 0.96
dnacopy[tmp$seqnames1 == "chr1" & tmp$start1 >= 143200000,  17:24] <- 0.95

dnacopy[tmp$seqnames1 == "chr1" & tmp$start1 >= 201000000,  1:8] <- 1.57
dnacopy[tmp$seqnames1 == "chr1" & tmp$start1 >= 201000000,  9:16] <- 0.87
dnacopy[tmp$seqnames1 == "chr1" & tmp$start1 >= 201000000,  17:24] <- 0.88

dnacopy[tmp$seqnames1 == "chr3" & tmp$start1 >= 69000000,  1:8] <- 0.95
dnacopy[tmp$seqnames1 == "chr3" & tmp$start1 >= 69000000,  9:16] <- 0.91
dnacopy[tmp$seqnames1 == "chr3" & tmp$start1 >= 69000000,  17:24] <- 1.43

dnacopy[tmp$seqnames1 == "chr8" & tmp$start1 >= 99000000,  1:8] <- 1.48
dnacopy[tmp$seqnames1 == "chr8" & tmp$start1 >= 99000000,  9:16] <- 0.89
dnacopy[tmp$seqnames1 == "chr8" & tmp$start1 >= 99000000,  17:24] <- 0.93

dnacopy[tmp$seqnames1 == "chr9" & tmp$end1 <= 39000000,  1:8] <- 0.81
dnacopy[tmp$seqnames1 == "chr9" & tmp$end1 <= 39000000,  9:16] <- 1.06
dnacopy[tmp$seqnames1 == "chr9" & tmp$end1 <= 39000000,  17:24] <- 1.43

dnacopy[tmp$seqnames1 == "chr9" & tmp$start1 >= 65000000,  1:8] <- 0.93
dnacopy[tmp$seqnames1 == "chr9" & tmp$start1 >= 65000000,  9:16] <- 0.99
dnacopy[tmp$seqnames1 == "chr9" & tmp$start1 >= 65000000,  17:24] <- 1.32

dnacopy[tmp$seqnames1 == "chr10", 1:8] <- 0.97
dnacopy[tmp$seqnames1 == "chr10", 9:16] <- 1.02
dnacopy[tmp$seqnames1 == "chr10", 17:24] <- 1.23

dnacopy[tmp$seqnames1 == "chr13" & tmp$start >= 51000000, 1:8] <- 1.19
dnacopy[tmp$seqnames1 == "chr13" & tmp$start >= 51000000, 9:16] <- 1.00
dnacopy[tmp$seqnames1 == "chr13" & tmp$start >= 51000000, 17:24] <- 1.07

dnacopy[tmp$seqnames1 == "chr19" & tmp$start >= 34000000, 1:8] <- 1.02
dnacopy[tmp$seqnames1 == "chr19" & tmp$start >= 34000000, 9:16] <- 1.22
dnacopy[tmp$seqnames1 == "chr19" & tmp$start >= 34000000, 17:24] <- 0.97

dnacopy[tmp$seqnames1 == "chr20", 1:8] <- 0.96
dnacopy[tmp$seqnames1 == "chr20", 9:16] <- 1.33
dnacopy[tmp$seqnames1 == "chr20", 17:24] <- 0.96

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

# Make list of all significant contacts
L = log(1.5, base=2)
p = 0.025

sigList = lapply(resList, function(res){
  rownames(res)[which(res$padj <= p &
                        abs(res$log2FoldChange) >= L)]
})

allSig = unique(unlist(sigList))
length(allSig)

hist(resList$resAT$pvalue, 
     breaks=seq(0,1,by=0.01),
     col = rep(c("red", "grey"), 
               times = c(sum(seq(0,1,by=0.01) <= p),
                         sum(seq(0,1,by=0.01) > p))),
     main = "Raw p-values")

## Compile contact info --------
# Add additional loop info to DF
contacts = as.data.frame(fullabcgi)
contacts$contactName = rownames(fullCountMatrix)

# Loop spans
contacts$span = contacts$end2 - contacts$start1

# Average + max raw counts
contacts$avgCount = rowMeans(fullCountMatrix)
contacts$maxCount = rowMaxs(fullCountMatrix)

# DESeq info
contacts$ATlfc = resList$resAT$log2FoldChange[match(contacts$contactName,
                                                    rownames(resList$resAT))]
contacts$AClfc = resList$resAC$log2FoldChange[match(contacts$contactName,
                                                    rownames(resList$resAC))]
contacts$TClfc = resList$resTC$log2FoldChange[match(contacts$contactName,
                                                    rownames(resList$resTC))]
contacts$padj = resList$resAT$padj[match(contacts$contactName,
                                         rownames(resList$resAT))]
contacts$ATdiff = contacts$contactName %in% sigList$resAT
contacts$ACdiff = contacts$contactName %in% sigList$resAC
contacts$TCdiff = contacts$contactName %in% sigList$resTC

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
contacts = cbind(contacts,
                 fullCountMatrix,
                 'maxLFC' = maxLFC[match(contacts$contactName,
                              rownames(resList$resAT))],
                 'maxLFCcomp' = maxLFCcomp[match(contacts$contactName,
                                  rownames(resList$resAT))])
diffContacts = contacts[contacts$ATdiff | contacts$ACdiff | contacts$TCdiff,]

# Add to ABC gInteractions object
pairs = pairs[match(paste0(fullabcgi$name, fullabcgi$TargetGene), 
                    paste0(pairs$name, pairs$TargetGene)),]
identical(pairs$TargetGene, contacts$TargetGene)
identical(fullabcgi$TargetGene, contacts$TargetGene)
identical(pairs$name, contacts$name)
identical(fullabcgi$name, contacts$name)

colnames(contacts) = paste0("HIC_", colnames(contacts))

pairs = cbind(pairs,
              contacts[,163:199])


# PLOT: EP distance -----------
pdf("./output/abcProcessing/FigS4A-EPdistanceDistr.pdf",
    width = 6, height = 6)

plot(density(pairs$distance), main = "E-P Distance")
abline(v = mean(pairs$distance), col = "red", lty = 2)
mean(pairs$distance)

dev.off()



# OUTPUT -------------
save(contacts, countMatrix, dds, dds.trans, allSig,
     fullCountMatrix, fullabcgi, abcgi,
     file = "./output/abcProcessing/diffContactObjects.rda")
saveRDS(abcMerge, file = "./output/abcProcessing/abcAll.rds")
saveRDS(pairs, file = "./output/abcProcessing/abcPairs.rds")
