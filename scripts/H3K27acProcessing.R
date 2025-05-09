
# INTIALIZE -----------------------
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(DESeq2)


# READ IN -----------------------

## ABC enhancer regions w/ H3K27ac counts -------
enh = readRDS("./output/atacProcessing/K27acCounts.rds")
enh$name = paste0("putativeEnh_", 1:nrow(enh))
enh = GRanges(enh)

## ENCODE denylist -------
# https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gzdenylist = "~/Tools/reference/hg"
deny = import("~/Tools/reference/hg38-blacklist.v2.bed.gz",
              seqinfo=seqinfo(BSgenome.Hsapiens.UCSC.hg38))

## H3K27ac peaks -------
k27peaks = list.files("./input/h3k27ac/peaks/", full.names = T)


# FILTER PEAKS -----------------------

## Format denylist ------
readLength <- 50
deny <- trim(resize(deny,
                    width=width(deny) + 2*readLength, 
                    fix="center"))

# 199k -> 182k regions
enh = enh[enh %outside% deny,]
enhMat = as.data.frame(enh)[,-c(1:5, 12)]
rownames(enhMat) = enh$name


# DIFFERENTIAL PEAKS -----------------------
## Run DESeq on count matrix -----
colData = data.frame(
  cell = rep(c("MCF10A", "MCF10AT1", "MCF10CA1a"), each = 2),
  rep = as.factor(rep(c(1, 2), times = 3))
)
rownames(colData) = colnames(enhMat)

dds = DESeqDataSetFromMatrix(
  countData = enhMat, 
  colData = colData, 
  design = ~ cell) 

dds <- DESeq(dds)

## Find results
resAC <- results(dds, contrast = c("cell", "MCF10CA1a", "MCF10A"),
                 independentFiltering = F)
resAT <- results(dds, contrast = c("cell", "MCF10AT1", "MCF10A"),
                 independentFiltering = F)
resTC <- results(dds, contrast = c("cell", "MCF10CA1a", "MCF10AT1"),
                 independentFiltering = F)

## Pull significantly differential peaks --------
findSig <- function(res, pThresh, lThresh){
  dat = res[!is.na(res$padj),]
  dat = dat[dat$padj <= pThresh & abs(dat$log2FoldChange) >= lThresh,]
  return(dat)
}

p = 0.05
L = 0
sigAC <- findSig(resAC, p, L) 
sigAT <- findSig(resAT, p, L) 
sigTC <- findSig(resTC, p, L) 

dim(sigAC)
dim(sigAT)
dim(sigTC)

## Compile into data.frame --------
diffInfo = cbind(resAC, 
                 AC.sig = rownames(resAC) %in% rownames(sigAC),
                 resAT[,-1],
                 AT.sig = rownames(resAT) %in% rownames(sigAT),
                 resTC[,-1], 
                 TC.sig = rownames(resAC) %in% rownames(sigTC))
colnames(diffInfo)[2:6] = paste0("AC.", colnames(diffInfo)[2:6])
colnames(diffInfo)[8:12] = paste0("AT.", colnames(diffInfo)[8:12])
colnames(diffInfo)[14:18] = paste0("TC.", colnames(diffInfo)[14:18])

identical(enh$name, rownames(diffInfo))
mcols(enh) = cbind(mcols(enh), diffInfo)

## Overlap with H3K27ac peaks ------
## Consensus peak list 
k27peaks = lapply(k27peaks, import, format = "broadPeak")
k27peaks = unlist(as(k27peaks, "GRangesList"))
k27peaks = reduce(k27peaks)

## Overlap enhancers (ATAC peaks) with H3K27ac
ov = findOverlaps(query = k27peaks,
                  subject = enh)
enhOverlap = enh[subjectHits(ov)] # 182k --> 52.9k


# OUTPUT --------------------
saveRDS(k27peaks, file = "./output/H3K27acProcessing/H3K27acPeaks.rds")
saveRDS(enhOverlap, file = "./output/H3K27acProcessing/enhancers.rds")
saveRDS(enh, file = "./output/H3K27acProcessing/putativeEnhancers.rds")

