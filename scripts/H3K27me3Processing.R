
# INTIALIZE -----------------------

# READ IN -----------------------

## H3K27me3 peaks ----------
k27peaks = list.files("./input/h3k27me3/peaks/", full.names = T)

## H3K27me3 alignments ----------
k27BAMs = list.files("./input/h3k27me3/alignments/", full.names = T,
                     pattern = ".*bam$")

## ENCODE denylist -------
library(BSgenome.Hsapiens.UCSC.hg38)
# https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gzdenylist = "~/Tools/reference/hg"
deny = import("~/Tools/reference/hg38-blacklist.v2.bed.gz",
              seqinfo=seqinfo(BSgenome.Hsapiens.UCSC.hg38))


# FILTER PEAKS ----------------
## Consensus peak list ------
k27peaks = lapply(k27peaks, read.table, 
                  col.names = c("chr", 'start', 'stop', 
                                'name', "score", "strand",
                                "signalValue", "pValue", "qValue", "summit"))
k27peaks = lapply(k27peaks, GRanges)

k27peaks = do.call("c", k27peaks) # 294259; 170-20490 bp
k27peaks = reduce(k27peaks) # 206584; 170-27680 bp

## Format denylist ------
readLength <- 50
deny <- trim(resize(deny,
                    width=width(deny) + 2*readLength, 
                    fix="center"))

# 206k -> 202k regions
k27peaks = k27peaks[k27peaks %outside% deny,]


# EXTRACT COUNTS ----------------
library(bedtoolsr)
k27Counts = bt.multicov(bams = paste0(k27BAMs, collapse = " "),
                        bed = as.data.frame(k27peaks)) # 3:21-3:31
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

# DIFFERENTIAL PEAKS ----------------
rownames(k27Counts) = paste0("putativeRep_", 1:nrow(k27Counts))
countMat = k27Counts[, -c(1:5)]

## Run DESeq on count matrix -----
colData = data.frame(
  cell = rep(c("MCF10A", "MCF10AT1", "MCF10CA1a"), each = 2),
  rep = as.factor(rep(c(1, 2), times = 3))
)
rownames(colData) = colnames(countMat)

dds = DESeqDataSetFromMatrix(
  countData = countMat, 
  colData = colData, 
  design = ~ cell) 

dds <- DESeq(dds)

## Find results ----------
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

identical(rownames(k27Counts), rownames(diffInfo))
repr = cbind(k27Counts, diffInfo)

# OUTPUT -----------------
saveRDS(repr, file = "./output/H3K27me3Processing/repressors.rds")

