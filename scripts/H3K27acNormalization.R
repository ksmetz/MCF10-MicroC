
# INTIALIZE -----------------------
library(rtracklayer)
library(edgeR)

# READ IN -----------------------

## H3K27ac signal tracks -------
sigFiles = list.files("./input/h3k27ac/signal/",
                      pattern = "*H3K27ac_norm.bw",
                      full.names = T)
names(sigFiles) = gsub("_H3K27ac_norm.bw",
                       "",
                       basename(sigFiles))

## Read in enhancers -------
enh = GRanges(readRDS(file = "./output/H3K27acProcessing/enhancers.rds"))

enh = resize(enh, 5000, fix = "center")

# EXTRACT COUNTS ----------------------
count <- c()
for (n in 1:length(sigFiles)) {
  meancov <- unlist(summary(BigWigFile(sigFiles[n]), 
                            enh, size = 1, type = "mean"))
  
  score = meancov$score*width(enh)
  
  count = cbind(count, score)
}
rownames(count) <- colnames(count) <- NULL
colnames(count) = names(sigFiles)


# NORM FACTORS ----------------------
normFactor = calcNormFactors(count, method = "TMM")
libSize <- colSums(count)
libSize = libSize/max(libSize)
normFactors <- normFactor * 1/libSize