
## Note:
## First ran juicer_tools dump using the 'mega' + cell type map .hic files
## extracting SCALE norms for each chromosome @ 5kb + 10kb:
##    for c in {1..22} X; 
##    do echo chr${c}; 
##    java -jar ~/Tools/juicer_tools_1.22.01.jar dump norm SCALE /Volumes/MCF10/CantataData/proc/mega/MCF10_merged_map.hic chr${c} BP 5000 ./chr${c}_5kb_SCALE-norms.txt; 
##    done

# INITIALIZE -----------------------------------------------------------------------------

library(readr)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

txdb = TxDb.Hsapiens.UCSC.hg38.knownGene


# FUNCTIONS    -----------------------------------------------------------------------------

## Function for making a deny list from a norm file at a given resolution
makeDenyList <- function(chroms = paste0("chr", c(1:22, "X")), res = 5000){
  
  # Read in normalization tables
  megaNorms = lapply(chroms, function(c){
    read.table(
      list.files(
        paste0("./input/norms/", res/1000, "kb/mega"), 
        pattern = paste0("^", c, "_"),
        full.names = T))
  })
  aNorms = lapply(chroms, function(c){
    read.table(
      list.files(
        paste0("./input/norms/", res/1000, "kb/A"), 
        pattern = paste0("^", c, "_"),
        full.names = T))
  })
  tNorms = lapply(chroms, function(c){
    read.table(
      list.files(
        paste0("./input/norms/", res/1000, "kb/T"), 
        pattern = paste0("^", c, "_"),
        full.names = T))
  })
  cNorms = lapply(chroms, function(c){
    read.table(
      list.files(
        paste0("./input/norms/", res/1000, "kb/C"), 
        pattern = paste0("^", c, "_"),
        full.names = T))
  })
  
  # Combine normalizations for each chromosome
  norms = Map(function(m, a, t, c){cbind(m, a, t, c)},
              m = megaNorms, a = aNorms, t = tNorms, c = cNorms)
  
  # Find sparse regions
  denyList = lapply(norms, function(n){
    p05 = 0.15
    p99 = 2
    naCheck = is.na(n)
    loCheck = t(apply(n, 1, function(r){r < p05}))
    hiCheck = t(apply(n, 1, function(r){r > p99}))
    check = naCheck | loCheck #| hiCheck
    dat = as.data.frame(apply(check, 1, any)) #| (n > p99)
    return(dat)
  })
  
  # Add coordinates
  allBins = lapply(chroms, function(c){
    seqlevels(txdb) <- c
    bins = tileGenome(seqlengths = seqinfo(txdb),
                      tilewidth = res,
                      cut.last.tile.in.chrom = TRUE)
    bins = as.data.frame(bins)
    bins = bins[, c("seqnames", "start", "end")]
    return(bins)
  })
  
  # Concatenate
  allBins = do.call(rbind, allBins)
  denyList = do.call(rbind, denyList)
  
  # Format
  denyList = cbind(allBins, denyList)
  colnames(denyList) = c("chr", "start", "stop", "denied")
  
  # Save
  return(denyList)
}


# RUN    -----------------------------------------------------------------------------

## Micro-C denylist ---------
# Make denylists at desired resolution
normDenyList = makeDenyList()
normDenyList = GRanges(normDenyList[normDenyList$denied == T,])

# Ignore single bins
normDenyList = reduce(normDenyList, min.gapwidth = 50000)
normDenyList = normDenyList[width(normDenyList) > 5000,]

# Merge areas within 100kb
normDenyList = reduce(normDenyList, min.gapwidth = 100000)

## ENCODE denylist -------
# https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gzdenylist = "~/Tools/reference/hg"
library(BSgenome.Hsapiens.UCSC.hg38)
encodeDenyList = rtracklayer::import("~/Tools/reference/hg38-blacklist.v2.bed.gz",
                                     seqinfo=seqinfo(BSgenome.Hsapiens.UCSC.hg38))

# Format denylist 
readLength <- 50
encodeDenyList <- trim(resize(encodeDenyList,
                              width=width(encodeDenyList) + 2*readLength, 
                              fix="center"))

## Combine denylists ------
denyList = c(normDenyList, encodeDenyList)
denyList = reduce(denyList)
denyList = trim(denyList)


# OUTPUT      -------------------------------------------------------
# Save as RDS files
saveRDS(denyList, file="./output/denyList/denylist.rds")
