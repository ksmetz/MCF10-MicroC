
# INITIALIZE  ---------------------------------------------
library(SpectralTAD)
library(dplyr)
library(strawr)


# READ IN     ---------------------------------------------

## .hic files ---------
hicFiles <- list.files(path = "/Volumes/MCF10/CantataData/proc/celltype/hicFiles/",
                       full.names = T)
names(hicFiles) = c("A",
                    "C",
                    "T")
hicFiles = hicFiles[c(1, 3, 2)]


# CALL TADs     ---------------------------------------------

chrm = paste0('chr', c(1:22, "X"))

## Read in data with strawr ----------
dat = lapply(hicFiles, function(hic){
  mat = lapply(chrm, function(c){
    straw(norm = "SCALE", 
          fname = hic, 
          chr1loc = c, 
          chr2loc = c, 
          unit = 'BP', 
          binsize = 10000, 
          matrix = "observed")
  })
})

## Call TADs ----------
tads = lapply(dat, function(d){
  Map(function(c, d){
    out = SpectralTAD(d,
                      chr = c,
                      resolution = 10000,
                      levels = 3,
                      qual_filter = T,
                      z_clust = F,
                      window_size = 300)
    out = bind_rows(out)
    return(out)
  }, c=chrm, d=d)
})

## Save ----------
saveRDS(tads, file = "./output/tadCalling/tads.RDS")
