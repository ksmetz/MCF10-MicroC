
# INITIALIZE --------------------------
library(ComplexHeatmap)


# READ IN --------------------------

## Read in compartments ------
EVfiles = file.path("./input/eigenvectors/")

brs = c("MCF10_A", "MCF10_T", "MCF10_C")
brs = paste0(brs, "_eigen_100kb.bed")
EV = lapply(file.path(EVfiles, brs), read.table)
EV = cbind(EV[[1]][,c(1:3,5)],
           EV[[2]][,5],
           EV[[3]][,5])
colnames(EV) = c("chr", "start", "end",
                 "MCF10_A",
                 "MCF10_T",
                 "MCF10_C")

## .hic files -------
hicFiles <- list.files(path = "/Volumes/MCF10/CantataData/proc/celltype/hicFiles/",
                       full.names = T)
names(hicFiles) = c("A",
                    "C",
                    "T")
hicFiles = hicFiles[c(1, 3, 2)]


## Denylist -------------
deny = readRDS("./output/denyList/denylist.rds")

# FILTER --------------------------

## Remove denylist regions --------
ov = findOverlaps(query = GRanges(EV),
                  subject = deny)
EV = EV[-unique(queryHits(ov)),]

## Subset for specific chromosomes -------------
chrom = c("chr2", "chr12", "chr17")
EV = EV[EV$chr %in% chrom,]
EVs = split(EV, EV$chr)


# BIN EVs --------------------------
EVs = lapply(EVs, function(EV){
  
  cuts = cut(EV$MCF10_A, 
             unique(quantile(EV$MCF10_A, seq(0, 1, inc))), 
             include.lowest = T)
  levels(cuts) = 1:(1/inc)
  EV$Abucket = cuts
  
  cuts = cut(EV$MCF10_T, 
             unique(quantile(EV$MCF10_T, seq(0, 1, inc))), 
             include.lowest = T)
  levels(cuts) = 1:(1/inc)
  EV$Tbucket = cuts
  
  cuts = cut(EV$MCF10_C, 
             unique(quantile(EV$MCF10_C, seq(0, 1, inc))), 
             include.lowest = T)
  levels(cuts) = 1:(1/inc)
  EV$Cbucket = cuts
  
  return(EV)
})
EV = do.call(rbind, EVs)



# EXTRACT COUNTS ---------
saddleMat <- function(binCol, hicFile){
  mat = list()
  for(n in 1:(1/inc)){
    b1 = levels(EV[,binCol])[n]
    print(b1)
    
    if(n != 1){
      b2list = levels(EV[,binCol])[-(1:n-1)]
    }else{
      b2list = levels(EV[,binCol])
    }
    
    dat = c()
    i = 0
    for(b2 in b2list){
      i = i+1
      regions = c()
      for(c in unique(EV$chr)){
        r = expand.grid(which(EV[,binCol] == b1 & EV$chr == c),
                        which(EV[,binCol] == b2 & EV$chr == c))
        regions[[c]] = r
      }
      regions = do.call(rbind, regions)
      regions = as_ginteractions(cbind(EV[regions$Var1, 1:3],
                                       EV[regions$Var2, 1:3]))
      regions = assignToBins(regions, binSize = 100000, 
                             pos1="center", pos2="center")
      
      counts = pullHicPixels(x = regions, 
                             files = hicFile, 
                             binSize = 100000, 
                             matrix = "oe",
                             norm = "SCALE")
      val = mean(counts(counts), na.rm=T)
      dat[i] = val
    }
    
    mat[[n]] = dat
  }
  
  ## Prep matrix ----------
  tmp = matrix(data = NA, nrow = 1/inc, ncol = 1/inc)
  tmp[lower.tri(tmp, diag = T)] = unlist(mat)
  
  tmp[upper.tri(tmp)] = t(tmp)[upper.tri(tmp)]
  rownames(tmp) = levels(EV$Abucket)
  colnames(tmp) = levels(EV$Abucket)
  
  return(tmp)
}

Amat = saddleMat(binCol = "Abucket", hicFile = hicFiles[1])
Tmat = saddleMat(binCol = "Tbucket", hicFile = hicFiles[2])
Cmat = saddleMat(binCol = "Cbucket", hicFile = hicFiles[3])

# PLOT: Saddleplots ---------
pdf("./output/EVsaddleplots/Fig1F-saddleplots.pdf",
    width = 8, height = 8)

col_fun = colorRamp2(c(-0.75, 0, 0.75), c("steelblue3", "white", "firebrick3"))

Heatmap(log(Amat[(1/inc):1, 1:(1/inc)]), 
        cluster_rows = F, cluster_columns = F, 
        col = col_fun,
        name = "log(O/E)",
        column_title = "MCF10A",
        # column_labels = c("inactive", rep("", times = 18), "active"),
        column_names_rot = 0,
        column_names_side = "bottom",
        # row_labels = c("active", rep("", times = 18), "inactive"),
        row_names_side = "left")
lapply(1:(1/inc), 
       function(n){mean(EV$MCF10_A[EV$Abucket == n], na.rm = T)}) |> unlist() |> 
  barplot(border = NA, las = 2)


Heatmap(log(Tmat[(1/inc):1, 1:(1/inc)]), 
        cluster_rows = F, cluster_columns = F, 
        col = col_fun,
        name = "log(O/E)",
        column_title = "MCF10AT1",
        # column_labels = c("inactive", rep("", times = 18), "active"),
        column_names_rot = 0,
        column_names_side = "bottom",
        # row_labels = c("active", rep("", times = 18), "inactive"),
        row_names_side = "left")
lapply(1:(1/inc), 
       function(n){mean(EV$MCF10_T[EV$Tbucket == n], na.rm = T)}) |> unlist() |> 
  barplot(border = NA, las = 2)


Heatmap(log(Cmat[(1/inc):1, 1:(1/inc)]), 
        cluster_rows = F, cluster_columns = F, 
        col = col_fun,
        name = "log(O/E)",
        column_title = "MCF10CA1a",
        # column_labels = c("inactive", rep("", times = 18), "active"),
        column_names_rot = 0,
        column_names_side = "bottom",
        # row_labels = c("active", rep("", times = 18), "inactive"),
        row_names_side = "left")
lapply(1:(1/inc), 
       function(n){mean(EV$MCF10_C[EV$Cbucket == n], na.rm = T)}) |> unlist() |> 
  barplot(border = NA, las = 2)

dev.off()



# OUTPUT ---------
save(Amat, Tmat, Cmat,
     EV,
     file = "./output/EVsaddleplots/saddleplots.rda")
