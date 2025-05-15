
# INITIALIZE --------------------------
library(GenomicRanges)


# READ IN --------------------------

## CALDER output (Cell level): Compartments ---------
dir = "./input/compartments/"
compartments = lapply(list.files(path = dir, recursive = T, 
                            pattern = "*all_sub_compartments.tsv", full.names = T),
                 read.table, comment.char = "\\", header = T)
names(compartments) = list.dirs(path = dir, recursive = F, full.names = F)

## Processing CALDER compartments -------
compCols = c("#0000FF",
             "#4848FF",
             "#9191FF",
             "#DADAFF",
             "#FFDADA",
             "#FF9191",
             "#FF4848",
             "#FF0000")
names(compCols) = c("B.2.2", 
                    "B.2.1",
                    "B.1.2",
                    "B.1.1",
                    "A.2.2",
                    "A.2.1",
                    "A.1.2",
                    "A.1.1")

compartments = lapply(compartments, function(c){
  tmp = as.factor(c$comp_rank)
  c$compartment = names(compCols[tmp])
  c$color = compCols[tmp]
  return(GRanges(c))
})

# RUN --------------------------

## CALDER cell type overview -------
pdf("./output/compartments/Fig1C-compPercentages.pdf",
    width = 8, height = 8)

par(mfrow=c(1,1))
par(mar=c(4,4,4,4))

dat = c(unlist(lapply(split(width(compartments$MCF10A), 
                            compartments$MCF10A$compartment), sum)),
        unlist(lapply(split(width(compartments$MCF10AT1), 
                            compartments$MCF10AT1$compartment), sum)),
        unlist(lapply(split(width(compartments$MCF10CA1a), 
                            compartments$MCF10CA1a$compartment), sum)))
dat = matrix(dat, ncol = 3, byrow = F)
dat = t(apply(dat, 1, function(c){c/colSums(dat)}))

barplot(dat, 
        col = colorRampPalette(c("steelblue3", "white", "firebrick3"))(8),
        border = NA,
        names.arg = c("MCF10A", "MCF10AT1", "MCF10CA1a"),
        las = 1,
        yaxt = 'n',
        main = "Compartments")

axis(side = 2, at = seq(0, 1, .25), las = 2, 
     labels = paste0(seq(0, 100, 25), "%"),
     xpd = T, lwd = 0, line = -0.8)
abline(h = seq(0, 1, .125), col = "grey")

dev.off()

## Percent A/B by cell type ---------
colSums(dat[1:4,])/colSums(dat) # %A
colSums(dat[5:8,])/colSums(dat) # %B

## Compartment shifts --------
## Limit to common bins
commonComps = lapply(compartments, function(c){
  ov = subsetByOverlaps(x = c,
                        ranges = compartments$MCF10A)
  ov = subsetByOverlaps(x = ov,
                        ranges = compartments$MCF10AT1)
  ov = subsetByOverlaps(x = ov,
                        ranges = compartments$MCF10CA1a)
})

identical(ranges(commonComps$MCF10A),
          ranges(commonComps$MCF10AT1))
identical(ranges(commonComps$MCF10A),
          ranges(commonComps$MCF10CA1a))


## Changes in compartments -------

## Consolidate compartment calls
comps = commonComps$MCF10A[,4]

ov = findOverlaps(query = comps,
                  subject = commonComps$MCF10AT1)
mcols(comps) = cbind(mcols(comps),
                     mcols(commonComps$MCF10AT1[subjectHits(ov), 4]))

ov = findOverlaps(query = comps,
                  subject = commonComps$MCF10CA1a)
mcols(comps) = cbind(mcols(comps),
                     mcols(commonComps$MCF10CA1a[subjectHits(ov), 4]))

colnames(mcols(comps)) = c("A", "T", "C")


## Find changes in compartments
AvT = table(comps$A, # rows
            comps$T) # cols
sum(diag(AvT))/sum(AvT) # % same bin in A + T
sum(AvT[1:4,5:8], AvT[5:8, 1:4])/sum(AvT) # % A>B and B>A bins
sum(AvT[1:4,5:8])/sum(AvT) # % A>B bins
sum(AvT[5:8, 1:4])/sum(AvT) # % B>A bins
sum(AvT[upper.tri(AvT)])/sum(AvT) # % more B-like in T vs A
sum(AvT[lower.tri(AvT)])/sum(AvT) # % more A-like in T vs A
heatmap(AvT, Rowv = NA, Colv = NA)

TvC = table(comps$T, # rows
            comps$C) # cols
sum(diag(TvC))/sum(TvC) # % same bin in T + C
sum(TvC[1:4,5:8], TvC[5:8, 1:4])/sum(TvC) # % A>B and B>A bins
sum(TvC[1:4,5:8])/sum(TvC)  # % A>B bins
sum(TvC[5:8, 1:4])/sum(TvC) # % B>A bins
sum(TvC[upper.tri(TvC)])/sum(TvC) # % more B-like in C vs T
sum(TvC[lower.tri(TvC)])/sum(TvC) # % more A-like in C vs T
heatmap(TvC, Rowv = NA, Colv = NA)

total = AvT + TvC
sum(diag(total))/sum(total) # % same bin
sum(total[1:4,5:8], total[5:8, 1:4])/sum(total) # % A>B and B>A bins
sum(total[1:4,5:8])/sum(total) # % A>B bins
sum(total[5:8, 1:4])/sum(total) # % B>A bins
sum(total[upper.tri(total)])/sum(total) # % more B-like bins
sum(total[lower.tri(total)])/sum(total) # % more A-like bins
heatmap(total, Rowv = NA, Colv = NA)

AvC = table(comps$A, # rows
            comps$C) # cols
sum(diag(AvC))/sum(AvC) # % same bin
sum(AvC[1:4,5:8], AvC[5:8, 1:4])/sum(AvC) # % A>B and B>A bins
sum(AvC[1:4,5:8])/sum(AvC)  # % A>B bins
sum(AvC[5:8, 1:4])/sum(AvC) # % B>A bins
sum(AvC[upper.tri(AvC)])/sum(AvC) # % more B-like in C vs A
sum(AvC[lower.tri(AvC)])/sum(AvC) # % more A-like in C vs A
heatmap(AvC, Rowv = NA, Colv = NA) 


# OUTPUT --------------------------
saveRDS(comps, "./output/compartments/compartments.RDS")


