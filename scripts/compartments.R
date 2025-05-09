
# READ IN --------------------------

## CALDER output (Cell level): Compartments ---------
dir = "./input/compartments/"
comps = lapply(list.files(path = dir, recursive = T, 
                            pattern = "*all_sub_compartments.bed", full.names = T),
                 read.table, comment.char = "\\")
names(comps) = list.dirs(path = dir, recursive = F, full.names = F)

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

comps = lapply(comps, function(c){
  tmp = factor(c$V9,
               levels = compCols)
  levels(tmp) = names(compCols)
  c$V9 = tmp
  return(c)
})


# RUN --------------------------

## CALDER cell type overview -------
pdf("./output/compartments/Fig1C-compPercentages.pdf",
    width = 8, height = 8)

par(mfrow=c(1,1))
par(mar=c(4,4,4,4))

dat = c(unlist(lapply(split(comps$MCF10A$V3 - comps$MCF10A$V2, comps$MCF10A$V9), sum)),
        unlist(lapply(split(comps$MCF10AT1$V3 - comps$MCF10AT1$V2, comps$MCF10AT1$V9), sum)),
        unlist(lapply(split(comps$MCF10CA1a$V3 - comps$MCF10CA1a$V2, comps$MCF10CA1a$V9), sum)))
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
colSums(dat[1:4,])/colSums(dat) # %B
colSums(dat[5:8,])/colSums(dat) # %A


# OUTPUT --------------------------
saveRDS(comps, "./output/compartments/compartments.RDS")


