
# INTIALIZE -----------------------
library(mariner)
library(InteractionSet)
library(GenomicRanges)
library(nullranges)


# READ IN -----------------------

## Read in ABC pairs ---------
abc = readRDS("./output/abcProcessing/abcPairs.rds")

## Read in DEGs  ------
genes = readRDS("./output/rnaProcessing/genes.rds")

## Read in diff enhancers -------
enh = GRanges(readRDS("./output/H3K27acProcessing/enhancers.rds"))

## Read in loops -------
loops = as_ginteractions(readRDS("./output/loopProcessing/loops.rds"))

# BOXPLOTS -----------------------
setBoxplots <- function(set, description){
  
  par(mfrow=c(1,4))
  
  description = paste0(description,
                       "\nn=", nrow(set))
  
  dat = list("MCF10A" = rowSums(set[, c("ENH_MCF10A_R1", "ENH_MCF10A_R2")]),
             "MCF10CA1a" = rowSums(set[, c("ENH_MCF10CA1a_R1", "ENH_MCF10CA1a_R2")]))
  boxplot(dat,
          main = paste0("Activity at ", description), outline = F,
          col = c("pink"))
  ENHstats = t.test(x = dat$MCF10A,
                    y = dat$MCF10CA1a)
  
  dat = list("MCF10A" = rowSums(set[, c("HIC_MCF10_A_1_1_contact_map.hic", "HIC_MCF10_A_1_2_contact_map.hic",
                                        "HIC_MCF10_A_1_3_contact_map.hic", "HIC_MCF10_A_1_4_contact_map.hic",
                                        "HIC_MCF10_A_2_1_contact_map.hic", "HIC_MCF10_A_2_2_contact_map.hic",
                                        "HIC_MCF10_A_2_3_contact_map.hic", "HIC_MCF10_A_2_4_contact_map.hic")]),
             "MCF10CA1a" = rowSums(set[, c("HIC_MCF10_C_1_1_contact_map.hic", "HIC_MCF10_C_1_2_contact_map.hic",
                                           "HIC_MCF10_C_1_3_contact_map.hic", "HIC_MCF10_C_1_4_contact_map.hic",
                                           "HIC_MCF10_C_2_1_contact_map.hic", "HIC_MCF10_C_2_2_contact_map.hic",
                                           "HIC_MCF10_C_2_3_contact_map.hic", "HIC_MCF10_C_2_4_contact_map.hic")]))
  boxplot(dat,
          main = paste0("Contact at ", description), outline = F,
          col = c("lightblue"))
  HICstats = t.test(x = dat$MCF10A,
                    y = dat$MCF10CA1a)
  
  boxplot(list("MCF10A" = set$A_ABC.Score,
               "MCF10CA1a" = set$C_ABC.Score),
          main = paste0("ABC Score at ", description), outline = F,
          col = c("plum3"))
  ABCstats = t.test(x = set$A_ABC.Score,
                    y = set$C_ABC.Score)
  
  boxplot(list(set$RNA_AC.log2FoldChange),
          main = paste0("Gene LFC at ", description), outline = F,
          col = c("gold2"))
  RNAstats = t.test(x = set$RNA_AC.log2FoldChange)
  
  abline(h=0, lty = 2)
  
  return(list("ENH" = ENHstats,
              "HIC" = HICstats,
              "ABC" = ABCstats,
              "RNA" = RNAstats,
              "geneLFC" = set$RNA_AC.log2FoldChange))
  
}

## PLOT: ABC Boxplots (AvC) -------
pdf("./output/abcBoxplots/Fig3AE-abcBoxplots.pdf",
    width = 10, height = 3)

## Changes in enhancers
set = abc[abc$ENH_AC.sig & abc$ENH_AC.log2FoldChange > 0,]
gainedE = setBoxplots(set = set, description = "Gained enhancers")

## Changes in contact
set = abc[abc$HIC_ACdiff & abc$HIC_AClfc > 0,]
gainedC = setBoxplots(set = set, description = "Gained contact")

## Staticly looped pairs
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
ov = findOverlaps(query = abc2gi(abc[abc$ENH_AC.sig & abc$ENH_AC.log2FoldChange > 0,]),
                  subject = loops)
set = unique(abc[abc$ENH_AC.sig & abc$ENH_AC.log2FoldChange > 0,][queryHits(ov),])
gainedLoop = setBoxplots(set = set, description = "Gained looped enh")

## For matching
set = unique(abc[abc$ENH_AC.sig & abc$ENH_AC.log2FoldChange > 0,])
ov = findOverlaps(query = abc2gi(set),
                  subject = loops)
looped = unique(set[queryHits(ov),])
unloop = unique(set[-queryHits(ov),])

## Contact matched unlooped pairs
set.seed(123)
contMatch = matchRanges(focal = looped,
                     pool = unloop,
                     covar = ~HIC_avgCount) 
set = as.data.frame(contMatch)
contMatchEnh = setBoxplots(set = set, description = "Unlooped Gained Enh")

## Distance matched unlooped pairs
set.seed(123)
distMatch = matchRanges(focal = looped,
                     pool = unloop,
                     covar = ~HIC_span) 
set = as.data.frame(distMatch)
distMatchEnh = setBoxplots(set = set, description = "Unlooped Gained Enh")

dev.off()


## PLOT: Gene LFC Boxplots -------
pdf("./output/abcBoxplots/Fig3H-geneLFCboxplots.pdf",
    width = 6, height = 6)

par(mfrow=c(1,1))
dat = list("diff E" = gainedE$geneLFC,
           "diff C" = gainedC$geneLFC,
           "looped" = gainedLoop$geneLFC,
           "contact" = contMatchEnh$geneLFC,
           "distance" = distMatchEnh$geneLFC)
boxplot(dat, col = "gold", outline = F)
abline(h=0, lty=2)
t.test(dat$`diff E`, dat$`diff C`)
t.test(dat$`diff E`, dat$looped)
t.test(dat$`diff C`, dat$looped)
t.test(dat$looped, dat$contact)
t.test(dat$looped, dat$distance)
t.test(dat$distance, dat$contact)

dev.off()



# MATCHING DISTRIBUTIONS --------
## PLOT: Contact frequency distributions -----
pdf("./output/abcBoxplots/Fig3F-contactDistr.pdf",
    width = 6, height = 6)

plot(density(looped$HIC_avgCount), main="Contact Distribution", 
     xlim=c(0, 40), ylim = c(0, 0.17),
     xaxt='n', bty='n', axes=F, xlab="")
lines(density(unloop$HIC_avgCount), lty=2)

polygon(density(contMatch$HIC_avgCount), 
        col = alpha("steelblue", .5), border=NA)
polygon(density(distMatch$HIC_avgCount), 
        col = alpha("darkgrey", .5), border=NA)

axis(2, lwd=0, lwd.ticks=1, las=2)
axis(side=1, at=seq(0, 40, by=10))

legend("topright",
       legend = c("All E-P pairs",
                  "Looped E-P pairs (C)",
                  "Contact-matched pairs (D)",
                  "Distance-matched pairs (E)"),
       bty = 'n', lty = c(1, 2, NA, NA),
       text.col = c("black", "black", "steelblue", "darkgrey"))

dev.off()

## PLOT: EP distance distributions -----
pdf("./output/abcBoxplots/Fig3G-distanceDistr.pdf",
    width = 6, height = 6)

plot(density(looped$HIC_span), main="Contact Distribution", 
     xlim=c(0, 1e6), ylim = c(0, 1.7e-5),
     xaxt='n', bty='n', axes=F, xlab="")
lines(density(unloop$HIC_span), lty=2)

polygon(density(contMatch$HIC_span), 
        col = alpha("steelblue", .5), border=NA)
polygon(density(distMatch$HIC_span), 
        col = alpha("darkgrey", .5), border=NA)

axis(2, lwd=0, lwd.ticks=1, las=2)
axis(side=1, at=seq(0, 1e6, by=2e5))

legend("topright",
       legend = c("All E-P pairs",
                  "Looped E-P pairs (C)",
                  "Contact-matched pairs (D)",
                  "Distance-matched pairs (E)"),
       bty = 'n', lty = c(1, 2, NA, NA),
       text.col = c("black", "black", "steelblue", "darkgrey"))

dev.off()


# OTHER COMPARISONS -----------
setBoxplotsAT <- function(set, description){
  
  par(mfrow=c(1,4))
  
  description = paste0(description,
                       "\nn=", nrow(set))
  
  dat = list("MCF10A" = rowSums(set[, c("ENH_MCF10A_R1", "ENH_MCF10A_R2")]),
             "MCF10AT1" = rowSums(set[, c("ENH_MCF10AT1_R1", "ENH_MCF10AT1_R2")]))
  boxplot(dat,
          main = paste0("Activity at ", description), outline = F,
          col = c("pink"))
  ENHstats = t.test(x = dat$MCF10A,
                    y = dat$MCF10AT1)
  
  dat = list("MCF10A" = rowSums(set[, c("HIC_MCF10_A_1_1_contact_map.hic", "HIC_MCF10_A_1_2_contact_map.hic",
                                        "HIC_MCF10_A_1_3_contact_map.hic", "HIC_MCF10_A_1_4_contact_map.hic",
                                        "HIC_MCF10_A_2_1_contact_map.hic", "HIC_MCF10_A_2_2_contact_map.hic",
                                        "HIC_MCF10_A_2_3_contact_map.hic", "HIC_MCF10_A_2_4_contact_map.hic")]),
             "MCF10AT1" = rowSums(set[, c("HIC_MCF10_T_1_1_contact_map.hic", "HIC_MCF10_T_1_2_contact_map.hic",
                                          "HIC_MCF10_T_1_3_contact_map.hic", "HIC_MCF10_T_1_4_contact_map.hic",
                                          "HIC_MCF10_T_2_1_contact_map.hic", "HIC_MCF10_T_2_2_contact_map.hic",
                                          "HIC_MCF10_T_2_3_contact_map.hic", "HIC_MCF10_T_2_4_contact_map.hic")]))
  boxplot(dat,
          main = paste0("Contact at ", description), outline = F,
          col = c("lightblue"))
  HICstats = t.test(x = dat$MCF10A,
                    y = dat$MCF10AT1)
  
  boxplot(list("MCF10A" = set$A_ABC.Score,
               "MCF10AT1" = set$T_ABC.Score),
          main = paste0("ABC Score at ", description), outline = F,
          col = c("plum3"))
  ABCstats = t.test(x = set$A_ABC.Score,
                    y = set$T_ABC.Score)
  
  boxplot(list(set$RNA_AT.log2FoldChange),
          main = paste0("Gene LFC at ", description), outline = F,
          col = c("gold2"))
  RNAstats = t.test(x = set$RNA_AT.log2FoldChange)
  
  abline(h=0, lty = 2)
  
  return(list("ENH" = ENHstats,
              "HIC" = HICstats,
              "ABC" = ABCstats,
              "RNA" = RNAstats,
              "geneLFC" = set$RNA_AT.log2FoldChange))
  
}

## PLOT: ABC Boxplots (AvT) -------
pdf("./output/abcBoxplots/FigS4B-abcAvTboxplots.pdf",
    width = 10, height = 3)

## Changes in enhancers
set = abc[abc$ENH_AT.sig & abc$ENH_AT.log2FoldChange > 0,]
gainedE = setBoxplotsAT(set = set, description = "Gained enhancers")

## Changes in contact
set = abc[abc$HIC_ATdiff & abc$HIC_ATlfc > 0,]
gainedC = setBoxplotsAT(set = set, description = "Gained contact")

## Staticly looped pairs
ov = findOverlaps(query = abc2gi(abc[abc$ENH_AT.sig & abc$ENH_AT.log2FoldChange > 0,]),
                  subject = loops)
set = unique(abc[abc$ENH_AT.sig & abc$ENH_AT.log2FoldChange > 0,][queryHits(ov),])
gainedLoop = setBoxplotsAT(set = set, description = "Gained looped enh")

## For matching
set = unique(abc[abc$ENH_AT.sig & abc$ENH_AT.log2FoldChange > 0,])
ov = findOverlaps(query = abc2gi(set),
                  subject = loops)
looped = unique(set[queryHits(ov),])
unloop = unique(set[-queryHits(ov),])

## Contact matched unlooped pairs
set.seed(123)
contMatch = matchRanges(focal = looped,
                        pool = unloop,
                        covar = ~HIC_avgCount) 
set = as.data.frame(contMatch)
contMatchEnh = setBoxplotsAT(set = set, description = "Unlooped Gained Enh")

## Distance matched unlooped pairs
set.seed(123)
distMatch = matchRanges(focal = looped,
                        pool = unloop,
                        covar = ~HIC_span) 
set = as.data.frame(distMatch)
distMatchEnh = setBoxplotsAT(set = set, description = "Unlooped Gained Enh")

dev.off()


# OTHER COMPARISONS -----------
setBoxplotsTC <- function(set, description){
  
  par(mfrow=c(1,4))
  
  description = paste0(description,
                       "\nn=", nrow(set))
  
  dat = list("MCF10AT1" = rowSums(set[, c("ENH_MCF10AT1_R1", "ENH_MCF10AT1_R2")]),
             "MCF10CA1a" = rowSums(set[, c("ENH_MCF10CA1a_R1", "ENH_MCF10CA1a_R2")]))
  boxplot(dat,
          main = paste0("Activity at ", description), outline = F,
          col = c("pink"))
  ENHstats = t.test(x = dat$MCF10AT1,
                    y = dat$MCF10CA1a)
  
  dat = list("MCF10AT1" = rowSums(set[, c("HIC_MCF10_T_1_1_contact_map.hic", "HIC_MCF10_T_1_2_contact_map.hic",
                                          "HIC_MCF10_T_1_3_contact_map.hic", "HIC_MCF10_T_1_4_contact_map.hic",
                                          "HIC_MCF10_T_2_1_contact_map.hic", "HIC_MCF10_T_2_2_contact_map.hic",
                                          "HIC_MCF10_T_2_3_contact_map.hic", "HIC_MCF10_T_2_4_contact_map.hic")]),
             "MCF10CA1a" = rowSums(set[, c("HIC_MCF10_C_1_1_contact_map.hic", "HIC_MCF10_C_1_2_contact_map.hic",
                                           "HIC_MCF10_C_1_3_contact_map.hic", "HIC_MCF10_C_1_4_contact_map.hic",
                                           "HIC_MCF10_C_2_1_contact_map.hic", "HIC_MCF10_C_2_2_contact_map.hic",
                                           "HIC_MCF10_C_2_3_contact_map.hic", "HIC_MCF10_C_2_4_contact_map.hic")]))
  boxplot(dat,
          main = paste0("Contact at ", description), outline = F,
          col = c("lightblue"))
  HICstats = t.test(x = dat$MCF10AT1,
                    y = dat$MCF10CA1a)
  
  boxplot(list("MCF10AT1" = set$T_ABC.Score,
               "MCF10CA1a" = set$C_ABC.Score),
          main = paste0("ABC Score at ", description), outline = F,
          col = c("plum3"))
  ABCstats = t.test(x = set$T_ABC.Score,
                    y = set$C_ABC.Score)
  
  boxplot(list(set$RNA_TC.log2FoldChange),
          main = paste0("Gene LFC at ", description), outline = F,
          col = c("gold2"))
  RNAstats = t.test(x = set$RNA_TC.log2FoldChange)
  
  abline(h=0, lty = 2)
  
  return(list("ENH" = ENHstats,
              "HIC" = HICstats,
              "ABC" = ABCstats,
              "RNA" = RNAstats,
              "geneLFC" = set$RNA_TC.log2FoldChange))
  
}

## PLOT: ABC Boxplots (TvC) -------
pdf("./output/abcBoxplots/FigS4B-abcTvCboxplots.pdf",
    width = 10, height = 3)

## Changes in enhancers
set = abc[abc$ENH_TC.sig & abc$ENH_TC.log2FoldChange > 0,]
gainedE = setBoxplotsTC(set = set, description = "Gained enhancers")

## Changes in contact
set = abc[abc$HIC_TCdiff & abc$HIC_TClfc > 0,]
gainedC = setBoxplotsTC(set = set, description = "Gained contact")

## Staticly looped pairs
ov = findOverlaps(query = abc2gi(abc[abc$ENH_TC.sig & abc$ENH_TC.log2FoldChange > 0,]),
                  subject = loops)
set = unique(abc[abc$ENH_TC.sig & abc$ENH_TC.log2FoldChange > 0,][queryHits(ov),])
gainedLoop = setBoxplotsTC(set = set, description = "Gained looped enh")

## For matching
set = unique(abc[abc$ENH_TC.sig & abc$ENH_TC.log2FoldChange > 0,])
ov = findOverlaps(query = abc2gi(set),
                  subject = loops)
looped = unique(set[queryHits(ov),])
unloop = unique(set[-queryHits(ov),])

## Contact matched unlooped pairs
set.seed(123)
contMatch = matchRanges(focal = looped,
                        pool = unloop,
                        covar = ~HIC_avgCount) 
set = as.data.frame(contMatch)
contMatchEnh = setBoxplotsTC(set = set, description = "Unlooped Gained Enh")

## Distance matched unlooped pairs
set.seed(123)
distMatch = matchRanges(focal = looped,
                        pool = unloop,
                        covar = ~HIC_span) 
set = as.data.frame(distMatch)
distMatchEnh = setBoxplotsTC(set = set, description = "Unlooped Gained Enh")

dev.off()


