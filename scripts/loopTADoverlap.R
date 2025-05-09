
# INTIALIZE -----------------------
library(GenomicRanges)
library(mariner)
library(eulerr)


# READ IN -----------------------

## Read in loops -------
loops = as_ginteractions(readRDS("./output/loopProcessing/loops.rds"))

## Read in TADs -------
tads = readRDS("./output/tadProcessing/tads.RDS")


# OVERLAP -----------------------
ov = findOverlaps(query = loops,
                  subject = tads, 
                  maxgap = 50000)
length(unique(queryHits(ov)))/length(loops) # 27.0% of loops overlap TADs
length(unique(subjectHits(ov)))/length(tads) # 52.4% of TADs overlap loops

# PLOT: Venn -----------------------
all = list("loops" = 1:length(loops),
           "TADs" = c(length(loops)+(1:length(tads))[-unique(subjectHits(ov))],
                      unique(queryHits(ov))))

pdf("./output/loopTADoverlap/FigS2G-loopTADvenn.pdf", width = 6, height = 6)
plot(euler(all), col = NA, fill = c("#41B6C4", "#225EA8"), quantities = T)
dev.off()


## PLOT: Loop/TAD sizes -----------------------
pdf("./output/loopTADoverlap/FigS2C-loopTADsizes.pdf", width = 6, height = 6)

plot(density(pairdist(loops)/1e6),
     col = "#41B6C4",
     main = "Feature sizes (Mb)",
     xlim = c(0, 4))
lines(density(pairdist(tads)/1e6),
      col = "#225EA8")
legend("topright", 
       legend = c(paste0("loops (mean = ", 
                         signif(mean(pairdist(loops)/1e3), 3), "kb)"), 
                  paste0("TADs (mean = ", 
                         signif(mean(pairdist(tads)/1e3), 3), "kb)")),
       text.col = c("#41B6C4", "#225EA8"),
       bty = 'n')

dev.off()

## PLOT: Loop/TAD size boxplots -----------------------
dat = list("other loops" = pairdist(loops[-unique(queryHits(ov))]),
           "TAD loops" = pairdist(loops[unique(queryHits(ov))]),
           "loop TADs" = pairdist(tads[unique(subjectHits(ov))]),
           "other TADs" = pairdist(tads[-unique(subjectHits(ov))]))

pdf("./output/loopTADoverlap/FigS2H-loopTADclassSizes.pdf",
    width = 6, height = 6)
par(mfrow=c(1,1))
boxplot(dat, outline = F,
        col = c("grey70", "grey30", "darkblue", "steelblue"))
dev.off()

