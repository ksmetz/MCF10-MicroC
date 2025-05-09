
# INITIALIZE -------------------------------------
library(tximeta)
library(DESeq2)
library(GenomicFeatures)
library(scales)
library(RColorBrewer)


# READ IN -------------------------------------

## Read in Salmon quantification outputs -------
samples = list.dirs("./input/rna/salmon_quant/", 
                    recursive = F, full.names = T)

## Parse names + relevant info --------
sampleNames = lapply(basename(samples), function(n){
  strsplit(n, "_")[[1]][[1]]
}) |> unlist()
repNames = lapply(basename(samples), function(n){
  strsplit(strsplit(n, "_")[[1]][[2]], "\\.")[[1]][[1]]
}) |> unlist()
uniqueNames = lapply(basename(samples), function(n){
  strsplit(n, "\\.")[[1]][[1]]
}) |> unlist()

## Find specific quant.sf files -------
files = list.files(samples, "quant.sf", full.names = T)
file.exists(files)



# PROCESS INPUTS -------------------------------------

## Importing to gene level -----------
## Make colData
coldata = data.frame(files,
                     names = uniqueNames,
                     sample = sampleNames,
                     rep = repNames)

## Because Salmon was run in alignment mode, there is no hash. Run as tximport
## Make TxDb from GENCODE GTF file
txdb = makeTxDbFromGFF(file = "~/Tools/reference/gencode.v36.annotation.gtf.gz")

## Generate tx2gene table
k = keys(txdb, keytype = "TXNAME")
tx2gene = select(txdb, k, "GENEID", "TXNAME")

## Create summarized experiment
se <- tximeta(coldata, skipMeta = T, txOut = F, tx2gene = tx2gene)


# DIFFERENTIAL GENES -------------------------------------
## Run DESeq
dds <- DESeqDataSet(se, ~sample)
dds <- DESeq(dds)

## Find results
resAC <- results(dds, name = "sample_MCF10CA1a_vs_MCF10A")
resAT <- results(dds, name = "sample_MCF10AT1_vs_MCF10A")

resAC <- lfcShrink(dds = dds, 
                   contrast = c("sample", "MCF10CA1a", "MCF10A"),
                   type = 'ashr')
resAT <- lfcShrink(dds = dds, 
                   contrast = c("sample", "MCF10AT1", "MCF10A"),
                   type = 'ashr')
resTC <- lfcShrink(dds = dds, 
                   contrast = c("sample", "MCF10CA1a", "MCF10AT1"),
                   type = 'ashr')


## Pull significantly differential genes
findSig <- function(res, pThresh, lThresh){
  dat = res[!is.na(res$padj),]
  dat = dat[dat$padj <= pThresh & abs(dat$log2FoldChange) >= lThresh,]
  return(dat)
}

p = 0.01
L = 0
sigAC <- findSig(resAC, p, L) 
sigAT <- findSig(resAT, p, L) 
sigTC <- findSig(resTC, p, L) 
allSig = unique(c(rownames(sigAC),
                  rownames(sigAT),
                  rownames(sigTC)))



# GENE INFO -------------------------------------
## Get gene symbols ------
ensemblIDs = keys(txdb, keytype = "GENEID")
trimEnsemblIDs = unlist(lapply(ensemblIDs, function(g){
  strsplit(g, "\\.")[[1]][[1]]
}))

## AnnotationHub: 60,419 out of 60,660 genes kept
library(AnnotationHub)
ah <- AnnotationHub()
ahDb = query(ah, c("EnsDb", "110", "sapiens"))[[1]]

keytypes(ahDb)
head(keys(ahDb, keytype = "GENEBIOTYPE"))

geneIDtable = select(ahDb, trimEnsemblIDs, columns = c("SYMBOL", "GENEBIOTYPE", "GENEID"), keytype="GENEID")

length(unique(geneIDtable$GENEID))

## Add in DESeq IDs for matching
geneIDtable$DESEQID = ensemblIDs[match(geneIDtable$GENEID, trimEnsemblIDs)] # Add IDs

## Adding gene coordinates -----------
## Get gene coordinates
geneInfo = GenomicFeatures::genes(txdb)
identical(geneInfo$gene_id, ensemblIDs)

## Match DESeq IDs to ensembl IDs
geneIDtable = geneIDtable[match(geneInfo$gene_id, geneIDtable$DESEQID),] # Reorder to full set
geneIDtable$DESEQID = ensemblIDs



## Make complete gene table ------------
## Prep differential info
identical(rownames(resAC), geneIDtable$DESEQID)
identical(rownames(resAT), geneIDtable$DESEQID)
identical(rownames(resTC), geneIDtable$DESEQID)

diffInfo = cbind(resAC, 
                 AC.sig = rownames(resAC) %in% rownames(sigAC),
                 resAT[,-1],
                 AT.sig = rownames(resAT) %in% rownames(sigAT),
                 resTC[,-1], 
                 TC.sig = rownames(resAC) %in% rownames(sigTC))
colnames(diffInfo)[2:5] = paste0("AC.", colnames(diffInfo)[2:5])
colnames(diffInfo)[7:10] = paste0("AT.", colnames(diffInfo)[7:10])
colnames(diffInfo)[12:15] = paste0("TC.", colnames(diffInfo)[12:15])

identical(rownames(diffInfo), geneIDtable$DESEQID)

## Prep coordinates
identical(geneInfo$gene_id, geneIDtable$DESEQID)

## Prep raw counts
identical(geneIDtable$DESEQID, rownames(counts(dds)))

## Add all to geneInfo object
mcols(geneInfo) = cbind(geneIDtable, diffInfo, counts(dds))


## PLOT: Number of DEGs ---------
DEGbyComp = matrix(c(sum(geneInfo$AT.log2FoldChange[geneInfo$AT.sig == T] > 0),
                     sum(geneInfo$AT.log2FoldChange[geneInfo$AT.sig == T] < 0),
                     sum(geneInfo$TC.log2FoldChange[geneInfo$TC.sig == T] > 0),
                     sum(geneInfo$TC.log2FoldChange[geneInfo$TC.sig == T] < 0),
                     sum(geneInfo$AC.log2FoldChange[geneInfo$AC.sig == T] > 0),
                     sum(geneInfo$AC.log2FoldChange[geneInfo$AC.sig == T] < 0)),
                   ncol = 3, byrow = F)
colnames(DEGbyComp) = c("A vs T", "T vs C", "A vs C")

pdf("./output/rnaProcessing/FigS3A-diffGeneNumbers.pdf",
    width = 6, height = 6)
barplot(DEGbyComp,
        col = c("#FD8D3C", "#FECC5C"),
        border = NA,
        main = "Number of differential genes")
legend("topleft",
       legend = c("Down-regulated", "Up-regulated"),
       text.col = c("#FECC5C", "#FD8D3C"),
       bty = 'n')
dev.off()



# CLUSTERING -------------------------------------
combineReps <- function(matrix, new.colnames)
{
  newMatrix = c()
  for (i in c(1, 4, 7))
  {
    tempMatrix = matrix[,i:(i+2)]
    newCol = apply(tempMatrix,1,mean)
    newMatrix  = cbind(newMatrix,newCol)
  }
  rownames(newMatrix) = rownames(matrix)
  
  if(!is.null(new.colnames)){
    colnames(newMatrix) = new.colnames
  }
  
  return(newMatrix)
}

# Build a matrix of transformed loop counts
countMat = assay(vsd)[,1:9]

# Center and scale data
countMat.norm <- (countMat - rowMeans(countMat))/rowSds(countMat + .5)

# Subset for significant loops
countMatSig.norm <- countMat.norm[which(rownames(countMat.norm) %in% allSig),]

# Combine replicates
countMat.combo <- combineReps(countMat, 
                              new.colnames=c("A", "T", "C"))
countMat.norm.combo <- combineReps(countMat.norm, 
                                   new.colnames=c("A", "T", "C"))
countMatSig.norm.combo <- combineReps(countMatSig.norm, 
                                      new.colnames=c("A", "T", "C"))

# Set seed to preserve manual ordering
seed = 32
set.seed(seed)

# Perform clustering
k = 4
clust = kmeans(countMatSig.norm, centers=k)
cut = clust$cluster

# Order clusters
clusterOrder = c(3, 1, 2, 4) 

names(clusterOrder) = c("up.early", "up.late", 
                        "down.early", "down.late")

## PLOT: Line Plots -----------
library(colorspace)
library(RColorBrewer)
library(scales)

k.pal = brewer.pal(n=4, name = "Spectral")


pdf("./output/rnaProcessing/FigS3B-clusterLineplots.pdf",
    width = 3, height = 12)

# Plot each cluster
par(mar=c(3,3,2,2))
par(mfrow=c(4, 1))

n=0
for (i in clusterOrder) {
  n=n+1
  
  # grab data
  dat = countMatSig.norm.combo[cut == i,]
  
  # make empty plot
  plot(x=1:3, 
       y=colMeans(countMatSig.norm.combo[cut==i,]), 
       type="n",
       main=paste(names(clusterOrder)[n],"\nn=",table(cut)[i],sep=""),
       ylim=c(-1,1), 
       xaxt="n", 
       xlab="Cancer progression", 
       ylab="Relative loop strength")
  
  # add fill (std dev)
  polygon(c(1:3, 3:1), 
          c(colMeans(dat)-colSds(dat), 
            rev(colMeans(dat)+colSds(dat))), 
          col=lighten(k.pal[n], amount=.85), border=NA)
  
  # and transparent grey lines
  #apply(countMatSig.norm.combo[cut==i,],1,lines,col=adjustcolor("grey",alpha.f=0.4), x=c(0, .25, .5, .75, 1,2,3,4))
  abline(h=0, col=alpha("black", .2))
  
  # add median line
  lines(x=1:3,
        y=colMeans(countMatSig.norm.combo[cut==i,]),
        col=k.pal[n],
        lwd=2)
  
  # add axis
  axis(1, at=1:3, 
       labels=c("A", "T", "C"))
}

dev.off()


## PLOT: Heatmap -----------
library(pheatmap)

png(filename="./output/rnaProcessing/FigS3B-clusterHeatmap.png",
    width=8, height=8, units="in", res=300)
par(mfrow=c(1,1))
par(mar=c(5,5,4,5))

pheatmap(countMatSig.norm.combo[order(match(cut, clusterOrder)),], 
         cluster_rows = F, show_rownames = F,
         cluster_cols = F,
         annotation_row = data.frame(cluster = match(cut, clusterOrder)[order(match(cut, clusterOrder))],
                                     row.names = rownames(countMatSig.norm.combo[order(match(cut, clusterOrder)),])),
         annotation_colors = list(cluster = k.pal),
         annotation_legend = F)

dev.off()


## TABLE: Gene info -----------
clusteredGenes = cbind(as.data.frame(geneInfo),
                       data.frame(
                         cluster = cut[match(x=rownames(countMat.norm.combo),
                                             table=names(cut))],
                         A_ZSCR = countMat.norm.combo[,1],
                         T_ZSCR = countMat.norm.combo[,2],
                         C_ZSCR = countMat.norm.combo[,3],
                         max_ZSCR = rowMaxs(countMat.norm.combo),
                         A_VST = countMat.combo[,1],
                         T_VST = countMat.combo[,2],
                         C_VST = countMat.combo[,3]))

clusteredGenes$cluster[which(is.na(clusteredGenes$cluster) == TRUE)] = "static"
for (i in 1:k){
  clusteredGenes$cluster[which(clusteredGenes$cluster == clusterOrder[i])] = names(clusterOrder[i])
}



# GO ENRICHMENT ----------------
library(gprofiler2) 

allGenes = clusteredGenes$GENEID
ATgenes = clusteredGenes$GENEID[clusteredGenes$AT.sig]
TCgenes = clusteredGenes$GENEID[clusteredGenes$TC.sig]
ACgenes = clusteredGenes$GENEID[clusteredGenes$AC.sig]

goAT = gost(ATgenes, organism="hsapiens", ordered_query=F, significant=T, 
            user_threshold=0.05, correction_method="fdr", sources=c("GO", "KEGG","REAC", "WP"),
            custom_bg=allGenes) # running GO
goTC = gost(TCgenes, organism="hsapiens", ordered_query=F, significant=T, 
            user_threshold=0.05, correction_method="fdr", sources=c("GO", "KEGG","REAC", "WP"),
            custom_bg=allGenes) # running GO
goAC = gost(ACgenes, organism="hsapiens", ordered_query=F, significant=T, 
            user_threshold=0.05, correction_method="fdr", sources=c("GO", "KEGG","REAC", "WP"),
            custom_bg=allGenes) # running GO

godat = function(goList, source, n = 30){ # plotting the result
  goList = lapply(goList, function(goterm){
    goterm = goterm$result
    goterm = goterm[goterm$term_size<800, ]
    goterm$geneRatio = goterm$intersection_size/goterm$term_size
    goterm = goterm[order(goterm$p_value),]
    goterm = goterm[goterm$source %in% source,]
    return(goterm)
  })
  
  topTerms = lapply(goList, function(goterm){
    goterm$term_name[1:n]
  })
  topTerms = unique(unlist(topTerms))
  
  topIDs = lapply(goList, function(goterm){
    goterm$term_id[1:n]
  })
  topIDs = unique(unlist(topIDs))
  
  goterm2plot = lapply(goList, function(goterm){
    goterm2plot = goterm[,c("term_name","p_value","intersection_size","term_id", "geneRatio")]
    goterm2plot = goterm2plot[match(topTerms, goterm2plot$term_name),]
    goterm2plot$log10P = -log10(goterm2plot$p_value)
    goterm2plot$term_name = topTerms
    goterm2plot$term_id = topIDs
    goterm2plot$term_name = factor(goterm2plot$term_name, levels=rev(goterm2plot$term_name))
    goterm2plot[is.na(goterm2plot)] = 0
    return(goterm2plot)
  })
  
  names(goterm2plot) = names(goList)
  return(goterm2plot)
}

goplot = function(dat, title, rFactor = 50){
  plot(x=0, y=0, type='n',
       yaxt = 'n', xaxt = 'n',
       ylab = NA, xlab = NA,
       bty = 'n',
       xlim = c(-16, length(dat)),
       ylim = c(1, nrow(dat[[1]])),
       xpd = T, main = title)
  for(i in 1:length(dat)){
    radius = sqrt(dat[[i]]$log10P/pi)
    symbols(x = rep(i, nrow(dat[[i]])),
            y = nrow(dat[[i]]):1,
            circles=radius/25,
            inches = F,
            bg = "black",
            fg = NA, add=T)
  }
  axis(side = 1, lwd = 0,
       at = 1:length(dat),
       labels = names(dat))
  text(x = rep(0.5, times = nrow(dat[[1]])),
       y = nrow(dat[[1]]):1,
       labels = dat[[1]]$term_name, adj = c(1, 0.5), 
       xpd = T)
}

datAll = godat(list("AvT" = goAT, 
                    "TvC" = goTC, 
                    "AvC" = goAC), source = "GO:BP", n = 10)

## PLOT: GO Term Enrichment ---------
pdf("./output/rnaProcessing/FigS3C-GOenrichment.pdf",
    width = 8, height = 6)
par(mfrow=c(1,1))
par(mar=c(3,3,4,3))
goplot(datAll, title = "GO term enrichment")
dev.off()

# GSEA INPUT -------------
dat = counts(dds, normalized = T)
dat = dat[,1:9]
identical(rownames(dat), geneInfo$DESEQID)

ids = rownames(dat)
ids = lapply(ids, function(i){strsplit(i, "\\.")[[1]][1]})
ids = unlist(ids)

dat = dat[!duplicated(ids),]

dat = data.frame(gene = ids[!duplicated(ids)],
                 description = NA,
                 dat)

write.table(dat, file = "./output/rnaProcessing/GSEA_input.txt", quote = F, 
            row.names = F, col.names = T, sep = "\t")


# OUTPUT ------------------
saveRDS(clusteredGenes, file = "./output/rnaProcessing/genes.rds")
