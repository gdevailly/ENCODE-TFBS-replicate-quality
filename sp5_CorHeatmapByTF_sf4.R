###
# This script generate correlation heatmaps for figure SX
# unfortunatly, it uses data made from another (yet unpublished) super project
# "Rdata/encodeTFmatrix.clostest.gencode.transcript.RData" contains precompute correlations
# of all ENCODE TFBS ChIP-seq (Uniform peaks) based on peak overlap
####

setwd("myDirectory/")
library(GenomicRanges)
library(Repitools)
load("Rdata/listOfAllRep_57.RData")

metaData<-read.table("TFtables/metadata_TFtables.txt", stringsAsFactors = FALSE)
metaData<-metaData[, 2:6]
colnames(metaData) <- c("lab", "cell", "treatment", "antibody", "name")

temp <- sapply(as.character(metaData$antibody), function(x) strsplit(x, "_"))
temp2 <- sapply(temp, function(x) x[[1]])
metaData <- cbind(metaData, "TF" = temp2)
rm(temp, temp2)

metaData$TF <- as.character(metaData$TF)
metaData[which(metaData$TF == "GATA-2"), "TF"] <- "GATA2"
metaData[which(metaData$TF == "GRp20"), "TF"] <- "GR"
metaData[which(metaData$TF == "PAX5-C20"), "TF"] <- "PAX5"
metaData[which(metaData$TF == "PAX5-N19"), "TF"] <- "PAX5"
metaData[which(metaData$TF == "Sin3Ak-20"), "TF"] <- "SIN3A"
metaData[which(metaData$TF == "USF-1"), "TF"] <- "USF1"
metaData[which(metaData$TF == "Pol2(b)"), "TF"] <- "Pol2_b"
metaData[which(metaData$TF == "Pol2(phosphoS2)"), "TF"] <- "Pol2_phosphoS2"

load("Rdata/encodeTFmatrix.clostest.gencode.transcript.RData") # missing from gitHub, mostly because it is an HEAVY RData object
# please contact me if needed, or wait a few more month until I finished this other project of mine

TFreg$filename <- sub("UniPk", "", sub("wgEncodeAwgTfbs" ,"" , metaData$name))
colnames(TFreg$mat) <- TFreg$filename

library(gplots)
library(RColorBrewer)
colors = c(seq(-0.5, 1,length = 100))
my_palette <- colorRampPalette(c("blue", "white", "red", "black"))(n = 99)
library(cba)
library(SDMTools)

indLetter = LETTERS[1:5]
names(indLetter) <- c("HDAC2", "p300", "ATF3", "NRSF",  "Pol2")

drawHeatmapForTFwithL <- function(TF) {
    if (TF == "HDAC2") {
        rTF <- which(metaData$TF %in% c("HDAC2", "HDAC6", "eGFP-HDAC8", "p300") &
                   metaData$cell %in% c("H1-hESC", "K562"))
    } else if (TF =="Pol2") {
        rTF <- which(metaData$TF == TF & metaData$cell == "K562")
    } else {
        rTF <- which(metaData$TF == TF)
    }
    
    mat <- TFreg$mat[,rTF]
    cmat <- cor(mat, method = "pearson")
    d <- as.dist(1 - cmat)
    hc <- hclust(d, method = "complete", members = NULL)        
    co <- order.optimal(d, hc$merge)
    hc$merge <- co$merge
    hc$order <- co$order
    cmat <- cmat[hc$order,hc$order]
    par(mai = c(1.2,0.25,0.25,1.2), las = 2, xpd = NA, cex.axis = 0.5, cex.main = 0.8)
    image(cmat, breaks = colors, col=my_palette, useRaster = FALSE, axes = FALSE)
    title(main = TF, cex.main = 1.2)
    axis(1, at = seq(0,1, length.out = ncol(cmat)), labels = colnames(cmat))        
    axis(4, at = seq(0,1, length.out = ncol(cmat)), labels = colnames(cmat))   
    mtext(indLetter[TF], adj = -0.2, padj = -0.15, cex = 1.35, las = 1) 
    legend.gradient(cbind(x = c(0.3,0.4,0.4,0.3) + 1, y = c(-0.8,-0.8,-1.6,-1.6)/1.2), my_palette, c("-0.5","1"),
                    title = "Pearson\ncorrelation", cex = 0.5)
}

intTF <- c("HDAC2", "p300", "ATF3", "NRSF",  "Pol2")

pdf(paste0("Rgraphs/corHeatmapsByTF/article_heatmapsTF_n1.pdf"), height = 7.2, width = 4.8)
layout(rbind(1:2,3:4,5:6))
empty <- lapply(intTF, drawHeatmapForTFwithL)
plot.new()
par(mar = c(0,0,0,0), las = 1)
mtext("Supplementary figure 4: Comparison of\nvarious peak lists with ChIP-seq of the\nsame factor done in other cell lines or\nconditions",
      side = 1, adj = 0, padj = -0.6, cex = 0.8, at = -0.6)
dev.off()

system("firefox Rgraphs/corHeatmapsByTF/article_heatmapsTF_n1.pdf &")
