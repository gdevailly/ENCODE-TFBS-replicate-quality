setwd("myDirectory/")
library(GenomicRanges)
library(readr)
library(magrittr)
load("Rdata/listOfAllRep_57.RData")
load("Rdata/binaryRepMatrix_57.RData")
load("Rdata/binaryRepMatrix_peakSeq.RData")
load("Rdata/MAnorm_results.RData")
peakSeq <- read_tsv("peakSeq_table.txt")
peakSeq$peakSeq_fn <- sub(".bb", "", peakSeq$peakSeq_fn, fixed = TRUE)
categ<-read.table("classificationOfEncodeReplicates_57.txt")
colnames(categ) <- c("name", "category")
rownames(categ) <- categ$name

titleline <- 0.3
letterCex <- 1.35
cexText <- 0.96


peakComparisonPlot <- function(idN) {
    # data preparation
    myName <- unique(peakSeq$name)[idN]
    metaD <- listOfAllRep[[myName]]
    metaD_ps <- subset(peakSeq, peakSeq[,"name"] == myName)
    metaD_ps[,1:4]
    
    # peak matrix for Unform/SPP peaks
    binMat<-expMatrix[[myName]][,5:ncol(expMatrix[[myName]])]
    Nrep <- ncol(binMat)
    for (i in 1:ncol(binMat)){
        binMat <- binMat[order(binMat[,i], decreasing = TRUE),]
    }
    binMatRed <- matrix(nrow = 200, ncol = ncol(binMat))
    breaks <- round(seq(from = 1, to = nrow(binMat), length.out = 201))
    for(i in 1:200) {
        binMatRed[i,] <- colMeans(binMat[breaks[i]:breaks[i+1],])
    }
    
    #peak natrix for peakSeq peaks
    binMat_ps <- expMatrix_peakSeq[[myName]][,5:ncol(expMatrix_peakSeq[[myName]])]
    Nrep_ps <- ncol(binMat_ps)
    for (i in 1:ncol(binMat_ps)){
        binMat_ps <- binMat_ps[order(binMat_ps[,i], decreasing = TRUE),]
    }
    binMatRed_ps <- matrix(nrow = 200, ncol = ncol(binMat_ps))
    breaks <- round(seq(from = 1, to = nrow(binMat_ps), length.out = 201))
    for(i in 1:200) {
        binMatRed_ps[i,] <- colMeans(binMat_ps[breaks[i]:breaks[i+1],])
    }
    
    mt <- paste0(gsub("_", " ", myName), " \n")
    
    image(as.matrix(binMatRed)[,Nrep:1], col = c("white", "black"), useRaster = FALSE, axes = FALSE)
    box()
    title(main = paste0(mt, "(SPP)"), line = titleline )
    title(xlab = "Number of peaks")
    axis(1, at = c(0, length(which(rowSums(binMatRed) == Nrep))/199, 1),
         labels = c(0, length(which(rowSums(binMat) == Nrep )), nrow(binMat)))
    axis(2, at = seq(0,1, length.out = Nrep), labels = Nrep:1, font = 2, cex.axis = 1.2)
    mtext(LETTERS[3*idN - 2], adj = -0.1, padj = -0.2, cex = letterCex) 
    
    image(as.matrix(binMatRed_ps)[,Nrep_ps:1], col = c("white", "black"), useRaster = FALSE, axes = FALSE)
    box()
    title(main = paste0(mt, "(peakSeq)"), line = titleline )
    title(xlab = "Number of peaks")
    axis(1, at = c(0, length(which(rowSums(binMatRed_ps) == Nrep))/199, 1),
         labels = c(0, length(which(rowSums(binMat_ps) == Nrep_ps)), nrow(binMat_ps)))
    axis(2, at = seq(0,1, length.out=Nrep_ps), labels = Nrep_ps:1, font = 2, cex.axis = 1.2)
    mtext(LETTERS[3*idN - 1], adj = -0.1, padj = -0.2, cex = letterCex) 
    
    myPalette <- c("black", "gray" ,"lightgray" ,"white")
    image(as.matrix(MAnorm_results[[myName]]$red_sorted_pval), col = myPalette, axes = FALSE, useRaster = FALSE)
    box()
    axis(1, at = MAnorm_results[[myName]]$at, labels = MAnorm_results[[myName]]$labels)
    title(main = paste0(mt, "(MAnorm)"), line = titleline)
    title(xlab = "Number of peaks")
    mtext(LETTERS[3*idN], adj = -0.1, padj = -0.2, cex = letterCex) 

}


pdf(file="Rgraphs/replicates/sup_fig_1_spp_vs_peakseq_vs_MAnorm.pdf", width = 8, height = 8)
layout(matrix(1:24, ncol = 3, byrow = TRUE), heights = c(0.5, rep(1,6), 0.25))
par(mar = c(0,0,0,0), las = 1, mgp = c(2,1,0), xpd = TRUE, cex.axis = cexText, cex.lab = cexText, cex.main = cexText, cex.sub = cexText)
plot.new()
title(main = "ENCODE SPP peaks", cex.main = 1.2, line = -1.5)
legend(x = "bottom", legend = c("Peak detected", "No peak detected"), fill = c("black", "white"), bty = "n")
plot.new()
title(main = "ENCODE peakSeq peaks", cex.main = 1.2, line = -1.5)
legend(x = "bottom", legend = c("Peak detected", "No peak detected"), fill = c("black", "white"), bty = "n")
plot.new()
title(main = "SPP peak comparison\nwith MAnorm", cex.main = 1.2, line = -2.6, adj = 0)
legend(x = "right",
       legend = c(as.expression(bquote(p-value > 10^-2)), as.expression(bquote(p-value <= 10^-2)), as.expression(bquote(p-value <= 10^-3)), as.expression(bquote(p-value <= 10^-5))),
       fill = c("black", "gray" ,"lightgray" ,"white"), bty="n")
par(mar = c(3,3,2.5,1.5), las = 1, mgp = c(2,1,0), xpd = TRUE, cex.axis = cexText, cex.lab = cexText, cex.main = cexText, cex.sub = cexText)
empty <- lapply(1:6, peakComparisonPlot)
par(mar = c(0,1.5,0,0))
plot.new()
mtext("Supplementary figure 1: ENCODE SPP peaks, ENCODE peakSeq peaks, and MAnorm SPP peaks comparison.",
      side = 3, adj = 0, padj = 2, cex = 0.8)
dev.off()

system("firefox Rgraphs/replicates/sup_fig_1_spp_vs_peakseq_vs_MAnorm.pdf &")




