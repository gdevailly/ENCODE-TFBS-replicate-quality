setwd("myDirectory/")
library(GenomicRanges)
library(Repitools)
load("Rdata/binaryRepMatrix_57.RData")
load("Rdata/listOfFPKMtable.repplicates_57.RData")
load("Rdata/listOfAllRep_57.RData")
load("Rdata/TSSproximity.replicates_57.RData")
load("Rdata/listOfCovMat.replicates_57.RData")
categ <- read.table("classificationOfEncodeReplicates_57.txt")
colnames(categ) <- c("name", "category")
rownames(categ) <- categ$name

library(parallel)
library(lsr)

doNiceAAD <- function(name){
    metaD <- listOfAllRep[[name]]
    binMat <- expMatrix[[name]][, 5:ncol(expMatrix[[name]])]
    Nrep <- ncol(binMat)
    pFPKM <- listOfFPKMtable[[name]]
    
    nPeaks <- colSums(binMat)
    nReads <- pFPKM$nReads
    medFPKM <- apply(pFPKM$FPKMtable[,5:(4 + Nrep)], 2, median)
    
    values <- c(
        "P" = aad(nPeaks)*2,
        "R" = aad(nReads)*2*ifelse(lm(nReads~nPeaks)$coefficients[2] >= 0, 1, -1),
        "F" = aad(medFPKM)*2*ifelse(lm(medFPKM~nPeaks)$coefficients[2] >= 0, 1, -1)
    )
    return(values)
}

a <- Sys.time()
values <- mclapply(names(listOfAllRep), doNiceAAD, mc.preschedule = FALSE, mc.cores = 16) # adpat N cores to you computer
Sys.time() - a 
values <- do.call(rbind, values)
rownames(values) <- names(listOfAllRep)
colnames(values) <- c("P", "R", "F")

size<-which(categ$category=="sensitivity")
good<-which(categ$category=="similar")
bad<-which(categ$category=="dissimilar")

pdf(file = "Rgraphs/replicates/NEWcorrlations.pdf", width = 5, height = 3, pointsize = 8)
layout(cbind(1,2))
par(mar = c(6,3,2,1),  mgp = c(2,1,0))
plot(values[,c("R", "P")], pch = 20, col = "gray75", xlab = expression(paste(Delta,"Reads")), ylab = expression(paste(Delta,"Peaks")), xlim = c(-5E7, 12E7))
points(values[size,c("R", "P")], pch = 20)
text(values[size,c("R", "P")], rownames(values)[size], pos = 4, col = "gray50", cex = 0.8)
abline(h = 0, lwd = 0.6)
abline(v = 0, lwd = 0.6)
abline(lm(values[, "P"]~values[,"R"]), lty = 2)
legend("topleft", legend = c("\"Sensitive\"", "Other"), col = c("black", "gray50"), pch = 20, bty = "o", bg = rgb(255,255,255,80, maxColorValue = 255))
mtext("A", at = -5E7, cex = 2)
mtext("Supplementary figure 2: Influence of sequencing depth and common\npeaks height on the number of peaks detected",
      side = 1, adj = 0, padj = 2.4, cex = 1.2)

plot(values[,c("F", "P")], pch = 20, col="gray75", xlab = expression(paste(Delta,"Median FPKM")), ylab = expression(paste(Delta,"Peaks")), xlim = c(-5,24))
points(values[size,c("F", "P")], pch = 20)
text(values[size,c("F", "P")], rownames(values)[size], pos = 4, col = "gray50", cex = 0.8)
abline(h = 0, lwd = 0.6)
abline(v = 0, lwd = 0.6)
abline(lm(values[, "P"]~values[,"F"]), lty = 2)
legend("topleft", legend = c("\"Sensitive\"", "Other"), col = c("black", "gray50"), pch = 20, bty = "o",  bg = rgb(255,255,255,80, maxColorValue = 255))
mtext("B", at = -5, cex = 2)
dev.off()

system("firefox Rgraphs/replicates/NEWcorrlations.pdf &")

