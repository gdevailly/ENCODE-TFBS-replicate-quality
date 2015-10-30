###
# This script generate Figure 1
# Due to poor par() initialistation, it needs to run twice in order to produce a nice figure...
###

setwd("myDirectory/")
library(GenomicRanges)
library(Repitools)
library(motifStack)
load("Rdata/binaryRepMatrix_57.RData")
load("Rdata/listOfFPKMtable.repplicates_57.RData")
load("Rdata/listOfAllRep_57.RData")
load("Rdata/TSSproximity.replicates_57.RData")
load("Rdata/listOfCovMat.replicates_57.RData")
categ<-read.table("classificationOfEncodeReplicates_57.txt")
colnames(categ) <- c("name", "category")
rownames(categ) <- categ$name

# graph parameters
layoutMat <- matrix(c(1,1,4,5,1,1,4,6,21,21,21,21,3,3,20,9,3,3,20,10,16:19,2,2,11,7,2,2,11,8,12:15),
                  nrow = 9, ncol = 4, byrow = TRUE)

layoutDim <- list("row" = c(1,1,0.2,1,1,2,1,1,2), "col" = c(1.2,0.8,1,1))
cexText <- 0.96
letterCex <- 1.35
titleline <- 0.3
lwd <- 1.4
lAdj <- -0.5
spacing <- 25
height <- 5
space <- 0.6

levels(categ$category) <- c("dissimilar", "sensitive", "similar")

#ploting 
pdf(file = "Rgraphs/replicates/new_figure_1_v.pdf", width = 6.69, height = 7.5)

layout(layoutMat, widths=layoutDim$col, heights=layoutDim$row)
par(mar = c(3,3,2.5,1.5), las = 1, mgp = c(2,1,0), xpd = TRUE, cex.axis = cexText, cex.lab = cexText, cex.main = cexText, cex.sub = cexText)
f1_names <- c("HepG2_None_Rad21", "K562_None_HDAC2", "HepG2_None_BHLHE40")

# peak heatmaps --------------
letter <- "A"
for (name in f1_names) {
    mt <- paste0(gsub("_", " ", name), " \n(", categ[name,"category"],")")
    binMat <- expMatrix[[name]][, 5:ncol(expMatrix[[name]])]
    Nrep <- ncol(binMat)
    for (i in 1:ncol(binMat)){
        binMat <- binMat[order(binMat[,i], decreasing = TRUE),]
    }
    binMatRed <- matrix(nrow = 200, ncol = ncol(binMat))
    breaks <- round(seq(from = 1, to = nrow(binMat), length.out = 201))
    for(i in 1:200) {
        binMatRed[i,]<-colMeans(binMat[breaks[i]:breaks[i+1],])
    }
    
    image(as.matrix(binMatRed)[,Nrep:1], col = c("white", "black"), useRaster = FALSE, axes = FALSE)
    box()
    title(main = mt, line = titleline )
    title(xlab = "Number of peaks")
    axis(1, at = c(0, length(which(rowSums(binMatRed) == Nrep))/199, 1),
         labels = c(0, length(which(rowSums(binMat) == Nrep )), nrow(binMat)))
    axis(2, at = seq(0,1, length.out=Nrep), labels = Nrep:1, font = 2, cex.axis = 1.2)
    mtext(letter, adj = -0.1, padj = -0.2, cex = letterCex) 
    if (letter == "A") letter = "C"
    else if (letter == "C") letter = "B"
}
par(mar = c(0,0,0,0), xpd = NA)
plot.new()
legend(x = "topleft", legend = c("Peak detected", "No peak detected"), fill = c("black", "white"), bty = "n", inset = c(-0.12, 0.16))
ml <- legend(x = "bottomright", legend = c("","",""), lty = c(1,5,1),
           title = "Mean coverage at:", title.adj = -10, inset = c(-0.2, 0),
           col = c("black", "black", "darkgray"), bty = "n", lwd = lwd, xpd = T)
text(ml$text$x - 0.17, ml$text$y - 0.005, c("Common peaks", "Sample specific peaks", "Undetected peaks"), pos = 2)

# coverages ----------------------
par(mar = c(1.8,3,1.2,1.5), las = 1, mgp = c(2,1,0), xpd = TRUE, cex.axis = cexText, cex.lab = cexText, cex.main = cexText, cex.sub = cexText)

for (name in f1_names) {
    mt <- paste0(gsub("_", " ", name), " (", categ[name,"category"],")")
    cov <- listOfCovMat[[name]]
    binMat <- expMatrix[[name]][,5:ncol(expMatrix[[name]])]
    Nrep <- ncol(binMat)
    
    for(i in 1:Nrep) {
        yMax <- ceiling(max(max(cov[[i]][,1]), max(cov[[i]][,2])))
        plot(cov[[i]][,1] ~ rownames(cov[[i]]), type = "l", lwd = lwd, axes = FALSE,
             xlab = "Distance to peak center (bp)", ylab = "FPKM", ylim = c(0, yMax))
        title(main = paste("Mean coverage for", i), line = titleline)
        lines(cov[[i]][,2] ~ rownames(cov[[i]]), lty = 5, lwd = lwd)
        lines(cov[[i]][,3] ~ rownames(cov[[i]]), col = "darkgray", lwd = lwd)
        axis(1, at = c(-2000, 0, 2000))
        axis(2, at = c(0, yMax), labels = c(0, floor(yMax*10/3)))
    }
}

par(mar = c(0,0,0,0), xpd = NA)
plot.new()
legend(x = "topleft", legend = c("Peak detected", "No peak detected"), fill = c("black", "white"), bty = "n", inset = c(-0.12, 0.16))
ml <- legend(x = "bottomright", legend = c("","",""), lty = c(1,5,1),
           title = "Mean coverage at:", title.adj = -10, inset = c(-0.2, 0),
           col = c("black", "black", "darkgray"), bty = "n", lwd = lwd, xpd = TRUE)
text(ml$text$x - 0.17, ml$text$y - 0.005, c("Common peaks", "Sample specific peaks", "Undetected peaks"), pos = 2)

# bad analysis ----------------------
name <- f1_names[2]
metaD <- listOfAllRep[[name]]
distTSS <- listOfPeakPostions[[name]]
binMat <- expMatrix[[name]][,5:ncol(expMatrix[[name]])]
Nrep <- ncol(binMat)
mt <- paste0(gsub("_", " ", name), " \n(", categ[name,"category"],")")
pref <- "homerMotifOutput/individualTF/"
suf <- ".narrowPeak/homerResults/motif1.motif"

homMat <- list("1" = read.table(paste0(pref, metaD[1,"name"], suf), skip = 1),
             "2" = read.table(paste0(pref, metaD[2,"name"], suf), skip = 1) 
)

par(mar = c(4,5,2.5,1.5), las = 1, mgp = c(2,1,0), xpd = FALSE, cex.axis = cexText, cex.lab = cexText, cex.main = cexText, cex.sub = cexText)

col <- colorRampPalette(c("black", "white", "black"))(n = 7)
bpd <- barplot(distTSS[nrow(distTSS):1,ncol(distTSS):1], col = col, axes = FALSE, axisnames = FALSE, horiz = TRUE)
title(main = "Peak distribution", line = titleline)
title(xlab = "Fraction of peaks")
axis(2, at = bpd, labels = c(Nrep:1, "Common"), font = 2)
axis(1, at = c(0, 0.5, 1))
par(mar = c(0,0,0,0))
plot.new()
legend(x = "left", legend = c("< -10 kb", "> -10 kb", "> -1 kb", "at TSS", "< +1 kb", "< +10 kb", "> +10 kb"),
       title = "Closest TSS:", fill = col, bty = "n", ncol = 1)

proportion <- function(x){
    rs <- sum(x);
    return(x / rs);
}

par(mar = c(4,2,3,1), xpd = FALSE, las = 0)
homMat[[1]] <- apply(homMat[[1]], 1, proportion)
rownames(homMat[[1]]) <- c("A", "C", "G", "T")
motif1 <- new("pcm", mat = as.matrix(homMat[[1]]), name = "Top motif for 1")
plot(motif1, xlcex = 0.6, ylcex = 0.6, ncex = 0.7)

homMat[[2]] <- apply(homMat[[2]], 1, proportion)
rownames(homMat[[2]]) <- c("A", "C", "G", "T")
motif2 <- new("pcm", mat = as.matrix(homMat[[2]]), name = "Top motif for 2")
plot(motif2, xlcex = 0.6, ylcex = 0.6, ncex = 0.7)

# size analysis --------------------------
name <- f1_names[3]
mt <- paste0(gsub("_", " ", name), " \n(", categ[name,"category"],")")

metaD <- listOfAllRep[[name]]
binMat <- expMatrix[[name]][,5:ncol(expMatrix[[name]])]
Nrep <- ncol(binMat)
mt <- paste0(gsub("_", " ", name), " \n(", categ[name,"category"],")")

fpkmMat <- listOfFPKMtable[[name]]
bpMat <- matrix(nrow = 2, ncol = Nrep, dimnames = list("row" = c("common", "other"), "col" = metaD$name))
bpMat["common",] <- length(which(rowSums(binMat[,1:Nrep]) == Nrep))
lim <- ceiling(max(mapply(summary, bxpMat)[5,])*1.8)
limRN <- ceiling(max(fpkmMat$nReads)/1E6)*1E6
for(i in 1:Nrep){
    bpMat["other", i] <- length(which(binMat[,i] == 1)) - bpMat["common", i]
}
bxpMat <- fpkmMat$FPKMtable[which(rowSums(binMat) == Nrep), 5:(Nrep+4)]

# Number of peaks
par(mar = c(3,6,3,1), las = 1, xpd = FALSE)
mp1 <- barplot(bpMat, col=c("black", "grey"), axes = FALSE, axisnames = FALSE, 
             ylim = c(-max(colSums(bpMat))/spacing, max(colSums(bpMat))), space = space)
title(main = "Number of peaks", line = titleline)
axis(1, at = mp1, labels=1:Nrep, font = 2)
axis(2, at = c(0, bpMat[1,1], max(colSums(bpMat))))
par(mar = c(0,0,0,0))
plot.new()
legend(x = "left", legend = c("Common", "Not common"), fill = c("black", "grey"),  bty = "n", inset = c(0,-0.84))

# Number of reads
par(mar = c(3,4,3,1), las = 1)
mp2 <- barplot(fpkmMat$nReads, col = c("lightgrey"), axes = FALSE, axisnames = FALSE, ylim = c(-limRN/spacing, limRN), space = space)
title(main = "Sequencing depth", line = titleline)
title(ylab = "Millions of reads")
axis(1, at = mp2, labels = 1:Nrep, font = 2)
axis(2, at = c(0, limRN), labels = c(0, limRN/1E6))

# FPKM by peaks
par(mar = c(3,4,3,1), las = 1, xpd = FALSE)
boxplot(bxpMat,  pch = 20, ylim = c(-lim/spacing, lim), axes = FALSE, ylab = "FPKM", boxwex = 0.8, add = FALSE) 
title(main = "Common\n peak height", line = titleline)
axis(1, at = 1:Nrep, font = 2)
axis(2, at = c(0, lim))

par(mar = c(0,0,0,0), xpd = NA)
plot.new()
legend(x = "topleft", legend = c("Peak detected", "No peak detected"), fill = c("black", "white"), bty = "n", inset = c(-0.12, 0.16))
ml <- legend(x = "bottomright", legend=c("","",""), lty = c(1,5,1),
           title = "Mean coverage at:", title.adj = -10, inset = c(-0.2, 0),
           col = c("black", "black", "darkgray"), bty = "n", lwd = lwd, xpd = TRUE)
text(ml$text$x - 0.17, ml$text$y - 0.005, c("Common peaks", "Sample specific peaks", "Undetected peaks"), pos = 2)

# end of plot ------------
par(xpd = NA)
abline(h = grconvertX(-0.0022, 'ndc', 'device'), lty = 2)
abline(h = grconvertX(0.0022, 'ndc', 'device'), lty = 2)
dev.off()



