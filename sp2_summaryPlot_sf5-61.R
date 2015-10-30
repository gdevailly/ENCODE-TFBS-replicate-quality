###
# This generate one summary multipanel plot for each figure
###

setwd("myDirectory/")
library(GenomicRanges)
library(Repitools)
load("Rdata/binaryRepMatrix_57.RData")
load("Rdata/listOfFPKMtable.repplicates_57.RData")
load("Rdata/listOfAllRep_57.RData")
load("Rdata/TSSproximity.replicates_57.RData")
load("Rdata/listOfCovMat.replicates_57.RData")
categ<-read.table("classificationOfEncodeReplicates_57.txt")
colnames(categ) <- c("name", "category")
rownames(categ) <- categ$name
levels(categ$category) <- c("dissimilar", "sensitive", "similar")

TF <- names(listOfAllRep)
TF <- strsplit(TF, "_")
TF <- sapply(TF, function(x) tail(x,1))
figT <- cbind("name"=names(listOfAllRep)[order(TF)], "index"= 5:(length(names(listOfAllRep)) + 4)) # this is for proper numerotation of the figure

doMegaPlot <- function(idN){
    # data preparation
    name <- figT[idN - 4, "name"]
    metaD <- listOfAllRep[[name]]
    
    # peak matrix
    # we summarised peak matrix into a 200 x N repplicate matrix, smaller and faster to plot
    binMat <- expMatrix[[name]][,5:ncol(expMatrix[[name]])]
    Nrep <- ncol(binMat)
    for (i in 1:ncol(binMat)){
        binMat <- binMat[order(binMat[,i], decreasing = TRUE),]
    }
    binMatRed <- matrix(nrow = 200, ncol = ncol(binMat))
    breaks <- round(seq(from = 1, to = nrow(binMat), length.out = 201))
    for(i in 1:200) {
        binMatRed[i,] <- colMeans(binMat[breaks[i]:breaks[i + 1],])
    }
    
    cov <- listOfCovMat[[name]]
    fpkmMat <- listOfFPKMtable[[name]]
    bpMat <- matrix(nrow = 2, ncol = Nrep, dimnames = list("row" = c("common", "other"), "col" = metaD$name))
    bpMat["common",] <- length(which(rowSums(binMat[,1:Nrep]) == Nrep))
    for(i in 1:Nrep){
        bpMat["other", i] <- length(which(binMat[,i] == 1)) - bpMat["common", i]
    }
    bxpMat <- fpkmMat$FPKMtable[which(rowSums(binMat) == Nrep), 5:(Nrep+4)]
    distTSS <- listOfPeakPostions[[name]]
    
    # plot parameters
    lwd <- 1.4
    lim <- ceiling(max(mapply(summary, bxpMat)[5,])*2)
    limRN <- ceiling(max(fpkmMat$nReads)/1E6)*1E6
    spacing <- 25
    mar <- c(3, 4, 2.2, 1)
    titleline <- 0.25
    letterCex <- 1.35
    lAdj <- -0.4
    cexText <- 0.96
    layoutMat <- matrix(c(1,3,8,13,1,4,9,13,2,5,10,13,2,6,11,13,2,7,12,13), nrow = 4, ncol = 5)
    layoutDim <- list("row" = c(1.6,1,1.3,0.2), "col" =c (1,1,1,1,1))
    height <- 5.4
    space <- 0.2
    if(Nrep == 2) space = 0.6
    
    #actual ploting
    pdf(file = paste0("Rgraphs/replicates/megaplot2/Supplementary figure ", idN, " ", gsub("_", " ", name), " (", categ[name,"category"], ").pdf"), width = 6.7, height = height)
    layout(layoutMat, widths = layoutDim$col, heights = layoutDim$row)
    
    # peak overlap heatmap
    par(mar = mar + c(0, -2, 0, 1), las = 1, mgp = c(2,1,0), xpd = TRUE, cex.axis = cexText, cex.lab = cexText, cex.main = cexText, cex.sub = cexText)
    plot.new()
    legend(x = "top", legend = paste0(1:ncol(bpMat), " - ", sub("wgEncodeAwgTfbs", "", colnames(bpMat))), inset = 0.3, cex = 1)
    mtext(paste0(gsub("_", " ", name), " \n(", categ[name,"category"],")"), line = -2.8, font = 2, cex = 1.1)
    ml <- legend(x = "bottomright", legend = c("", ""), fill = c("black", "white"), bty = "n", inset = c(-0.19, -0.24), trace = TRUE)
    text(ml$text$x-0.07, ml$text$y-0.012, c("Peak detected", "No peak detected"), pos=2)
    mtext("A", adj = -0.06, padj = 0, cex = letterCex)            
    
    image(as.matrix(binMatRed)[, Nrep:1], col = c("white", "black"), useRaster = FALSE, axes = FALSE)
    box()
    title(main = "Peak overlap", line = titleline, cex.main = 1.2)
    title(xlab = "Number of peaks", cex.lab = 1.2)
    axis(1, at = c(0, length(which(rowSums(binMatRed) == Nrep))/199, 1),
         labels = c(0, length(which(rowSums(binMat) == Nrep )), nrow(binMat)), cex.axis = 1.2)
    axis(2, at = seq(0,1, length.out = Nrep), labels = Nrep:1, font = 2,  cex.axis = 1.2)
    mtext("B", adj = -0.1, padj = 0, cex = letterCex) 
    
    # mean peak coverage
    lAdj <- -0.8
    par(mar = mar, las = 1, mgp = c(2,1,0), xpd = T)
    for(i in 1:Nrep) {
        yMax <- ceiling(max(max(cov[[i]][,1]), max(cov[[i]][,2])))
        plot(cov[[i]][,1] ~ rownames(cov[[i]]), type = "l", lwd = lwd, axes = FALSE,
             xlab = "Distance to peak center (bp)", ylab = "FPKM", ylim = c(0, yMax))
        title(main = paste("Mean\n coverage for", i), line = titleline)
        lines(cov[[i]][,2] ~ rownames(cov[[i]]), lty = 5, lwd = lwd)
        lines(cov[[i]][,3] ~ rownames(cov[[i]]), col = "darkgray", lwd = lwd)
        axis(1, at = c(-2000, 0, 2000))
        axis(2, at = c(0, yMax), labels=c(0, floor(yMax*10/3)))
        if (i == 1) mtext("C", adj = lAdj, padj = 0, cex = letterCex) 
    }
    if (Nrep != 5) {
        par(mar = c(0,0,0,0), xpd = NA)
        plot.new()
        legend(x = "left", legend = c("Common peaks", "Sample specific peaks", "Undetected peaks"), lty = c(1,5,1),
               title = "Mean coverage at:",
               col = c("black", "black", "darkgray"), bty = "n", lwd = lwd, xpd = NA)
        par(mar = mar, xpd = TRUE)
    }
    if(Nrep == 3){
        plot.new()
    } else if (Nrep == 2){
        plot.new()
        plot.new()
    } 
    
    # Number of peaks
    par(mar = mar + c(1.5,0,0.5,0), las = 1, xpd = TRUE)
    mp1 <- barplot(bpMat, col = c("black", "grey"), axes = FALSE, axisnames = FALSE, 
                 ylim = c(-max(colSums(bpMat))/spacing, max(colSums(bpMat))), space = space)
    title(main = "Number of peaks", line = titleline)
    axis(1, at = mp1, labels = 1:Nrep, font = 2)
    axis(2, at = c(0, bpMat[1,1], max(colSums(bpMat))))
    legend(x = "bottom", legend = c("Common", "Not common"), fill = c("black", "grey"),  bty = "n", inset = c(0, -0.84))
    mtext("D", adj = lAdj, padj = lAdj/2, cex = letterCex) 
    
    # Number of reads
    par(las = 1)
    mp2 <- barplot(fpkmMat$nReads, col = c("lightgrey"), axes = FALSE, axisnames = FALSE, ylim = c(-limRN/spacing, limRN), space = space)
    title(main = "Sequencing\n depth", line = titleline)
    title(ylab = "Millions of reads")
    axis(1, at = mp2, labels = 1:Nrep, font = 2)
    axis(2, at = c(0, limRN), labels = c(0, limRN/1E6))
    mtext("E", adj = lAdj, padj = lAdj/3, cex = letterCex) 
    
    # FPKM by peaks
    par(las = 1, xpd = FALSE, mar = mar + c(1.5,0,0.5,-0.6))
    boxplot(bxpMat,  pch = 20, ylim = c(-lim/spacing, lim), axes = FALSE, ylab = "FPKM", boxwex = 0.8) 
    title(main="Common\n peak height", line = titleline)
    axis(1, at = 1:Nrep, font = 2)
    axis(2, at = c(0, lim))
    mtext("F", adj = lAdj, padj = lAdj/3, cex = letterCex) 
    
    # distance bias
    par(xpd = T, mar = mar + c(1.5,0,0.5,0))
    col <- colorRampPalette(c("black", "white", "black"))(n = 7)
    bpd <- barplot(distTSS[nrow(distTSS):1,ncol(distTSS):1], col = col, axes = FALSE, axisnames = FALSE, horiz = TRUE)
    title(main = "Peak distribution", line = titleline)
    title(xlab = "Fraction of peaks")
    axis(2, at = bpd, labels = c(Nrep:1, "Common"), font = 2)
    axis(1, at = c(0, 0.5, 1))
    mtext("G", adj = lAdj, padj = lAdj/3, cex = letterCex) 
    par(xpd = FALSE)
    par(xpd = TRUE)
    plot.new()
    legend(x = -0.8, y = 1.2, legend = c("< -10 kb", "> -10 kb", "> -1 kb", "at TSS", "< +1 kb", "< +10 kb", "> +10 kb"),
           title = "Closest TSS:", fill = col, bty = "n", ncol = 1)
    #subtitle
    par(mar = c(0,0.2,0,0), xpd = F)
    plot.new()
    mtext(paste0("Supplementary figure ", idN, ": ", gsub("_", " ", name), " (", categ[name,"category"], ")"), adj = 0.02, padj = 2)
    dev.off()
}

a <- Sys.time()
empty <- lapply(5:61, doMegaPlot)
Sys.time() - a # 14s

