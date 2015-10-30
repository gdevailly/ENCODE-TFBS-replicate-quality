# Figure for irreproducible motives

library(motifStack)
setwd("myDirectory/")
load("Rdata/listOfAllRep_57.RData")
categ<-read.table("classificationOfEncodeReplicates_57.txt")
colnames(categ) <- c("name", "category")
rownames(categ) <- categ$name
levels(categ$category) <- c("dissimilar", "sensitive", "similar")

discrepencies <- c("K562_None_ATF3", "H1-hESC_None_CHD1", "H1-hESC_None_c-Myc",  "K562_None_HDAC2",
                 "K562_None_p300", "GM12878_None_PAX5", "K562_None_Pol2", "K562_None_Pol2_phosphoS2")
names(discrepencies) <- discrepencies

proportion <- function(x){
    rs <- sum(x);
    return(x / rs);
}

fomatMotifData <- function(repName){
    tempT <- listOfAllRep[[repName]]
    homMat <- list("replicate 1" = read.table(paste0("homerMotifOutput/individualTF/", tempT[1,"name"], ".narrowPeak/homerResults/motif1.motif"),
                                skip = 1),
                 "replicate 2" = read.table(paste0("homerMotifOutput/individualTF/", tempT[2,"name"], ".narrowPeak/homerResults/motif1.motif"),
                                skip = 1),
                 "common peaks" = read.table(paste0("homerMotifOutput/repCommonPeaks/", repName, "/homerResults/motif1.motif"),
                                skip = 1)
    )
    if (repName ==  "K562_None_Pol2") names(homMat) <- c("replicates 1, 3, 4", "replicate 2", "common peaks")
    homMat <- lapply(homMat, function(x) apply(x, 1, proportion))
    rownames(homMat[[1]]) <- c("A", "C", "G", "T")
    rownames(homMat[[2]]) <- c("A", "C", "G", "T")
    rownames(homMat[[3]]) <- c("A", "C", "G", "T")
    
    homMot <- lapply(names(homMat), function(x) new("pcm", mat = as.matrix(homMat[[x]]), name = x, color = colorset("DNA")))
    names(homMot) <- names(homMat)
    return(homMot)
}

listOfMot <- lapply(discrepencies, fomatMotifData)

listOfMot$K562_None_p300[["common peaks"]] <- matrixReverseComplement(listOfMot$K562_None_p300[["common peaks"]])

layoutMat <- matrix(c(1:27), nrow = 9, ncol = 3, byrow = T)
cexText <- 0.96
letterCex <- 1.35

pdf(file = "Rgraphs/replicates/motives_u.pdf", width = 5, height = 6.68)
layout(layoutMat)
par(mar = c(0.6,2,3,0.3), cex.axis = cexText, cex.lab = cexText, cex.main = cexText, cex.sub = cexText)
for(i in 1:length(listOfMot)) {
    plot(listOfMot[[i]][[1]],  xlcex = 0.6, ylcex = 0.6, ncex = 0.8, xaxis = FALSE, yaxis = FALSE, xlab = NA, ylab = NA)
    axis(side = 2, at = c(0, 0.5, 1), labels = c(0,1,2))
    mtext(LETTERS[3*i - 2], adj = -0.2, padj = -0.15, cex = letterCex) 
    mtext(paste0(gsub("_", " ", names(listOfMot)[i]), " (", categ[names(listOfMot)[i],"category"],"):"),
          adj = -0, padj = -2.1, cex = 0.9)
    
    plot(listOfMot[[i]][[2]],  xlcex = 0.6, ylcex = 0.6, ncex = 0.8, xaxis = FALSE, yaxis = FALSE, xlab = NA, ylab = NA)
    axis(side = 2, at = c(0, 0.5, 1), labels = c(0,1,2))
    mtext(LETTERS[3*i - 1], adj = -0.2, padj = -0.15, cex = letterCex) 
    
    plot(listOfMot[[i]][[3]],  xlcex = 0.6, ylcex = 0.6, ncex = 0.8, xaxis = FALSE, yaxis = FALSE, xlab = NA, ylab = NA)
    axis(side = 2, at = c(0, 0.5, 1), labels = c(0,1,2))
    mtext(LETTERS[3*i], adj = -0.2, padj = -0.15, cex = letterCex) 
}
par(mar = c(0,0.1,0,0))
plot.new()
mtext("Supplementary figure 3: Differences in de novo motif discovery from different\nreplicate experiments of 8 experimental conditions.",
            side = 3, adj = 0, padj = 1.2, cex = 0.8)
dev.off()

system("firefox Rgraphs/replicates/motives_u.pdf &")


