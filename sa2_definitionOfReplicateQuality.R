# classification of replicates quality based on % overlap
setwd("myDirectory/")
library(GenomicRanges)
files <- list.files("TFtables/mergedDuplicates/")
files <- files[grep(".merged.bed$", files)]
merged.data <- lapply(files, function(x) read.table(paste0("TFtables/mergedDuplicates/",x)))
for (i in 1:length(merged.data)) colnames(merged.data[[i]]) <- c("chr", "start", "end", "N")
load("Rdata/listOfAllRep_57.RData")
names(merged.data) <- names(listOfAllRep)
merged.data.gr <- lapply(merged.data, function(x) GRanges(x$chr, IRanges(x$start, x$end)))

doExpMatrix <- function(name){
    table <- listOfAllRep[[name]]
    TFfiles <- table$name
    TFfiles <- paste0(TFfiles, ".narrowPeak.s")
    TFdata <- lapply(TFfiles, function(x) read.table(paste0("TFtables/sortedTFtables/" ,x)))
    names(TFdata) <- table$name
    TFdata.gr <- lapply(TFdata, function(x) GRanges(x[,1], IRanges(x[,2], x[,3])))
    hits <- lapply(TFdata.gr, function(x) findOverlaps(x, merged.data.gr[[name]], select="first"))
    expMatrix <- matrix(data=0, nrow=nrow(merged.data[[name]]), ncol=nrow(table))
    colnames(expMatrix) <- table$name
    for(i in 1:nrow(table)){
        expMatrix[hits[[i]],i] <- 1
    }
    return(cbind(merged.data[[name]], expMatrix))
}

library(parallel)
# change "mc.cores" depending on how many of them you want to use
a <- Sys.time()
expMatrix <- mclapply(names(listOfAllRep), doExpMatrix, mc.preschedule=F, mc.cores =32)
Sys.time() - a # 12s <3 <3 <3
names(expMatrix) <- names(listOfAllRep)

thresholds <- c("good" = 0.7, "bad" = 0.5, "size" = 0.5)

classifyRepplicates <- function(expName){
    mat <- expMatrix[[expName]]
    mat <- mat[,5:ncol(mat)]
    nPeaks <- colSums(mat)
    ratios <- sapply(colnames(mat), function(x) {
        length(which(rowSums(mat) == ncol(mat)))/nPeaks[x] 
    })
    names(ratios) <- colnames(mat)
    category <- "similar"
    if (any(ratios < thresholds["bad"])) {
        category <- "dissimilar"
    }
    if (min(nPeaks)/max(nPeaks) < thresholds["size"]) {
        nMin<-names(nPeaks[which(nPeaks == min(nPeaks))])
        if (ratios[nMin] > thresholds["good"]) {
            category <- "sensitivity"
        }
    }
    return(category)
}

categ <- sapply(names(listOfAllRep), classifyRepplicates)
table(categ)
categ <- cbind(names(categ), categ)
write.table(categ, file = "classificationOfEncodeReplicates_57.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
