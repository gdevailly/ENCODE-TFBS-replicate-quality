# Create a binarymatrix peaks x replicates for each conditions
setwd("myDirectory/")
library(GenomicRanges)
library(Repitools)
files <- list.files("TFtables/mergedDuplicates/")
files <- files[grep(".merged.bed$", files)]
merged.data <- lapply(files, function(x) read.table(paste0("TFtables/mergedDuplicates/",x)))
for (i in 1:length(merged.data)) colnames(merged.data[[i]]) <- c("chr", "start", "end", "N")
load("Rdata/listOfAllRep_57.RData")
names(merged.data) <- names(listOfAllRep)
merged.data.gr <- lapply(merged.data, function(x) GRanges(x$chr, IRanges(x$start, x$end)))

doExpMatrix<-function(name){
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
a<-Sys.time()
# change mc.cores depending on your number of cores
expMatrix <- mclapply(names(listOfAllRep), doExpMatrix, mc.preschedule=F, mc.cores =32)
Sys.time()-a 
names(expMatrix)<-names(listOfAllRep)

# we will use that later
save(expMatrix, file="Rdata/binaryRepMatrix_57.RData")



