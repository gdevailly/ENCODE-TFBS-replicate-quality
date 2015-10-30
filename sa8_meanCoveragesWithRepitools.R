# we compute mean coverage at various peaks from .bam files, using Repitools package
setwd("myDirectory/")
library(GenomicRanges)
library(Repitools) # Bioconductor package
load("Rdata/binaryRepMatrix_57.RData")
# load("Rdata/listOfFPKMtable.repplicates_57.RData")
load("Rdata/listOfAllRep_57.RData")
categ<-read.table("classificationOfEncodeReplicates_57.txt")
colnames(categ) <- c("name", "category")
rownames(categ) <- categ$name

calcCoverage <- function(name){
    metaD <- listOfAllRep[[name]]
    binMat <- expMatrix[[name]]
    coord <- with(binMat, GRanges(chr, IRanges(start, end)))
    fsList <- list()
    for(i in 1: nrow(metaD)){
        bam <- paste0("TFtables/bamReps/", metaD[i,"name"], ".bam")
        fsList[[i]] <- featureScores(bam, coord, up=2000, down=2000, freq=20, s.width=300) #longlonglong
    }
    names(fsList) <- metaD$name
    fsTable <- lapply(fsList, tables)
    for(i in 1:length(fsTable)){
        fsTable[[i]] <- fsTable[[i]][[1]]
    }
    fsTable <- lapply(fsTable, function(x) x*1E6)
    covMat <- list()
    for(i in 1:length(fsTable)){
        covMat[[i]] <- matrix(nrow = ncol(fsTable[[i]]), ncol = 3)
        rownames(covMat[[i]]) <- colnames(fsTable[[i]])
        colnames(covMat[[i]]) <- c("common", "not common", "undetected")
        covMat[[i]][,"common"] <- colMeans(subset(fsTable[[i]], rowSums(binMat[5:ncol(binMat)]) == nrow(metaD)))
        covMat[[i]][,"not common"] <- colMeans(subset(fsTable[[i]], binMat[, 4+i] == 1 & rowSums(binMat[5:ncol(binMat)]) != nrow(metaD)))
        covMat[[i]][,"undetected"] <- colMeans(subset(fsTable[[i]], binMat[, 4+i] == 0))
    }
    return(covMat)
}

library(parallel)
a <- Sys.time()
listOfCovMat <- mclapply(names(listOfAllRep), calcCoverage, mc.preschedule = FALSE, mc.cores = 10)
Sys.time() - a # ~ 1 hour
# Warning message:
# In mclapply(names(listOfAllRep), calcCoverage, mc.preschedule = FALSE,  :
#                 1 function calls resulted in an error
# this is again due to H1-hESC_None_JunD .bam file ill formating             
names(listOfCovMat) <- names(listOfAllRep)

# we load precompute data for H1-hESC_None_JunD (see sa6_countingFPKMperPeaks.R)
load(file="Rdata/temp.listOfCovMat.repplicates.ESCJun.RData")
listOfCovMat[["H1-hESC_None_JunD"]] <- covMat
save(listOfCovMat, file="Rdata/listOfCovMat.replicates_57.RData")



