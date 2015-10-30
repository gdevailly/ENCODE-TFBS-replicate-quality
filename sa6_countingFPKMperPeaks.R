# for each peaks in each peak files, we extract FPKM value from .bam files.
# the script is a bit messy, notably because one .bam fil from ENCODE seems ill formated and had to be treaten differently
# (H1hescJundv0416102)

setwd("myDirectory/")
library(GenomicRanges)
files <- list.files("TFtables/mergedDuplicates/")
files <- files[grep(".merged.bed$", files)]
merged.data <- lapply(files, function(x) read.table(paste0("TFtables/mergedDuplicates/", x)))
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
a<-Sys.time()
expMatrix <- mclapply(names(listOfAllRep), doExpMatrix, mc.preschedule=F, mc.cores = 32)
Sys.time() - a # 10s <3 <3 <3
names(expMatrix) <- names(listOfAllRep)

library(Repitools)

loadBam <- function(expName, repMatrix){
    bamGR <- BAM2GRanges(paste0("TFtables/bamReps/", expName))
    chr <- levels(factor(repMatrix$chr))
    bamGR <- bamGR[seqnames(bamGR) %in% chr]
    counts <- annotationBlocksCounts(bamGR, repMatrix[,1:3], seq.len=300)
    counts <- cbind(counts, (1E6*counts)/length(bamGR))
    colnames(counts) <- c("nReads", "FPM")
    counts <- cbind(counts, "FPKM" = (1000*counts[,"FPM"])/(repMatrix$end-repMatrix$start))
    dataList <- list("countsTable" = counts, "nReads" = length(bamGR))
    return(dataList)
}

extractDataFromBams <- function(dupName){
    tempT <- listOfAllRep[[dupName]]
    bamNames <- paste0(tempT$name, ".bam")
    repMatrix <- cbind(expMatrix[[dupName]])
    listOfDataList <- lapply(bamNames, function(x) loadBam(x, repMatrix))
    names(listOfDataList) <- tempT$name
    nReads <- vector(length = length(listOfDataList), mode = "numeric")
    names(nReads) <- tempT$name
    for(i in tempT$name) {
        nReads[i] <- listOfDataList[[i]][["nReads"]]
        repMatrix[,i] <- listOfDataList[[i]][["countsTable"]][,"FPKM"]
    }
    message(paste0(dupName, " done!"))
    return(list("nReads" = nReads, "FPKMtable" = repMatrix))
}

a <- Sys.time()
listOfFPKMtable <- mclapply(names(listOfAllRep), extractDataFromBams, mc.preschedule = FALSE, mc.cores = 10)
Sys.time() - a # 45 minutes, 1 error

names(listOfFPKMtable) <- names(listOfAllRep)


listOfFPKMtable[["H1-hESC_None_JunD"]]
# [1] "Error in GAlignments(seqnames = bamcols$rname, pos = bamcols$pos, cigar = bamcols$cigar,  : \n  'seqnames' cannot have NAs\n"
# attr(,"class")
# [1] "try-error"
# attr(,"condition")
# <simpleError in GAlignments(seqnames = bamcols$rname, pos = bamcols$pos, cigar = bamcols$cigar,     strand = bamcols$strand, seqlengths = seqlengths): 'seqnames' cannot have NAs>
#     

save(listOfFPKMtable, file = "Rdata/temp.listOfFPKMtable.repplicates_57.RData")


####
# H1hescJundv0416102
# the bam file is somehow ill formated and had to be treaten differently from the rest
# wgEncodeAwgTfbsHaibH1hescJundV0416102UniPk.sam.cut was generated from wgEncodeAwgTfbsHaibH1hescJundV0416102UniPk.bam
# by converting the .bam file in .sam zith samtools
# and keeping only flag, chr and start (comuns 2, 3 and 4) using awk
####

setwd("myDirectory/")
library(GenomicRanges)
library(Repitools)

Fsam <- read.table("TFtables/bamReps/wgEncodeAwgTfbsHaibH1hescJundV0416102UniPk.sam.cut", stringsAsFactors = FALSE, sep = "\t")
colnames(Fsam) <- c("flag", "chr", "start")
Fsam$end <- Fsam$start
Fsam$strand <- "*"
table(Fsam$flag)
plus <- which(Fsam$flag==0)
minus <- which(Fsam$flag==16)
Fsam[plus, "end"] <- Fsam[plus, "start"] + 50
Fsam[plus, "strand"] <- "+"
Fsam[minus, "start"] <- Fsam[minus, "end"] - 50
Fsam[minus, "strand"] <- "-"

FsamGood <- subset(Fsam, Fsam$flag != 4)
FsamGR <- with(Fsam, GRanges(chr, IRanges(start, end), "strand" = strand))
FotherSample <- BAM2GRanges("TFtables/bamReps/wgEncodeAwgTfbsSydhH1hescJundIggrabUniPk.bam")
Fexp <- listOfAllRep[["H1-hESC_None_JunD"]]
FGRL <- GRangesList(FsamGR, FotherSample)
names(FGRL) <- Fexp$name

### FPKMtable

FrepMatrix <- cbind(expMatrix[["H1-hESC_None_JunD"]])
chr <- levels(factor(FrepMatrix$chr))
FGRL[[1]] <- FGRL[[1]][seqnames(FGRL[[1]]) %in% chr]
FGRL[[2]] <- FGRL[[2]][seqnames(FGRL[[2]]) %in% chr]
Fcounts <- annotationBlocksCounts(FGRL, FrepMatrix[,1:3], seq.len = 300)
Fcounts <- cbind(Fcounts, 1E6*Fcounts[,1]/length(FGRL[[1]]), 1E6*Fcounts[,2]/length(FGRL[[2]]))              
colnames(Fcounts) <- c("nReads_1", "nReads_2",  "FPM_1", "FPM_2")
Fcounts <- cbind(Fcounts, "FPKM_1" = (1000*Fcounts[,"FPM_1"])/(FrepMatrix$end - FrepMatrix$start), "FPKM_2" = (1000*Fcounts[,"FPM_2"])/(FrepMatrix$end - FrepMatrix$start))
FFPKMtable <- list(
    "nReads" = sapply(FGRL, length),
    "FPKMtable" = cbind(FrepMatrix[,1:4], Fcounts[,5:6])
)
colnames(FFPKMtable$FPKMtable) <- c("chr", "start", "end", "N", "wgEncodeAwgTfbsHaibH1hescJundV0416102UniPk", "wgEncodeAwgTfbsSydhH1hescJundIggrabUniPk")

listOfFPKMtable[["H1-hESC_None_JunD"]] <- FFPKMtable
save(listOfFPKMtable, file = "Rdata/temp.listOfFPKMtable.repplicates_57_gh.RData")

### List of CovMat
# it will be merged later: see sa8_meanCoverageWithRepitools

coord <- with(FrepMatrix, GRanges(chr, IRanges(start, end)))
FfsList <- featureScores(FGRL, coord, up = 2000, down = 2000, freq = 20, s.width = 300) #longlonglong
FfsTable <- tables(FfsList)
FfsTable <- lapply(FfsTable, function(x) x*1E6)

covMat<-list()
for(i in 1:length(FfsTable)){
    covMat[[i]] <- matrix(nrow=ncol(FfsTable[[i]]), ncol = 3)
    rownames(covMat[[i]]) <- colnames(FfsTable[[i]])
    colnames(covMat[[i]]) <- c("common", "not common", "undetected")
    covMat[[i]][,"common"] <- colMeans(subset(FfsTable[[i]], rowSums(FrepMatrix[5:ncol(FrepMatrix)]) == 2))
    covMat[[i]][,"not common"] <- colMeans(subset(FfsTable[[i]], FrepMatrix[, 4+i] == 1 & rowSums(FrepMatrix[5:ncol(FrepMatrix)]) != 2))
    covMat[[i]][,"undetected"] < -colMeans(subset(FfsTable[[i]], FrepMatrix[, 4+i] == 0))
}

names(covMat) <-c ("wgEncodeAwgTfbsHaibH1hescJundV0416102UniPk", "wgEncodeAwgTfbsSydhH1hescJundIggrabUniPk")
save(covMat, file="Rdata/temp.listOfCovMat.repplicates.ESCJun.RData")
