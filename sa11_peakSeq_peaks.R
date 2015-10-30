####
# comparison of SPP peaks with peakSeq peaks
###
setwd("myDirectory")
library(readr)
library(magrittr)
load("Rdata/listOfAllRep_57.RData")
reptable <- do.call(rbind, listOfAllRep)

# copied from http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/byDataType/peaks/jan2011/peakSeq/optimal/hub/ webpage
peakSeq <- read_tsv("peakSeq_table.txt")

commands <- paste0("wget ", peakSeq$peakSeq_link)

write.table(as.data.frame(commands), file = "TFtables/peakSeq_peaks/download_peaks.sh", quote = FALSE, row.names = FALSE, col.names = FALSE)
# wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed

commands <- paste0("./bigBedToBed ", peakSeq$peakSeq_fn, " ", sub(".bb", ".bed", peakSeq$peakSeq_fn, fixed = TRUE))
write.table(as.data.frame(commands), file = "TFtables/peakSeq_peaks/convert2bed.sh", quote = FALSE, row.names = FALSE, col.names = FALSE)

peakSeq$peakSeq_fn <- sub(".bb", "", peakSeq$peakSeq_fn, fixed = TRUE)

listOfRep <- lapply(unique(peakSeq$name), function(x) subset(peakSeq, peakSeq$name == x))
names(listOfRep) <- unique(peakSeq$name)

printMergeCommands<-function(ln) {
    x <- listOfRep[[ln]]$peakSeq_fn
    vofn <- paste0(x, ".bed")
    files<-""
    for (i in 1:length(vofn)) {
        files<-paste(files, vofn[i])
    }
    command1<-paste0("cat ", files, " > mergedDuplicates/", ln, ".peakSeq.cat.bed")
    command2<-paste0("sortBed -i mergedDuplicates/", ln, ".peakSeq.cat.bed > mergedDuplicates/", ln, ".peakSeq.s.bed")
    command3<-paste0("bedtools merge -c 7 -o count -i mergedDuplicates/", ln, ".peakSeq.s.bed > mergedDuplicates/", ln, ".peakSeq.merged.bed")
    return(rbind(command1, command2, command3))
}

commands<-lapply(names(listOfRep), printMergeCommands)
formatedCommands<-""
for (i in 1:length(commands)) {
    formatedCommands<-rbind(formatedCommands, commands[[i]])
}

write.table(formatedCommands, file="TFtables/peakSeq_peaks/merge_peakSeq_replicates.sh", quote = FALSE, row.names = FALSE, col.names = FALSE)

###
library(GenomicRanges)
files <- list.files("TFtables/peakSeq_peaks/mergedDuplicates/")
files <- files[grep(".merged.bed$", files)]

merged.data <- lapply(files, function(x) read.table(paste0("TFtables/peakSeq_peaks/mergedDuplicates/", x)))
for (i in 1:length(merged.data)) colnames(merged.data[[i]]) <- c("chr", "start", "end", "N")
names(merged.data) <- sub(".peakSeq.merged.bed", "", files, fixed = TRUE)

merged.data.gr <- lapply(merged.data, function(x) GRanges(x$chr, IRanges(x$start, x$end)))

doExpMatrix <- function(name){
    table <- listOfRep[[name]]
    TFfiles <- table$peakSeq_fn
    TFfiles <- paste0(TFfiles, ".bed")
    TFdata <- lapply(TFfiles, function(x) read.table(paste0("TFtables/peakSeq_peaks/" ,x)))
    names(TFdata) <- table$ID_spp
    TFdata.gr <- lapply(TFdata, function(x) GRanges(x[,1], IRanges(x[,2], x[,3])))
    hits <- lapply(TFdata.gr, function(x) findOverlaps(x, merged.data.gr[[name]], select="first"))
    expMatrix <- matrix(data = 0, nrow = nrow(merged.data[[name]]), ncol = nrow(table))
    colnames(expMatrix) <- table$name
    for(i in 1:nrow(table)){
        expMatrix[hits[[i]],i] <- 1
    }
    return(cbind(merged.data[[name]], expMatrix))
}

library(parallel)
a <- Sys.time()
expMatrix_peakSeq <- mclapply(names(listOfRep), doExpMatrix, mc.preschedule = FALSE, mc.cores = 12)
Sys.time() - a 
names(expMatrix_peakSeq) <- names(listOfRep)

save(expMatrix_peakSeq, file = "Rdata/binaryRepMatrix_peakSeq.RData")


