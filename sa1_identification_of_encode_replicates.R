# Go to your working directory
setwd("myDirectory/")
# The following file was manually buid :-(
# it should be easy to be rebuild from:
# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/
# and
# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/
metaData <- read.table("TFtables/metadata_TFtables.txt")
metaData <- metaData[,2:6]
colnames(metaData) <- c("lab", "cell", "treatment", "antibody", "name")
levels(metaData$lab)
levels(metaData$treatment)
levels(metaData$cell)
levels(metaData$antibody)

temp <- sapply(as.character(metaData$antibody), function(x) strsplit(x, "_"))
temp2 <- sapply(temp, function(x) x[[1]])
metaData <- cbind(metaData, "TF"=temp2)
rm(temp, temp2)

# manual curration of ENCODE at UCSC names...
metaData[which(metaData$TF=="GATA-2"), "TF"] <- "GATA2"
metaData[which(metaData$TF=="GRp20"), "TF"] <- "GR"
metaData[which(metaData$TF=="PAX5-C20"), "TF"] <- "PAX5"
metaData[which(metaData$TF=="PAX5-N19"), "TF"] <- "PAX5"
metaData[which(metaData$TF=="Sin3Ak-20"), "TF"] <- "SIN3A"
metaData[which(metaData$TF=="USF-1"), "TF"] <- "USF1"
metaData[which(metaData$TF=="Pol2(b)"), "TF"] <- "Pol2_b"
metaData[which(metaData$TF=="Pol2(phosphoS2)"), "TF"] <- "Pol2_phosphoS2"

for (i in 1:nrow(metaData)) {
    metaData[i,"trueRep"] <- paste0(metaData[i,"cell"],"_", metaData[i,"treatment"], "_", metaData[i,"antibody"])
    metaData[i,"allRep"] <- paste0(metaData[i,"cell"],"_", metaData[i,"treatment"], "_", metaData[i,"TF"])
}

rep<-list(
    "trueRep" = table(metaData$trueRep),
    "allRep" = table(metaData$allRep)
)
sapply(rep, length)
rep <- lapply(rep, function(x) subset(x, x >= 2))
sapply(rep, length)

# trueRep  allRep
# 41      57

namesRep <- lapply(rep, names)

listOfAllRep <- lapply(namesRep$allRep, function(x) subset(metaData, metaData[,"allRep"]==x))
sapply(listOfAllRep, nrow)
peakNamesList <- lapply(listOfAllRep, function(x) as.character(x[,"name"]))

sapply(peakNamesList, length)
levels(factor(sapply(peakNamesList, length)))

tempNames <- namesRep$allRep
tempNames[which(tempNames=="K562_None_Pol2(phosphoS2)")] <- "K562_None_Pol2_phosphoS2"
names(peakNamesList) <- tempNames

# This function print a bedtools command for merging repplicates .bed files
# It requires SORTED .bed files of peaks in the "sortedTFtables/" directory (gunzip + sortBed)
# These bed files are sorted version of the files found in:
# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/
# It also require an "mergedDuplicates/" directory for writing output files
printMergeCommands <- function(ln) {
    x<-peakNamesList[[ln]]
    vofn <- paste0("sortedTFtables/" , x, ".narrowPeak.s")
    files<-""
    for (i in 1:length(vofn)) {
        files <- paste(files, vofn[i])
    }
    command1 <- paste0("cat ", files, " > mergedDuplicates/", ln, ".cat.bed")
    command2 <- paste0("sortBed -i mergedDuplicates/", ln, ".cat.bed > mergedDuplicates/", ln, ".s.bed")
    command3 <- paste0("bedtools merge -c 7 -o count -i mergedDuplicates/", ln, ".s.bed > mergedDuplicates/", ln, ".merged.bed")
    return(rbind(command1, command2, command3))
}

commands<-lapply(tempNames, printMergeCommands)
formatedCommands <- ""
for (i in 1:length(commands)) {
    formatedCommands<-rbind(formatedCommands, commands[[i]])
}
write.table(formatedCommands, file = "TFtables/merge_encode_replicates.sh", quote = FALSE, row.names = FALSE, col.names = FALSE)
# now we give the .sh file executable rights and we execute it
# require bedtools

names(listOfAllRep) <- names(rep$allRep)
save(listOfAllRep, file="Rdata/listOfAllRep_57.RData")
# this RData file will be needed for the rest of the analysis
