# writing commands for finding de novo motif discovery
# require homer
# http://homer.salk.edu/homer/
setwd("myDirectory/")
load("Rdata/listOfAllRep_57.RData")

prefixInput <- "TFtables/sortedTFtables/"
prefixOutput <- "homerMotifOutput/individualTF/"

# require "homerMotifOutput/individualTF/" directory
# modify -p parameter depending on number of cores you ant to use

makeHomerCommands<-function(name) {
    tempT <- listOfAllRep[[name]]
    commands <- data.frame()   
    commands[1:nrow(tempT),1] <- paste0("findMotifsGenome.pl ", prefixInput, tempT$name,".narrowPeak.s hg19 ", prefixOutput, tempT$name,
                                      ".narrowPeak/ -size given -p 20 -S 20 -preparsedDir tempHomer")
    return(commands)
}

listOfCommands <- lapply(names(listOfAllRep), makeHomerCommands)

com <- listOfCommands[[1]]
for (i in 2:length(listOfCommands)) {
    com<-rbind(com, listOfCommands[[i]])
}

write.table(com, file="homerOnOthersEncodeRep.sh", quote=F, row.names=F, col.names=F)
