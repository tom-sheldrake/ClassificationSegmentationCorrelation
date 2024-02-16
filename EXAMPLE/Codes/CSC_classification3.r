rm(list=ls())

source("INPUT3_classification.r")

setwd("../Phase-classification")

dfIn <- read.csv(list.files(pattern = ".csv"))
sample_name <- unlist(strsplit(list.files(pattern = ".csv"),split=".csv"))

dfOut <- dfIn[which(dfIn$N==phaseNumber),]

if(calibration==TRUE) {
    invisible(lapply(1:length(newVars),function(V) {
        tempFunc <- get(paste0(newVars[[V]],".function"))
        tempVarID <- match(calibrationVars[[V]],colnames(dfOut))
        tempVar <- dfOut[,tempVarID]
        assign(newVars[[V]],tempFunc(tempVar),envir = .GlobalEnv)
    }
    ))
}

if(calibration==TRUE) {
  dfOut <- cbind.data.frame(dfOut,lapply(newVars,get))
  a <- ncol(dfOut)-length(newVars) + 1
  colnames(dfOut)[a:ncol(dfOut)] <- unlist(newVars)
}

write.csv(dfOut,file=paste0(sample_name,"_",phaseName,".csv"))
file.remove(paste0(sample_name,".csv"))

summaryOut <- list(phaseName,phaseNumber)
names(summaryOut) <- c("Phase name","Phase number")
sink(paste0(sample_name,"_summary2.txt"))
print(summaryOut)
if(calibration==TRUE){
    summaryTemp <- list(newVars,calibrationVars,lapply(ls(pattern = ".function"),get))
    names(summaryTemp) <- c("New variables", "Calibration variables","Calibration function")
    print(summaryTemp)
}
sink()

dir.create(paste0("Samples/",sample_name),showWarnings = F)
filesSample <- list.files(pattern = sample_name)
invisible(lapply(filesSample,file.copy,to=paste0("Samples/",sample_name),overwrite=TRUE))
invisible(lapply(filesSample,file.remove))