rm(list=ls())

require(terra)

source("INPUT2_classification.r")

setwd("../Phase-classification/")

file.remove(list.files(pattern=".csv"))
file.remove(list.files(pattern=".tif"))
file.remove(list.files(pattern=".txt"))
file.remove(list.files(pattern=".png"))


if(is.na(folderChoice)==T){
dfTemp <- file.info(list.files())
dfTemp <- dfTemp[-which(rownames(dfTemp)=="Samples"),]
classFolder <- rownames(dfTemp)[which.max(dfTemp$mtime)]
rm(dfTemp)
} else {
classFolder <- folderChoice
}
setwd(classFolder)

load("input.Rdata")

Results <- read.csv("Results.csv",row.names=1)

setwd("../../Data/")
setwd(Sample_name)

elementsFiles <- list.files(pattern = ".txt")
elementsNames <- unlist(strsplit(elementsFiles,split = ".txt"))
invisible(lapply(elementsNames,function(x) {assign(x,read.table(paste0(x,".txt"),header=FALSE),envir = .GlobalEnv)}))
dfOut <- as.data.frame(lapply(elementsNames,function(x) {as.vector(unlist(t(get(x))))}))
colnames(dfOut) <- elementsNames

ResultsFinal <- Results[,match(modelChoice,colnames(Results))]
ResultsFinalUnique <- unique(ResultsFinal)
ResultsFinalUnique$N <- 1:nrow(ResultsFinalUnique)
ResultsFinal$X <- 1:nrow(ResultsFinal)
distinctClus <- merge(ResultsFinal,ResultsFinalUnique,by=modelChoice)
distinctClus <- distinctClus[order(distinctClus$X),]

clusID <- matrix(distinctClus$N,ncol = ncols)
clusIDRas <- rast(clusID)
clusIDdf <- as.data.frame(clusIDRas,xy=TRUE)
colnames(clusIDdf)[3] <- "N"

dfOut <- cbind.data.frame(clusIDdf,dfOut)

varRas <- focal(clusIDRas,w=c(3,3),fun=var)
components <- as.vector(unlist(unique(clusIDRas[which(varRas[]==0)])))

dfOut$N <- match(dfOut$N,components)
clusID <- t(matrix(dfOut$N,ncol = ncols))

setwd("../../Phase-classification/")

png(filename = paste0(Sample_name,"_phases.png"),width=1000,height=1000,units="px")
plot(rast(clusID),col=rainbow(length(components)))
dev.off()

writeRaster(rast(clusID),paste0(Sample_name,"phases.tif"),datatype = 'INT2S',overwrite=TRUE)

summaryOut <- list(classFolder,modelChoice)
names(summaryOut) <- c("Input folder","Chosen models")
sink(paste0(Sample_name,"_summary1.txt"))
print(summaryOut)
sink()

write.csv(dfOut,file=paste0(Sample_name,".csv"),row.names=F)
