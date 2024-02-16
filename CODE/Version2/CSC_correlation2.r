rm(list=ls())

library(terra)

source("INPUT2_correlation.r")

setwd("../Zone-correlation/")

file.remove(list.files(pattern=".csv"))
file.remove(list.files(pattern=".tif"))
file.remove(list.files(pattern=".txt"))

if(is.na(folderChoice)==T){
dfTemp <- file.info(list.files())
correlFolder <- rownames(dfTemp)[which.max(dfTemp$mtime)]
} else {
correlFolder <- folderChoice
}
setwd(correlFolder)
rm(dfTemp)

load("summary.Rdata")

SamplesAll <- summaryOut$SamplesAll
PhaseName <- summaryOut$PhaseName
Variables <- summaryOut$Variables

#Rearranges Samples in alphabetical order
SamplesAll <- sort(SamplesAll,decreasing = FALSE)

correlMat <- read.table("correlation_matrix.txt",header=T)
setwd("../")

distMat<- dist(correlMat, method = "euclidean", diag=FALSE)
hc <- hclust(d = distMat,method = "ward.D2")
ht <- cutree(hc,k=nClus)

setwd("../Phase-Classification/Samples/")
for (i in 1:length(SamplesAll)){
  assign(paste0(SamplesAll[i],"_",PhaseName),read.csv(paste0(SamplesAll[i],"/",SamplesAll[i],"_",PhaseName,".csv")))
  print(i)
}

zg_df <- list()
setwd("../../Crystal-segmentation//")
for (i in 1:length(SamplesAll)){
  assign("rasTemp",rast(paste0(SamplesAll[i],"zones.tif")))
  assign("labelsTemp",read.csv(paste0(SamplesAll[i],"zones.csv")))
  labelsTemp <-cbind(1:nrow(labelsTemp),labelsTemp)
  colnames(labelsTemp) <- c("ID","Zone")
  levels(rasTemp) <- labelsTemp
  names(rasTemp) <- "Labels"
  labels.id <- labelsTemp$x
  dfTemp1 <- as.data.frame(rasTemp,xy=TRUE)
  dfTemp2 <- get(paste0(SamplesAll[i],"_",PhaseName))
  dfTemp3 <- merge(dfTemp2,dfTemp1,by = c("x","y"))
  dfTemp3$Zone <- as.vector(ht[match(dfTemp3$Labels,names(ht))])

  zg_df[[i]] <- dfTemp3

  xy_points <- as.matrix(dfTemp3[,c("x","y")])
  rasOut <- rast(ext(rasTemp), ncol=ncol(rasTemp), nrow=nrow(rasTemp))
  rasOut <- rasterize(xy_points, rasOut, dfTemp3$Zone)
  
  setwd("../Zone-correlation/")
  writeRaster(rasOut,paste0(SamplesAll[i],"_ZG",nClus,".tif"),datatype = 'INT2S',overwrite=TRUE)
  setwd("../Crystal-segmentation//")
  rm(rasTemp,labelsTemp,dfTemp1,dfTemp2,dfTemp3,xy_points,rasOut)
  print(i)
}

setwd("../Zone-correlation/")
zg_df <- as.data.frame(do.call(rbind,zg_df))
write.csv(x = zg_df,file = paste0("ZG_",nClus,".csv"),row.names = F)

summaryOut <- append(summaryOut,nClus)
names(summaryOut)[6] <- "nClus"

sink(paste0("ZG_",nClus,"_summary.txt"))
print(summaryOut)
sink()