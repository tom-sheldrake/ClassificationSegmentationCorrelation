rm(list=ls())

require(terra)

source("INPUT1_correlation.r")

#Rearranges Samples in alphabetical order
SamplesAll <- sort(SamplesAll,decreasing = FALSE)

setwd("../Phase-classification/Samples/")
for (i in 1:length(SamplesAll)){
  assign(paste0(SamplesAll[i],"_",PhaseName),read.csv(paste0(SamplesAll[i],"/",SamplesAll[i],"_",PhaseName,".csv")))
  print(i)
}

setwd("../../Crystal-segmentation//")
for (i in 1:length(SamplesAll)){
  assign("rasTemp",rast(paste0(SamplesAll[i],"zones.tif")))
  assign("labelsTemp",read.csv(paste0(SamplesAll[i],"zones.csv")))
  labelsTemp <-cbind(1:nrow(labelsTemp),labelsTemp)
  colnames(labelsTemp) <- c("ID","Zone")
  levels(rasTemp) <- labelsTemp
  names(rasTemp) <- "Labels"
  labels.id <- labelsTemp$x
  assign(paste0(SamplesAll[i],"_textures"),rasTemp)
  dfTemp1 <- as.data.frame(rasTemp,xy=TRUE)
  dfTemp2 <- get(paste0(SamplesAll[i],"_",PhaseName))
  assign(paste0(SamplesAll[i],"_",PhaseName),merge(dfTemp2,dfTemp1,by = c("x","y")))
  rm(rasTemp,labelsTemp)
  print(i)
}


all_zones <- lapply(ls(pattern = "_textures",envir = sys.frame()), function(x) {temp <- levels(get(x))[[1]] ; return(temp)})
names(all_zones) <- SamplesAll
N <- unlist(lapply(all_zones,nrow))
zone_ID <- lapply(seq_along(SamplesAll), function(n,i) {lapply(eval(parse(text=paste0("all_zones$",n[[i]],"$Zone"))), function(x,n2) {which(eval(parse(text=paste0(n2,"_",PhaseName,"$Labels")))==x)},n2=n[[i]])},n=SamplesAll)

zone_ID <- unlist(zone_ID,recursive = FALSE)
names(zone_ID) <- as.vector(unlist(lapply(1:length(all_zones), function(x) {all_zones[[x]][,2]})))

correl_vars <- lapply(seq_along(SamplesAll), function(S) {
    dfTemp <- eval(parse(text=paste0(SamplesAll[S],"_",PhaseName)))
    idTemp <- match(Variables,colnames(dfTemp))
    return(as.matrix(dfTemp[,idTemp],ncol=length(Variables)))
    rm(dfTemp,idTemp)
})
names(correl_vars) <- SamplesAll


zoneVars <- lapply(1:length(zone_ID),function(Z){
    zoneTemp <- names(zone_ID)[Z]
    sampleTemp <- strsplit(zoneTemp,split = "Zone")[[1]][1]
    matTemp <- eval(parse(text=paste0("correl_vars$",sampleTemp)))
    corVar <- matTemp[as.vector(unlist(zone_ID[Z])),]
    return(matrix(corVar,ncol=length(Variables)))
    rm(zoneTemp,sampleTemp,matTemp,corVar)
    })
names(zoneVars) <- names(zone_ID)

minTemp <- lapply(1:length(zoneVars), function(Z) {apply(zoneVars[[Z]],2,min)})
minTemp <- matrix(unlist(minTemp),ncol=length(Variables),byrow = T)
maxTemp <- lapply(1:length(zoneVars), function(Z) {apply(zoneVars[[Z]],2,max)})
maxTemp <- matrix(unlist(maxTemp),ncol=length(Variables),byrow = T)
limits <- rbind(apply(minTemp,2,min),apply(maxTemp,2,max))
if(rel.scale==T){
    limits[1,] <- min(limits[1,])
    limits[2,] <- max(limits[2,])
}

correlVar <- lapply(1:length(zoneVars), function(Z) {
  matTemp <- zoneVars[[Z]]
  matOut <- sweep(matTemp,MARGIN=2,STATS = limits[1,],FUN="-")
  matOut <- sweep(matOut,MARGIN=2,STATS = (limits[2,]-limits[1,]),FUN="/")
  matOut <- apply(matOut,1,mean)
  matOut <- round(100*matOut,0)
  return(matOut)
})
names(correlVar) <- names(zone_ID)

xSeq <- seq(0,100,by=0.5)
correlMat <- matrix(,nrow=length(zoneVars),ncol = length(xSeq))
rownames(correlMat) <- names(zoneVars)

for (i in 1:length(zoneVars)){
  correlMat[i,] <-  unlist(lapply(1:length(xSeq), function(x) {length(which(correlVar[[i]] < xSeq[x]))})) / length(correlVar[[i]])
}

setwd("../Zone-correlation")
results.dir <- format(Sys.time(), "%d%b%Y_%Hhr%M")
dir.create(results.dir)

#Saves distance matrix
setwd(results.dir)
write.table(correlMat,"correlation_matrix.txt",row.names=TRUE,col.names=TRUE)

distMat<- dist(correlMat, method = "euclidean", diag=FALSE)
hc <- hclust(d = distMat,method = "ward.D2")
C_IND <- NULL
Smin <- NULL
Smax <- NULL
SW <- NULL
a <- 3 ; b <- 20
for(i in a:b){
  ind.tmp <- i-a+1
  memb <- cutree(hc, k = i)
  N_memb <- as.numeric(table(memb))
  NW <- (N_memb*(N_memb-1))/2
  NW <- sum(NW)
  S <- sort(as.numeric(distMat),decreasing = FALSE)[1:NW] ; Smin[ind.tmp] <- sum(S)
  S <- sort(as.numeric(distMat),decreasing = TRUE)[1:NW] ; Smax[ind.tmp] <- sum(S)
  within.clus <- NULL
  for(j in 1:i){
    id.temp <- which(memb==j)
    coords<- expand.grid(id.temp,id.temp)
    within.clus[j] <- sum(distMat[as.matrix(coords)])
  }
  SW[ind.tmp] <- sum(within.clus)
  C_IND[ind.tmp] <- (SW[ind.tmp] - Smin[ind.tmp]) / (Smax[ind.tmp] - Smin[ind.tmp])
}

png(filename = paste0("C-Index.png"),width = 1000,height = 500,units = "px")
par(mfrow=c(1,1))
plot(seq(a,b,1),C_IND,type="b",lwd=2,cex=2,main="C-Index",ylab="C-Index",xlab="Zoning groups")
dev.off()

summaryOut <- list(SamplesAll,PhaseName,Variables,rel.scale,results.dir )
names(summaryOut) <- c("SamplesAll","PhaseName","Variables","rel.scale","Results Dir")

save(summaryOut,file="summary.Rdata")






