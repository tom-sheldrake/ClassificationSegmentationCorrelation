rm(list=ls())

require(flexmix)
require(mvtnorm)
require(terra)

source("INPUT1_classification.r")
elementsAll <- expand.grid(Elements1,Elements2)

setwd(paste0("../Data/",Sample_name))

elementsFiles <- list.files(pattern = ".txt")
elementsNames <- unlist(strsplit(elementsFiles,split = ".txt"))
elementsInID <- match(unique(c(Elements1,Elements2)),elementsNames)
elementsFiles <- elementsFiles[elementsInID]
elementsNames <- elementsNames[elementsInID]

invisible(lapply(elementsNames,function(x) {assign(x,read.table(paste0(x,".txt"),header=FALSE),envir = .GlobalEnv)}))
invisible(lapply(elementsNames, function(x) {assign(x,value = matrix(unlist(get(x)),nrow = nrow(get(x)),byrow = FALSE),envir = .GlobalEnv)}))

Results <- list()
Converged <- rep(NA,nrow(elementsAll))

for (i in 1:nrow(elementsAll)){
  matTemp1 <- get(as.vector(elementsAll$Var1[i]))
  matTemp2 <- get(as.vector(elementsAll$Var2[i]))
  X1 <- (1*matTemp1)+(2*matTemp2)
  X2 <- (2*matTemp1)+(1*matTemp2)
  dfAll <- cbind(as.vector(X1),as.vector(X2))
  phases <- flexmix(dfAll~1, k = nK, model = FLXmclust(diagonal = FALSE))
  Results[[i]] <-  eval(parse(text = "phases@cluster"))
  Converged[i] <- ifelse(phases@converged==TRUE,"Convergence","No convergence")
  print(paste(i,"of",nrow(elementsAll)))
  rm(matTemp1,matTemp2,X1,X2,dfAll,phases)
}

names(Results) <- paste(elementsAll$Var1,elementsAll$Var2,sep="-")
names(Converged) <- paste(elementsAll$Var1,elementsAll$Var2,sep="-")

Results <- as.data.frame(Results)

nrows <- nrow(get(elementsNames[1]))
ncols <- ncol(get(elementsNames[1]))

setwd("../../Phase-classification")
results.dir <- format(Sys.time(), "%d%b%Y_%Hhr%M")
dir.create(results.dir)
setwd(results.dir)

write.table(Converged,file="Convergence.txt")
write.csv(Results,"Results.csv")

save(list=c("Sample_name","Elements1","Elements2","nK","nrows","ncols"),file="input.Rdata")

plotRows <- round(nrow(elementsAll)^0.5,0)
plotCols <- ceiling(nrow(elementsAll)^0.5)

png(filename = paste0(Sample_name,"_fmm.png"),height = plotRows*1000,width = plotCols*1000)
layout(matrix(1:(plotRows*plotCols), nrow=plotRows, byrow=TRUE))
for(i in 1:ncol(Results)){
    matIn <- matrix(Results[,i],nrow=nrows,ncol=ncols)
    plot(rast(matIn),main=colnames(Results)[i])
}
dev.off()

inputOut <- list(Sample_name,Elements1,Elements2,nK)
names(inputOut) <- c("Sample","Elements 1","Elements 2","max clusters")
sink(paste0("input.txt"))
print(inputOut)
sink()