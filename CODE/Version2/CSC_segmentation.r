rm(list=ls())

source("FUNCTIONS_segmentation.r")
source("INPUT1_segmentation.r")

# loads packages
library(terra)
library(apcluster)

setwd("../Zone-correlation/")
file.remove(list.files())

set.seed(100)
setwd("../Phase-Classification/Samples/")
########## 0. PRELIMINARY CHECKS AND INFORMATION ##########

for(S_number in 1:length(Samples)){
  
Sample_name <- Samples[S_number] ; print(Sample_name)
setwd(Sample_name)
########## Extracting crystals of minimum size  ##########

table_in <- paste0(Samples[S_number],"_",PhaseName,".csv")
map_in <- paste0(Samples[S_number],"phases.tif")

phase_in <- read.csv(table_in)
mineral_map_in <- rast(map_in)

xy_points <- as.matrix(phase_in[,c("x","y")])
values <- phase_in[,Variables]

ncols <- ncol(mineral_map_in)
nrows <- nrow(mineral_map_in)

rasIn <- rast(ext(mineral_map_in), ncol=ncols, nrow=nrows)
rasIn <- rasterize(xy_points, rasIn, values)

p2 <- patches(which.lyr(rasIn>0),allowGaps=F)
minPixels <- max(minArea/(pix_um^2),(minS*2)^2)
crystals <- table(p2[])
tooSmall <- as.vector(which(crystals < minPixels))
if(length(tooSmall)>0){
tooSmall_ID <- which(p2[] %in% tooSmall)
rasIn[tooSmall_ID] <- NA
crystals <- crystals[-tooSmall]
}

areaDecreasing <- order(crystals,decreasing = T)
areaDecreasing_ID <- as.numeric(names(crystals)[areaDecreasing])
p2[] <- match(p2[],areaDecreasing_ID)

rm(crystals,tooSmall,tooSmall_ID,areaDecreasing,areaDecreasing_ID)

########## Iteratively processing each crystal  ##########

labelsZone <- NULL
textures <- rast(ext(mineral_map_in), ncol=ncols, nrow=nrows)
textures[] <- 0

X <- max(p2[],na.rm=T)

for(N in 1:X){

cat(paste("\nCrystal",N,"\n"))

########## Crop raster  ##########

xstalCoords <- xyFromCell(p2, which(p2[] == N))

xmin <- min(xstalCoords[,1]) ; xmax <- max(xstalCoords[,1])
ymin <- min(xstalCoords[,2]) ; ymax <- max(xstalCoords[,2])
xmin <- xmin - 0.5 ; ymin <- ymin - 0.5
xmax <- xmax + 0.5 ; ymax <- ymax + 0.5

e <- ext(x = c(xmin,xmax),y=c(ymin,ymax))
xstalInv <- crop(rasIn,e)
maskInv <- crop(p2,e)

rmBgXstal <- which(maskInv[]!=N | is.na(maskInv[])==TRUE)
if(length(rmBgXstal)>0) {
    xstalInv[rmBgXstal] <- NA
    maskInv[rmBgXstal] <- NA
    }

nonXstal <- which(is.na(maskInv[])==T)

rm(xmin,xmax,ymin,ymax,xstalCoords,rmBgXstal)

########## Calculate gradient, S & M  ##########

nrowsInv <- nrow(xstalInv)
ncolsInv <- ncol(xstalInv)

xstalInv <- multiScale(xstalInv)

gradient <- xstalInv
gradient <- focal(gradient,w=c(3,3),fun=mean,na.rm=T)
gradient <- xstalInv - gradient
gradient <- app(gradient, mean)
gradient <- (gradient^2)^0.5

gradientMean <- global(gradient, mean,na.rm=T)
gradientMean <- min(gradientMean*10,0.45)

S1 <- floor((nrowsInv+ncolsInv)^gradientMean)
S2 <- round((min(nrowsInv,ncolsInv))/3,0)
S <- min(S1,S2)
S <- max(S,minS)

searchRadius <- round(S*2,0)

xstalPix <- (nrowsInv*ncolsInv)-length(nonXstal)
M <- ((nrowsInv*ncolsInv) / xstalPix)^gradientMean

#Overwrite S if manual spacing is chosen
S <- ifelse(is.na(S_manual)==TRUE,max(S,minS),S_manual)

#Overwrite M if manual weighting is chosen
M <- ifelse(is.na(M_manual)==TRUE,M,M_manual)

rm(gradientMean,S1,S2,xstalPix)

########## Convert rasters to input matrices (flipped on y-axis)  ##########

gradient[nonXstal] <- 99

matIn <- array(xstalInv[],dim = c(ncolsInv,nrowsInv,nlyr(xstalInv)))
gradIn <- matrix(gradient[],nrow=ncolsInv,ncol=nrowsInv,byrow = F)

matIn <- aperm(matIn,c(2,1,3))
gradIn <- t(gradIn)

########## Set up grid of centroids  ##########

nr <- seq(0,(floor((ncolsInv-2)/S))*S,by = S); nr <- nr + floor((ncolsInv-max(nr))/2)
nc <- seq(0,(floor((nrowsInv-2)/S))*S,by = S); nc <- nc + floor((nrowsInv-max(nc))/2)
centroidGrid <- expand.grid(nc,nr)
centroidGrid <- matrix(unlist(centroidGrid),ncol=2,byrow = FALSE)

centroidGrid <- matrix(unlist(lapply(1:nrow(centroidGrid), centroidGrid.lowest.grad)),ncol=2,byrow=TRUE)
centroidGrid <- unique.matrix(centroidGrid)

nonXstalMat <- which(gradIn==99,arr.ind = TRUE)
nonXstalMat2 <- paste(nonXstalMat[,1],nonXstalMat[,2],sep="-")
centroidGridMat <- paste(centroidGrid[,1],centroidGrid[,2],sep="-")

nonXstalCentroid <- unique(match(nonXstalMat2,centroidGridMat))
nonXstalCentroid <- nonXstalCentroid[-which(is.na(nonXstalCentroid)==TRUE)]
if(length(nonXstalCentroid)!=0){
centroidGrid <- centroidGrid[-nonXstalCentroid,]}

rm(nr,nc,centroidGridMat,nonXstalMat2,nonXstalCentroid)

########## Initiate SLIC algorithm  ##########

#Defines the initial superpixel matrix (which superpixel each pixel belongs to)
l <- matrix(-1,nrowsInv,ncolsInv)
#Defines the initial minimum distance between each pixel and a centroid
d <- matrix(Inf,nrowsInv,ncolsInv)
#Defines all coordinates in grid
coords <- expand.grid(seq(1,nrowsInv,1),seq(1,ncolsInv,1))
xstalMat <- which(gradIn!=99)
coords <- coords[xstalMat,]

dSum <- rep(NA,nIters)

for (j in 1:nIters){
    tempM <- max(c(M/j,minM))
  for (i in 1:nrow(centroidGrid)){
    center <- centroidGrid[i,] ; rC <- center[1] ; cC <- center[2]
    r <- seq(rC-searchRadius,rC+searchRadius,1) ; r <- r[which(r>0)] ; r <- r[which(r<=nrowsInv)]
    c <- seq(cC-searchRadius,cC+searchRadius,1) ; c <- c[which(c>0)] ; c <- c[which(c<=ncolsInv)]
    sp <- expand.grid(r,c) ; sp <- as.matrix(sp)
    rP <- as.numeric(sp[,1])
    cP <- as.numeric(sp[,2])
    d1 <- as.numeric(unlist(lapply(seq(1,nrow(sp),1), function(x) { sum(unlist(lapply(1:length(Variables), function(N) {score.funct(arr = matIn,r1 = rC,c1 = cC,r2 = rP[x],c2 = cP[x],N=N) })))})))
    d1 <- d1^0.5
    d2 <- as.numeric(unlist(lapply(seq(1,nrow(sp),1), function(x) { ((rP[x]-rC)^2 + (cP[x]-cC)^2)^0.5 })))
    D <-  (d1^2 + ((d2^2)/S)*((tempM)^2))^0.5 
    test <- D < d[sp] 
    testTrue <- which(test==TRUE)
    l[as.matrix(sp[testTrue,])] <- i
    d[as.matrix(sp[testTrue,])] <- D[testTrue]
    rm(center,sp,rP,cP,d1,d2,D,test,testTrue)
  }

l[nonXstalMat] <- NA
l[as.matrix(coords)] <- unlist(lapply(1:nrow(coords), single.func,Mat=l))

duplicates <- which(duplicated.matrix(centroidGrid)==TRUE)
if (length(duplicates)>0){
  for (t in 1:length(duplicates)){
    d.replace <- which(centroidGrid[,1]==centroidGrid[duplicates[t],1] & centroidGrid[,2]==centroidGrid[duplicates[t],2])[1]
    l[which(l==duplicates[t],arr.ind = TRUE)] <- d.replace
    rm(d.replace)
  }
}

uniqueL <- sort(unique(as.vector(l)))
centroidGrid <- centroidGrid[uniqueL,]
l <- matrix(match(l,uniqueL),ncol=ncol(l),byrow = FALSE)  #This ensures empty indices are removed
l[nonXstalMat] <- -1

vectorL <- as.vector(l)
newCen <- lapply(1:nrow(centroidGrid),function(x) {
    idTemp <- which(vectorL[xstalMat]==x)
    coordsL <- as.matrix(coords[idTemp,])
    newVarTemp <- matrix(NA,nrow=nrow(coordsL),ncol=nlyr(xstalInv))
    for (o in 1:nrow(coordsL)){
        newVarTemp[o,] <- matIn[coordsL[o,1],coordsL[o,2],]
    }
    newCenVal <- apply(newVarTemp,2,mean) ; rm(newVarTemp)
    newCenCoord <- apply(coordsL,2,mean)
    newCenCoord <- ((sweep(coordsL,2,STATS=newCenCoord,FUN="-"))^2)
    newCenCoord <- (apply(newCenCoord,1,sum))^2
    newCenCoord <- coordsL[which.min(newCenCoord),]
    return(list(newCenCoord,newCenVal))
    rm(idTemp,coordsL)
    })

newCen <- matrix(unlist(newCen), ncol = 2+length(Variables), byrow = TRUE)
centroidGrid <- newCen[,c(1,2)]
meanSpsValue <- matrix(newCen[,-c(1,2)],ncol=length(Variables))
rm(newCen)

d[nonXstalMat] <- NA
dSum[j] <- sum(d,na.rm = TRUE)
if(j<4) {dDiff<- c(4,3,2,1)} else {dDiff <- c(NA,NA,NA,diff(diff(log(dSum))))}
stop.test <- dDiff[j+1]>dDiff[j]
if(stop.test==TRUE) {break}

# print(paste0(j," of ",N.Iters))

rm(duplicates,uniqueL,dDiff,vectorL,tempM)
}
rm(coords,d,xstalMat,dSum)

########## Perform AP algorithm & saves output ##########

l[which(l==-1)] <- NA
sps <- rast(l)
ext(sps) <- ext(xstalInv)
superpixelID <- sort(unique(as.vector(l)))
superpixelID <- superpixelID[which(superpixelID>0)]

meanSpsValue <- round(100*meanSpsValue,0)
diss <- matrix(,nrow=length(superpixelID),ncol=length(superpixelID))
for(d1 in 1:length(superpixelID)){
    cen1 <- meanSpsValue[d1,]
    for(d2 in 1:length(superpixelID)){
        cen2 <- meanSpsValue[d2,]
        cenDiff <- -1*(sum((cen1 - cen2)^2))
        diss[d1,d2] <- cenDiff
        rm(cen2,cenDiff)
}
    rm(cen1)
}

Ap <- apcluster(diss,q=q_choice) 
clusN <- length(Ap@exemplars)

Zone <- rast(ext(sps),ncol=ncolsInv, nrow=nrowsInv)
for (Ap_i in 1:clusN){
  Zone[which(sps[] %in% Ap@clusters[[Ap_i]])] <- Ap_i + minmax(textures)[2]
}

Zone <- extend(Zone,ext(textures)) ; crs(Zone) <- ""

Zone[which(is.na(Zone[])==T)] <- 0
textures <- textures + Zone 
for (j in 1:clusN){
  labelsZone <- c(labelsZone,paste0(Sample_name,"Zone",N,letters[j]))
}

rm(Ap,clusN,sps,superpixelID,Zone)

}

setwd("../../../Crystal-segmentation/")
writeRaster(textures,paste0(Sample_name,"zones.tif"),datatype = 'INT2S',overwrite=TRUE)
write.csv(x = labelsZone,file = paste0(Sample_name,"zones.csv"),row.names = F)

summaryOut <- list(PhaseName,pix_um,minArea,Variables,nIters,rel.scale,S_manual,minS,M_manual,minM,q_choice)
names(summaryOut) <- c("PhaseName","pix_um","minArea","Variables","nIters","rel.scale","S_manual","minS","M_manual","minM","q_choice")

if(S_number==length(Samples)){
sink(paste0("Summary.txt"))
print(summaryOut)
sink()
}

setwd("../Phase-Classification/Samples/")

}


