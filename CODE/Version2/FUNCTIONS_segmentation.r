multiScale <- function(x) {
  n <- seq(nlyr(x))
  limits <- minmax(x)
  if(rel.scale==T){
    limits[1,] <- min(limits[1,])
    limits[2,] <- max(limits[2,])
  }
  list1 <- lapply(X=n, FUN=function(ind) {
    ((x[[ind]]) - limits[1,ind])/(limits[2,ind]-limits[1,ind])
    })
  out <- rast(list1)
  return(out)
}

centroidGrid.lowest.grad <- function(x,M){
  i <- centroidGrid[x,1] ; j <- centroidGrid[x,2]
  searchRad <- floor(S/2)
  xr <- seq(i-searchRad,i+searchRad,1) ; xr <- xr[which(xr > 0 & xr <= nrowsInv)]
  xc <- seq(j-searchRad,j+searchRad,1) ; xc <- xc[which(xc > 0 & xc <= ncolsInv)]
  if (gradIn[i,j] > min(gradIn[xr,xc])) {
  return(as.numeric(expand.grid(xr,xc)[which.min(gradIn[xr,xc]),]))
  } else
  return(centroidGrid[x,])
}

score.funct <- function(N,arr,r1,c1,r2,c2) {((arr[r1,c1,N] - arr[r2,c2,N])^2)}

  single.func <- function(x,Mat) {
    r <- coords[x,1] ; c <- coords[x,2]
    xr <- seq(r-1,r+1,1) ; xr <- xr[which(xr > 0 & xr <= nrowsInv)]
    xc <- seq(c-1,c+1,1) ; xc <- xc[which(xc > 0 & xc <= ncolsInv)]
    diff <- Mat[xr,xc] - Mat[r,c]
    test <- length(which(diff==0)) > 1
    if(test==FALSE) {
        aroundPix <- table(Mat[xr,xc])
        newVal <- which.max(aroundPix)
      return(as.numeric(names(newVal)))} else {return(l[r,c])}
  } 
