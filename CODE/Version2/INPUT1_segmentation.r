########### 1. USER INPUT REQUIRED  - DEFINES THIN SECTION AND PARAMETERS FOR SEGMENTATION ########### 

### Choose the samples ###
Samples <- c("SK394A", "SK394C")

### Name of mineral phases being segmented ###
PhaseName <- "plag"

### Pixel size equivalent (in um) ###
pix_um <- 20

### Minimum crystal size to segment (in um) ###
minArea <- 180^2

### Choose which elements to perform the segmentation on. This could also be BSE, for example ###
Variables <- c("An") #This can be more than one variable

### Set maximum number of iterations the SLICAP algorithm is run for convergence ###
nIters <- 50

### Scale multiple parameters on relative scale (FALSE or TRUE) ###
rel.scale <- FALSE

### Set spacing of interval. Leave as NA to use automatic setting ###
S_manual <- NA

### Set minimum spacing interval ###
minS <- 5

### Set chemical weighting parameter. Leave as NA to use automatic setting ###
M_manual <- NA

### Set minimum chemical weighting parameter ###
minM <- 0

### Set value of q for AP algorithm ###
q_choice <- 0.001

