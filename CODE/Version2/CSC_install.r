df <- installed.packages() ; df <- as.data.frame(df)

requiredPackages <- c("flexmix","mvtnorm","terra","apcluster")
packageCheck <- match(requiredPackages,df$Package)
packageCheckID <- which(is.na(packageCheck)==TRUE)

if(length(packageCheckID)>0){
    install.packages(requiredPackages[packageCheckID],repos="https://cloud.r-project.org")
}

dir.create("../Phase-classification")
dir.create("../Phase-classification/Samples")
dir.create("../Crystal-segmentation")
dir.create("../Zone-correlation")
