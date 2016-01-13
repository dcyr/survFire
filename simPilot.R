####################################################################
rm(list=ls())
####################################################################
# setwd("/media/dcyr/Windows7_OS/Travail/SCF/fcEstimationExp")
wwd <- paste(getwd(), Sys.Date(), sep="/")
dir.create(wwd)
setwd(wwd)
rm(wwd)

####################################################################
####################################################################
######
######      initial landscape
library(raster)
####################################################################
####################################################################
initialLandscape <- raster(extent(0, 100000, 0, 100000), res=1000)
initialLandscape[] <- 50

## convertion from pixel to hectares
scaleFactor <- prod(res(initialLandscape))/10000

####################################################################
####################################################################
######
######fire size distribution
library(MASS)
####################################################################
####################################################################
### empirical reference
fireObs <- read.csv("../data/fireObs.csv", header=TRUE)
fireSizeObs <- fireObs[, "SIZE"]
fireSizeObs <- fireSizeObs[order(fireSizeObs)]
### fire size distribution
fireSizeLogNormalFit <- fitdistr(fireSizeObs, "lognormal")
fireSizeExpFit <- fitdistr(fireSizeObs, "exponential")

####################################################################
####################################################################
######
######      sim Parameters +
# sourcing function
source("../scripts/simFunc.R")
####################################################################
####################################################################

nRep <- 10
require(doSNOW)
clusterN <- 2  ### choose number of nodes to add to cluster.

cl = makeCluster(clusterN)
registerDoSNOW(cl)

for (fc in c(1000, 500, 250, 125, 75)) {
    output <- foreach(i = 1:nRep) %dopar%  { # i <- 1
        output <- sim(initialLandscape, simDuration = 4 * fc,
                      fireCycle = fc, fireSizeFit = fireSizeLogNormalFit)
        return(output)
    }
    save(output, file = paste0("simOutput_", fc, ".RData"))
    rm(output)
}

stopCluster(cl)




####################################################################
####################################################################
######
######      post-processing simulations
# sourcing function
#source("../simFunc.R")
####################################################################
####################################################################
#samplingEffort <- 100
#binWidth <- c(5,25,50,75,100,125,150,175,200,250,300)
#censInt <-c(100,300)
#spinOffDuration <- 2 * max(fireCycle)
