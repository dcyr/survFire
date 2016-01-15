####################################################################
####################################################################
####################################################################
rm(list=ls())
####################################################################
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
######
library(raster)
####################################################################
####################################################################
initialLandscape <- raster(extent(0, 120000, 0, 120000), res=1000)
initialLandscape[] <- NA

## convertion from pixel to hectares
scaleFactor <- prod(res(initialLandscape))/10000

####################################################################
####################################################################
######
######      fire size distribution
######
library(MASS)
####################################################################
####################################################################
### empirical reference
fireObs <- read.csv("../data/fireObs.csv", header=TRUE)
fireSizeObs <- fireObs[, "SIZE"]
fireSizeObs <- fireSizeObs[order(fireSizeObs)]
### fire size distribution
fireSizeFit <- fitdistr(fireSizeObs, "lognormal")
#fireSizeFit <- fitdistr(fireSizeObs, "exponential")


####################################################################
####################################################################
######
######      simulations
######
source("../scripts/simFunc.R")
####################################################################
####################################################################
nRep <- 10
require(doSNOW)
clusterN <- 3  ### choose number of nodes to add to cluster.

cl = makeCluster(clusterN)
registerDoSNOW(cl)

for (fc in c(1000, 500, 250, 125, 75)) {
    for (corr in c(-.3, 0, .3)) {
        ## create landscape at equilibrium according to simulated fire cycle
        initialLandscape[] <- round(rexp(ncell(initialLandscape), rate = 1/fc))
        initialLandscape[initialLandscape == 0] <- fc # removing zeros, remplace by fc
        ##
        output <- foreach(i = 1:nRep) %dopar%  {
            output <- sim(initialLandscape, simDuration = 300,
                          fireCycle = fc, fireSizeFit = fireSizeFit, corr = corr)
            return(output)
        }
        save(output, file = paste0("simOutput_", fc, "_", corr, ".RData"))
        rm(output)
    }

}

stopCluster(cl)
