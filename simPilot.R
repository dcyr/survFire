####################################################################
####################################################################
###### Simulation pilot
######
###### Prepare simulation inputs
###### Deploy experimental designS
######
###### Dominic Cyr
####################################################################
####################################################################
rm(list=ls())
####################################################################
####################################################################
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
initialLandscape <- raster(extent(0, 145000, 0, 145000), res=1000)
initialLandscape[] <- NA

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
source("../scripts/simFnc.R")
####################################################################
####################################################################
nRep <- 100
simDuration <- 300
fireCycles <- c(62.5, 125, 250, 500)
corrPAAB <- c(-0.5, 0, 0.5)
require(doSNOW)
require(parallel)

clusterN <-  detectCores() - 1  ### choose number of nodes to add to cluster.
#######
cl = makeCluster(clusterN)
registerDoSNOW(cl)
for (fc in fireCycles) {
    for (corr in corrPAAB) {
        #### initial landscape with tsf == 1
        # initialLandscape[] <- 1
        #### or initial landscape at equilibrium according to simulated fire cycle
        initialLandscape[] <- round(rexp(ncell(initialLandscape), rate = 1/fc))
        initialLandscape[initialLandscape == 0] <- fc # removing zeros, remplace by fc

        ##
        t1 <- Sys.time()
        output <- foreach(i = 1:nRep) %dopar%  {
            output <- sim(initialLandscape, simDuration = simDuration,
                          fireCycle = fc, fireSizeFit = fireSizeFit, corr = corr)
            return(output)
        }
        t2 <- Sys.time()
        print(paste(fc, corr, round(t2-t1, 1), "sec"))
        save(output, file = paste0("simOutput_", fc, "_", corr, ".RData"))
        rm(output)
    }
}
stopCluster(cl)
#######