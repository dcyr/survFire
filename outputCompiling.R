####################################################################
####################################################################
####################################################################
rm(list=ls())
####################################################################
####################################################################
setwd("/media/dcyr/Windows7_OS/Travail/SCF/fcEstimationExp")
outputFolder <- paste(getwd(), "outputs", sep="/")
wwd <- paste(getwd(), Sys.Date(), sep="/")
dir.create(wwd)
setwd(wwd)
rm(wwd)

####################################################################
####################################################################
######
######      compiling simulation outputs
######      +
######      simulated sampling
# source("../scripts/censFnc.R")
# samplingEffort <- c(10, 20, 30, 50, 75, 94, 250)
# resamplingN <- 10000
# binWidth <- c(5, 25, 50, 75, 100, 125, 150, 175, 200, 250, 300)
# censMin <- 100
# censMax <- 300

####################################################################
####################################################################

x <- list.files(outputFolder)
simInfo <- gsub(".RData", "", x)
simInfo <- strsplit(simInfo, "_")
fc <- as.numeric(lapply(simInfo, function(x) x[2]))
corr <- as.numeric(lapply(simInfo, function(x) x[3]))

require(doSNOW)
clusterN <- 2  ### choose number of nodes to add to cluster.

cl = makeCluster(clusterN)
registerDoSNOW(cl)

###########################################
## realization
fireCycle <- trueFC <- meanTSF <- year <- replicate <- areaBurned <- treatment <- list()
## simulated sampling
###########################################
outputList <- foreach(i = seq_along(x)) %dopar%  {#
    require(raster)
    ##
    output <- get(load(paste(outputFolder, x[i], sep="/")))
    fireCycle <- meanTSF <- year <- replicate <- areaBurned <- treatment <-  list()
    ##
    for (r in seq_along(output)) {

        tsfStack <- output[[r]][["tsf"]]
        convFactor <- prod(res(tsfStack))/10000### to convert to hectares
        ## remove edge effet
        e <- extent(tsfStack, 11, nrow(tsfStack)-10, 11, ncol(tsfStack)-10)
        tsfStack <- crop(tsfStack, e)
        w <- tsfStack[[1]]
        w[] <- 1
        ## compiling 'true' statistics
        meanTSF[[r]] <- round(as.numeric(zonal(tsfStack, w, mean))[-1], 1)
        areaBurned[[r]] <- as.numeric(zonal(tsfStack==0, w, sum))[-1] * convFactor
        year[[r]] <- 1:nlayers(tsfStack)
        replicate[[r]] <- rep(r, nlayers(tsfStack))
        fireCycle[[r]] <- rep(fc[i], nlayers(tsfStack))
        treatment[[r]] <- rep(corr[i], nlayers(tsfStack))
        ##

        ##
        rm(tsfStack)
        print(r)
    }
    ##
    fireCycle <- as.numeric(do.call("cbind", fireCycle))
    treatment <- as.numeric(do.call("cbind", treatment))
    replicate <- as.numeric(do.call("cbind", replicate))
    year <- as.numeric(do.call("cbind", year))
    areaBurned_ha <- as.numeric(do.call("cbind", areaBurned))
    meanTSF <- as.numeric(do.call("cbind", meanTSF))

    ##
    output <- data.frame(fireCycle, treatment, replicate, year, areaBurned_ha, meanTSF)
    return(output)
    #print(paste(i, "of", length(x)))
}
stopCluster(cl)
output <- do.call("rbind", outputList)
save(output, file = "simOutputCompiled.RData")


## simulated sampling
###########################################
cl = makeCluster(clusterN)
registerDoSNOW(cl)

tsfFinalList <- foreach(i = seq_along(x)) %dopar%  {#
    require(raster)
    ##
    output <- get(load(paste(outputFolder, x[i], sep="/")))
    tsfFinal <- NULL
    ##
    for (r in seq_along(output)) {

        tsfStack <- output[[r]][["tsf"]]
        tsfFinal <- cbind(tsfFinal, values(tsfStack[[nlayers(tsfStack)]]))
        rm(tsfStack)
    }
    return(tsfFinal)
    #print(paste(i, "of", length(x)))
}


stopCluster(cl)
names(tsfFinalList) <- paste0(gsub("simOutput_|.RData", "", x))
save(tsfFinalList, file = "tsfFinal.RData")

