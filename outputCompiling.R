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


####################################################################
####################################################################
######
######      compiling simulation outputs
######
## Uncomment the following line if compiling 2000-yrs simulations outputs
outputFolder <- paste(outputFolder, "2000YrsSims", sep="/")
####################################################################
####################################################################

x <- list.files(outputFolder)
x <- x[grep(".RData", x)]
simInfo <- gsub(".RData", "", x)
simInfo <- strsplit(simInfo, "_")
fc <- as.numeric(lapply(simInfo, function(x) x[2]))
corr <- as.numeric(lapply(simInfo, function(x) x[3]))
###########################################
###########################################
outputList <- list()
require(raster)
for (i in seq_along(x)) {#
    ##
    output <- get(load(paste(outputFolder, x[i], sep="/")))
    fireCycle <- meanTSF <- year <- replicate <- areaBurned <- treatment <- list()
    ##
    for (r in seq_along(output)) {
        tsfStack <- output[[r]][["tsf"]]
        convFactor <- prod(res(tsfStack))/10000### to convert to hectares
        ## removing edge (10km)
        e <- extent(tsfStack, 11, nrow(tsfStack)-10, 11, ncol(tsfStack)-10)
        tsfStack <- crop(tsfStack, e)
        w <- tsfStack[[1]]
        w[] <- 1
        ## compiling 'true' FC statistics
        meanTSF[[r]] <- round(as.numeric(zonal(tsfStack, w, mean))[-1], 1)
        areaBurned[[r]] <- as.numeric(zonal(tsfStack==0, w, sum))[-1] * convFactor
        year[[r]] <- 1:nlayers(tsfStack)
        replicate[[r]] <- rep(r, nlayers(tsfStack))
        fireCycle[[r]] <- rep(fc[i], nlayers(tsfStack))
        treatment[[r]] <- rep(corr[i], nlayers(tsfStack))
        ## recording final TSF for downstream analyses
    }
    ##
    fireCycle <- as.numeric(do.call("cbind", fireCycle))
    treatment <- as.numeric(do.call("cbind", treatment))
    replicate <- as.numeric(do.call("cbind", replicate))
    year <- as.numeric(do.call("cbind", year))
    areaBurned_ha <- as.numeric(do.call("cbind", areaBurned))
    meanTSF <- as.numeric(do.call("cbind", meanTSF))
    ##
    outputList[[i]] <- data.frame(fireCycle, treatment, replicate, year, areaBurned_ha, meanTSF)
    print(i)
    rm(output)
}
outputCompiled <- do.call("rbind", outputList)
save(outputCompiled, file = "simOutputCompiled.RData")


# tsfFinalList <- list()
# require(raster)
# for (i in seq_along(x)) {#
#     ##
#     output <- get(load(paste(outputFolder, x[i], sep="/")))
#     fireCycle <- replicate <- treatment <- tsfFinal <- list()
#     ##
#     for (r in seq_along(output)) {
#
#         tsfFinalRaster <- output[[r]][["tsf"]]
#         tsfFinalRaster <- tsfFinalRaster[[nlayers(tsfFinalRaster)]]
#         convFactor <- prod(res(tsfFinalRaster))/10000### to convert to hectares
#         ## removing edge (10km)
#         e <- extent(tsfFinalRaster, 11, nrow(tsfFinalRaster)-10, 11, ncol(tsfFinalRaster)-10)
#         tsfFinalRaster <- crop(tsfFinalRaster, e)
#         nCells <- ncell(tsfFinalRaster)
#         ## compiling 'true' FC statistics
#         replicate[[r]] <- rep(r, nCells)
#         fireCycle[[r]] <- rep(fc[i], nCells)
#         treatment[[r]] <- rep(corr[i], nCells)
#         ## recording final TSF for downstream analyses
#         tsfFinal[[r]] <- values(tsfFinalRaster)
#     }
#     ##
#     fireCycle <- as.numeric(do.call("cbind", fireCycle))
#     treatment <- as.numeric(do.call("cbind", treatment))
#     replicate <- as.numeric(do.call("cbind", replicate))
#     tsfFinal <- as.numeric(do.call("cbind", tsfFinal))
#     ##
#     tsfFinalList[[i]] <- data.frame(fireCycle, treatment, replicate, tsfFinal)
#     print(i)
#     rm(output)
# }
# ##
# tsfFinal <-  do.call("rbind", tsfFinalList)
# save(tsfFinal, file = "tsfFinal.RData")

