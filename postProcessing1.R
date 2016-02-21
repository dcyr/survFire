################################################################################
################################################################################
######
######      Memory-demanding steps to do once with clean workspace
######
################################################################################
################################################################################

################################################################################
rm(list=ls())
require(dplyr)
require(raster)
require(reshape2)
require(dplyr)
require(zoo)
################################################################################
################################################################################
# setwd("/media/dcyr/Windows7_OS/Travail/SCF/fcEstimationExp")
outputFolder <- paste(getwd(), "compiledOutputs", sep="/")
wwd <- paste(getwd(), Sys.Date(), sep="/")
dir.create(wwd)
setwd(wwd)
rm(wwd)

# ### fetching total area from an output raster stack
# refRaster <- get(load("../outputs/simOutput_125_0.RData"))
# refRaster <- refRaster[[1]]$tsf[[1]]
# ### removing 10-km edge (should check if that's coherent with upstream output compiling)
# e <- extent(refRaster, 11, nrow(refRaster)-10, 11, ncol(refRaster)-10)
# refRaster <- crop(refRaster, e)
totalArea <- 1562500
# totalArea <- ncell(refRaster)
# totalArea <- ncell(refRaster) * prod(res(refRaster))/10000 ### total area in hectare




####################################################################
######  true fire cycle statistics computation
output <- get(load(paste(outputFolder, "simOutputCompiled.RData", sep="/")))
####################################################################
output <- output %>%
    mutate(id = paste0(fireCycle, replicate, treatment))
output$treatment <- as.factor(output$treatment)
### ordering (following loop depends or ordered time series)
output <- output %>%
    arrange(year, replicate, treatment, fireCycle)

trueFCNames <- c("trueFC50", "trueFC150", "trueFC300")
output[, trueFCNames] <- NA

trueFC <- list()
for (fc in unique(output$fireCycle)) {
    trueFC[[fc]] <- list()
    for (treat in unique(output$treatment)) {
        trueFC[[fc]][[treat]] <- list()
        for (r in unique(output$replicate)) {
            trueFC[[fc]][[treat]][[r]] <- output %>%
                filter(fireCycle == fc,
                       treatment == treat,
                       replicate == r)
            year <- trueFC[[fc]][[treat]][[r]]$year
            aab <-trueFC[[fc]][[treat]][[r]]$areaBurned_ha


            trueFC[[fc]][[treat]][[r]][,"trueFC50"] <- totalArea/rollmean(aab, 50, align = "right", fill = NA)
            trueFC[[fc]][[treat]][[r]][,"trueFC150"] <- totalArea/rollmean(aab, 150, align = "right", fill = NA)
            trueFC[[fc]][[treat]][[r]][,"trueFC300"] <- totalArea/rollmean(aab, 300, align = "right", fill = NA)
        }
    }

}
### Unlisting everything (a "foreach" construct in the previous loop could avoid that)
tmp <- list()
for (fc in unique(output$fireCycle)) {
    tmp2 <- list()
    for (treat in unique(output$treatment)) {
        tmp2[[treat]] <- do.call("rbind", trueFC[[fc]][[treat]])
    }
    tmp[[fc]] <- do.call("rbind", tmp2)
}
trueFC <-  do.call("rbind", tmp)

trueFC <- merge(trueFC, trueFC %>%
                    group_by(fireCycle, treatment, replicate) %>%
                    summarise(meanAAB = mean(areaBurned_ha)))

### Keep only complete cases (year == 300)
trueFC <- trueFC[complete.cases(trueFC),]
### further house cleaning
rownames(trueFC) <-  1:nrow(trueFC)
trueFC <- trueFC[,-which(colnames(trueFC) %in% c("id", "areaBurned_ha"))]
trueFC <- trueFC %>%
    arrange(fireCycle, treatment, replicate)
###
rm(output)



####################################################################
######  bootstrap estimates
survivalBootstrap <- get(load(paste(outputFolder, "survivalBootstrap.RData", sep="/")))
################################################################################
################################################################################


################################################################################
#### Filtering results that were deemed uninteresting through a trial and error process
##### (Usually because it creates CI so wide they were useless)
survivalBootstrap <- survivalBootstrap %>%
    #filter(propNonFinite > 0.995) %>%
    filter(is.finite(estimate)) %>%
    filter((fireCycle >= 1000) == F) %>%
    filter((fireCycle >= 125 & sampleSize <= 10) == F) %>%
    filter((fireCycle >= 500 & method == "weib") == F) %>%
    filter((fireCycle >= 250 & sampleSize <= 25 & method == "weib") == F) %>%
    filter((fireCycle >= 500 & sampleSize <= 25 &  treatment == "-0.5") == F)
#     filter((sampleSize == 25 & fireCycle >= 1000) == F) %>%
#     filter((method == "exp" & fireCycle >= 1000 & treatment == "-0.5") == F) %>%
#     filter((sampleSize <= 94 & fireCycle >= 1000 & treatment == "-0.5") == F)

### the following takes a while ...
survivalBootstrap <- merge(survivalBootstrap, trueFC)
save(survivalBootstrap, file = "survivalBootstrap.RData")

################################################################################
##### loading FC estimation obtained from simulated field sampling experiment
survivalEstimates <- get(load(paste(outputFolder, "survivalEstimatesFullDF.RData", sep ="/")))
survivalEstimates$sampleSize <- as.numeric(survivalEstimates$sampleSize)
################################################################################
################################################################################


##### Filtering results that were deemed uninteresting through a trial and error
##### process (Usually because it creates CI so wide they're useless)
survivalEstimates <- survivalEstimates %>%
    #filter(propNonFinite > 0.995) %>%
    filter(is.finite(estimate)) %>%
    filter((fireCycle >= 1000) == F) %>%
    filter((fireCycle >= 125 & sampleSize <= 10) == F) %>%
    filter((fireCycle >= 500 & method == "weib") == F) %>%
    filter((fireCycle >= 250 & sampleSize <= 25 & method == "weib") == F) %>%
    filter((fireCycle >= 500 & sampleSize <= 25 &  treatment == "-0.5") == F)
#     filter((sampleSize == 25 & fireCycle >= 1000) == F) %>%
#     filter((method == "exp" & fireCycle >= 1000 & treatment == "-0.5") == F) %>%
#     filter((sampleSize <= 94 & fireCycle >= 1000 & treatment == "-0.5") == F)

### the following takes a while ...
survivalEstimates <- merge(survivalEstimates, trueFC)
save(survivalEstimates, file = "survivalEstimates.RData")

