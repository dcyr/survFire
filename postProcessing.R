####################################################################
####################################################################
####################################################################
rm(list=ls())
require(raster)
####################################################################
####################################################################
setwd("/media/dcyr/Windows7_OS/Travail/SCF/fcEstimationExp")
outputFolder <- paste(getwd(), "compiledOutputs", sep="/")
wwd <- paste(getwd(), Sys.Date(), sep="/")
dir.create(wwd)
setwd(wwd)
rm(wwd)

### fetching total area
refRaster <- get(load("../outputs/simOutput_75_0.RData"))
refRaster <- refRaster[[1]]$tsf[[1]]
totalArea <- ncell(refRaster) * prod(res(refRaster))/10000
###
rm(refRaster)


####################################################################
####################################################################
######
######      2000-yrs constant fc simulations
######
require(dplyr)
require(ggplot2)
require(zoo) ### for 'rollmean'
####################################################################
####################################################################
output <- get(load(paste(outputFolder, "simOutputCompiled2000.RData", sep="/")))

output <- output %>%
    mutate(id = paste0(fireCycle, replicate, treatment))
output$treatment <- as.factor(output$treatment)

uniqueSimID <- unique(output$id)
output$simID <- match(output$id, sample(uniqueSimID))

head(output)
constantFC <- output %>%
    arrange(year, fireCycle)

for (fc in unique(constantFC$fireCycle)) {
    trueTSF <- constantFC %>%
        filter(fireCycle == fc)
    aab <- aab$areaBurned_ha


    trueTSF <- matrix(NA, ncol = 3, nrow = length(aab))
    trueTSF[,1] <- totalArea/rollmean(aab, 50, align = "right", fill = NA)
    trueTSF[,2] <- totalArea/rollmean(aab, 150, align = "right", fill = NA)
    trueTSF[,3] <- totalArea/rollmean(aab, 300, align = "right", fill = NA)
    trueTSFNames <- c("trueTSF50", "trueTSF150", "trueTSF300")
    colnames(trueTSF) <- trueTSFNames


    constantFC[match(year, constantFC$year),] <-trueTSF
}
head(constantFC)




trueTSF <- ggplot(constantFC, aes(x=year, y=meanTSF,
                                  group=id, color = treatment)) +
    #scale_colour_manual(values = colTreatment) +
    #stat_summary(fun.data = "median_hilow", geom = "smooth", size=5) +
    geom_line(aes(group=id),
              size=1,
              type=1,
              alpha=0.5) +
    facet_wrap( ~ fireCycle, ncol = 1,
                scales = "free")

trueTSF

####################################################################
####################################################################
######
######      sampling landscape
######
source("../scripts/censFnc.R")
source("../scripts/fcEstSurvFnc.R")
require(dplyr)
require(ggplot2)
colTreatment <- c("dodgerblue2", "black", "red3")
####################################################################
####################################################################

### a little tidying up
tsfFinal <- get(load(paste(outputFolder, "tsfFinal.RData", sep="/")))
rm(tsfFinalList)
output <- get(load(paste(outputFolder, "simOutputCompiled.RData", sep="/")))

require(dplyr)
output <- output %>%
    mutate(id = paste0(fireCycle, replicate, treatment))
output$treatment <- as.factor(output$treatment)

df0 <- output %>% distinct(fireCycle, treatment, replicate) %>%
    mutate(meanTSF = fireCycle,
           areaBurned_ha = NA,
           year = 0)

output <-  rbind(output, df0)



head(output)
uniqueSimID <- unique(output$id)
output$simID <- match(output$id, sample(uniqueSimID))



####################################################################
####################################################################
######
######      simulation example (one replicate)
######
require(zoo) ## for rollmean
####################################################################
####################################################################

constantFC <- output %>%
    # select only constant fire cycle, one replicate,
    filter(replicate == 1,
           treatment == "0") %>%
    arrange(year, fireCycle)

for (fc in unique(constantFC$fireCycle)) {
   aab <- constantFC %>%
       filter(fireCycle == fc)
   aab <- aab$areaBurned_ha

    trueFC <- rollmean()
}




head(constantFC)


meanTSF <- ggplot(constantFC, aes(x=year, y=meanTSF,
                              group=id, color = treatment)) +
    #scale_colour_manual(values = colTreatment) +
    #stat_summary(fun.data = "median_hilow", geom = "smooth", size=5) +
    geom_line(aes(group=id),
              size=1,
              type=1,
              alpha=0.5) +
    facet_wrap( ~ fireCycle, ncol = 1,
               scales = "free")

meanTSF


meanTSF <- ggplot(output, aes(x=year, y=meanTSF,
                    group=id, color = treatment)) +
    scale_colour_manual(values = colTreatment) +
    #stat_summary(fun.data = "median_hilow", geom = "smooth", size=5) +
    geom_line(aes(group=id),
              size=1,
              type=1,
               alpha=0.5) +
    facet_wrap(fireCycle ~ treatment, ncol = 3,
               scales = "free")


png(filename = "simTSF.png", width = 800, height = 1024, units = "px",
    pointsize = 32)

    print(meanTSF)
dev.off()


aab <- ggplot(output, aes(x=year, y=areaBurned_ha,
                              group=id, color = treatment)) +
    scale_colour_manual(values = colTreatment) +
    #stat_summary(fun.data = "median_hilow", geom = "smooth", size=5) +
    geom_line(aes(group=id),
              size=1,
              type=1,
              alpha=0.5) +
    facet_wrap(fireCycle ~ treatment, ncol = 3,
                scales = "fixed")


png(filename = "simAAB.png", width = 800, height = 1024, units = "px",
    pointsize = 32)

print(aab)
dev.off()
