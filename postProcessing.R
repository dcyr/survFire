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
#  refRaster <- get(load("../outputs/simOutput_125_0.RData"))
#  refRaster <- refRaster[[1]]$tsf[[1]]
totalArea <- 1000000
###



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
    filter(replicate == 1) %>%
    mutate(id = paste0(fireCycle, replicate, treatment))
output$treatment <- as.factor(output$treatment)

output <- output %>%
    arrange(year, fireCycle)

trueFCNames <- c("trueFC50", "trueFC150", "trueFC300")
output[, trueFCNames] <- NA

trueFC <- list()
for (fc in unique(output$fireCycle)) {
    trueFC[[fc]] <- output %>%
        filter(fireCycle == fc)
    year <- trueFC[[fc]]$year
    aab <- trueFC[[fc]]$areaBurned_ha


    trueFC[[fc]][,"trueFC50"] <- totalArea/rollmean(aab, 50, align = "right", fill = NA)
    trueFC[[fc]][,"trueFC150"] <- totalArea/rollmean(aab, 150, align = "right", fill = NA)
    trueFC[[fc]][,"trueFC300"] <- totalArea/rollmean(aab, 300, align = "right", fill = NA)


}

## unlisting
#trueFC <- do.call("rbind", trueFC)
## 'melting' dataframe
# require(reshape2)
# trueFC <- melt(trueFC, id.vars = c("fireCycle", "treatment", "replicate", "id", "year", "areaBurned_ha"),
#      meas.vars = c("meanTSF", "trueTSF50", "trueTSF150", "trueTSF300"))
#
# head(trueFC)
#
# unique(trueFC$id)
# ## recreating unique IDs
# trueFC <- trueFC %>%
#     mutate(id = as.factor(as.numeric(as.factor(paste(id, variable)))))
# summary(trueFC)
#
# ###### Plotting examples
# colors <- c("black", "darkred", "darkblue", "darkgreen")


#####################
#####################

############################################
##### graphiques - statistiques "r?elles"
############################################

yLimAAB <- c(0, max(as.numeric(lapply(trueFC, function(x) max(x$areaBurned_ha)))))

for (fc in c(62.5, 125, 250, 500, 1000)) {
    df <- trueFC[[fc]]
    x <- df$year
    meanTSF <- df$meanTSF
    FC50 <- df$trueFC50
    FC150 <- df$trueFC150
    FC300 <- df$trueFC300
    aab <- df$areaBurned_ha
    #
    realizedFC <- totalArea/mean(aab)


    png(filename = paste0("simTrueTSF", fc, ".png"),
        #width = 1800, height = (160*length(unique(df$species))+200), units = "px",
        width = 8, height = 5, units = "in", res = 300, pointsize = 8,
        bg = "white")

        #########
        par(mfrow=c(2,1), mar=c(0,4,3,1))
        options(scipen=7) ###force fixed notation (no scientific notation unless > 10E7 )
        #
        plot(x = x, y = aab, type="l", ylab="Annual area burned (ha)", xaxt="n", lwd=0.5, xlim=c(0,2000), ylim=yLimAAB)
        abline(v = seq(from=0, to=2000, by= 100), col = "gray", lty=3, lwd=.5)
        grid(col="gray", lty=3, lwd=.8)
        #
        text(max(x), y = yLimAAB[2]*.98, paste("Realized fire cycle (entire simulation) :", round(realizedFC, 1), "years"), adj = c(1,1), cex=1)
        text(max(x), y = yLimAAB[2]*.92, paste("Mean percent annual area burned :", round(100*1/realizedFC, 2), "%"), adj = c(1,1), cex=1)
        ###
        par(mar=c(5,4,0,1))
        ylim <- c(0,(realizedFC*3.5))
        #ylimits <- c(0, 300)
        plot(meanTSF, type="l", lwd=2, lty=1, col="black", ylab="Mean time since fire / True fire cycle", xlab="Simulated years", ylim=ylim)
        abline(v = seq(from=0, to=2000, by= 100), col = "gray", lty=3, lwd=.5)
        grid(col="gray", lty=3, lwd=.8)
        abline(h=realizedFC, lty=1, lwd=1, col="black")
        lines(FC50, lty=1, lwd=1, col="red")
        lines(FC150, lty=1, lwd=1, col="blue")
        lines(FC300, lty=1, lwd=1, col="darkgreen")
        legend(0, ylim[2],
               c("True mean TSF", "True FC (past 50 years)", "True FC (past 150 years)", "True FC (past 300 years)"),
               lty = 1, lwd = c(2,1,1,1), col = c("black", "red", "blue", "darkgreen"),
               cex=1,
               bg="white")

    dev.off()
}


#####################
#####################
#####################



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
