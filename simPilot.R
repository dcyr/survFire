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
fireObs <- read.csv("../Data/fireObs.csv", header=TRUE)
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
source("../simFunc.R")
####################################################################
####################################################################

nRep <- 15
require(doSNOW)
clusterN <- 3  ### choose number of nodes to add to cluster.

cl = makeCluster(clusterN)
registerDoSNOW(cl)

for (fc in c(1000, 500, 250, 125, 75)) {
    output <- foreach(i = 1:nRep) %dopar%  { # i <- 1
        output <- sim(initialLandscape, simDuration = 5 * fc,
                      fireCycle = fc, fireSizeFit = fireSizeLogNormalFit)
        return(output)
    }
    save(output, file = paste0("simOutput_", fc, ".RData"))
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



#
#
# ####################################################################
# ####################################################################
#
# require(RColorBrewer)
# brks <- c(seq(from = 0, to = 200, by = 25), 1000)
# col <- brewer.pal(length(brks), "RdYlGn")
#
#
# hist(tsf)
# plot(tsf, main = "TSF")
#

# ###
# # ploting obs with fitted distrib
# ###
# require(EnvStats)
# fireCumsum <- cumsum(fireSize)
# fireCumProp <- fireCumsum/max(fireCumsum)
# fitCdfLogNorm <- plnorm(1:max(fireSize),
#                      meanlog = fireSizeLogNormalFit$estimate[1],
#                      sdlog = fireSizeLogNormalFit$estimate[2])
# fitPdfLogNorm <- dlnorm(1:max(fireSize),
#                         meanlog = fireSizeLogNormalFit$estimate[1],
#                         sdlog = fireSizeLogNormalFit$estimate[2])
# #fitExp <- pexp(1:max(fireSize), rate = fireSizeExpFit$estimate)
#
#
# options(scipen = 5)
#
# png(filename = "fireSizeDistrib.png",
#     width = 7, height = 4, res = 300, units = "in",
#     pointsize = 10)
#
#     par(mfrow=c(1,2))
#
#     ####### pdf plot
#     plot(density(fireSize, bw = "SJ"), ylim = c(0, max(fitPdfLogNorm)),
#          type = "l", lwd = 1, lty = 1, col = "black",
#          xlab = "Fire size (ha)",
#          ylab = "Probability density",
#          main = "Probability density function")
#     rug(fireSize, col='darkred')
#     ## add lognormal fit
#     lines(x = 1:length(fitPdfLogNorm), y = fitPdfLogNorm,
#           lty = 2, col = "black")
#     # cdf plot
#     plot(x = fireSize, y = (1:length(fireSize)/length(fireSize)),
#              type = "l", lwd = 1, lty = 1, col = "black",
#              xlab = "Fire size (ha)",
#              ylab = "Cumulative probability",
#              main = "Cumulative distribution function")
#     rug(fireSize, col='darkred') # add a rugplot of y3
#     ## add lognormal fit
#     lines(x = 1:length(fitCdfLogNorm), y = fitCdfLogNorm,
#           lty = 2, col = "black")
#
#
# dev.off()


# #################
# ###### lognorm.fit de la distribution des tailles de feux
# #################
#
# library(MASS)
# fire.size.lognorm.fit <- fitdistr(ind.fire.size.LTG, "lognormal")
#
# # ###### goodness of fit
# # shapiro.test(log(ind.fire.size.LTG))              ####fit pas bon, consid?rer ?chantillonner directement la table d'observations empiriques.
# # library(nortest)
# # ad.test(log(ind.fire.size.LTG))    ###### fit pas bon non plus avec Anderson-Darling
# ######
#
#
# ###### histogramme
# breaks <- seq(from=0, to=150000, by=10000)
# h <- hist(fireSizeDistEmp, breaks = break.seq)
# xfit <- seq(min(fireSize),max(ind.fire.size.LTG),length=40)
# yfit <- dlnorm(xfit, meanlog = fire.size.lognorm.fit$estimate[1], sdlog = fire.size.lognorm.fit$estimate[2], log = FALSE)
# yfit <- yfit*diff(h$mids[1:2])*length(ind.fire.size.LTG)
# lines(xfit, yfit, col="black", lwd=1)





head(fireObs)
fire.year <- seq(from=1959, to=1999)
annual.area.burned <- rep(0.0001, length(fire.year))
annual.area.burned[fire.year %in% ind.fire.year.LTG] <- by(ind.fire.size.LTG, ind.fire.year.LTG, sum)

fireSizeDist <- fitdistr(fireObs, "lognormal")
fireNum <-

tsf <- initTsf

# source(fireModelPath)

for (i in 1:simDuration) {
    ## first burn forest
    tsf <- fireModel(tsf)
    ## forest aging
    tsf <- tsf + 1

    ## sampling
    # tsfSample <- sample(tsf, sampSize)


    ## censoring
    # status <- NA #defining status
    # tsfSampleCens <- tsfSample[which ]

}

library(survival)



#########################################
######## linear censoring
#########################################
lin.cens.fn <- function(int) { # where 'int' interval within which censoring occur
    cens.prob <- rep(0, simulation.length)
    cens.prob[sim.year > min(cens.lin.lim) & sim.year <= max(cens.lin.lim)] <- 1/(max(cens.lin.lim)-min(cens.lin.lim))
    cumul.cens.prob <- rep(NA, simulation.length)
    for (i in 1:simulation.length) {
        cumul.cens.prob[i] <- sum(cens.prob[1:i])
    }
    status <- rep(NA, length(int))
    int.cens <- as.numeric(status)

    for (j in 1:max(int))      {
        status[int==j] <- ifelse(j<=min(cens.lin.lim), TRUE, ifelse(j>=max(cens.lin.lim), FALSE,
                                                                    sample(c(FALSE,TRUE), prob=c(cumul.cens.prob[j], (1-cumul.cens.prob[j]))))) #### probability of censoring (losing information about past fire, 0=> censored (info lost), 1 => not censored (info remains))
        int.cens[int==j] <- ifelse(status[int==j]==0 & int[int==j]>=max(cens.lin.lim), sample(seq(min(cens.lin.lim), max(cens.lin.lim)), 1),
                                   ifelse(status[int==j]==0 & int[int==j]<max(cens.lin.lim), sample(seq(min(cens.lin.lim), j), 1), j))  ###atttribution d'un nouveau tsf lorsque censur?
    }
    list(int.cens=int.cens, status=status)
}
