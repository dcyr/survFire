####################################################################
####################################################################
######
######  Simulation of  dendroecological sampling and FC estimation
######
######  Indicates which simulation output files to process,
######  with which methods, sampling effort and number of trials.
######
######  Set censoring function inputs
######
######  Dominic Cyr
####################################################################
####################################################################
rm(list=ls())
####################################################################
####################################################################
setwd("/media/dcyr/Windows7_OS/Travail/SCF/fcEstimationExp")
outputFolder <- paste(getwd(), "compiledOutputs", sep="/")
wwd <- paste(getwd(), Sys.Date(), sep="/")
dir.create(wwd)
setwd(wwd)
rm(wwd)

############################################
##### simulating field sampling and FC estimation
require(survival)
require(dplyr)
require(reshape2)
require(doSNOW)
require(parallel)
source("../scripts/censFnc.R")
source("../scripts/fcEstSurvFnc.R")
############################################
tsfFinal <- get(load(paste(outputFolder, "tsfFinal.RData", sep="/")))
###
sampleSize <- c(25, 50, 75, 94, 150, 250, 500)
nTrials <- 100
###
clusterN <- detectCores()-1  ### choose number of nodes to add to cluster.
sysName <- Sys.info()["sysname"]
opts <- list(chunkSize=clusterN)
#
if (sysName=="Windows") {
    cl = makeCluster(clusterN, rscript="Rscript.exe", type='SOCK')
}
if (sysName=="Linux") {
    cl = makeCluster(clusterN)
}

registerDoSNOW(cl)
survivalEstimates <- foreach(fc = unique(tsfFinal$fireCycle), .combine="rbind") %do% {
    # filtering by fireCycle
    trueTSF1 <- tsfFinal %>%
        filter(fireCycle == fc)
    foreach(treat = unique(trueTSF1$treatment), .combine="rbind") %do% {
        # filtering by treatment
        trueTSF2 <- trueTSF1 %>%
            filter(treatment == treat)
        foreach(r = unique(tsfFinal$replicate), .combine="rbind") %do% {#
            # filtering by replicate
            trueTSF3 <- trueTSF2 %>%
                filter(replicate == r)
            trueTSF <- trueTSF3$tsfFinal
            tmp <- foreach(ss = sampleSize, .combine="rbind",
                            .packages = c("reshape2", "doSNOW")) %do%  {
                tmp <- foreach(i = 1:nTrials, .combine="rbind",
                        .packages = c("survival")) %dopar%  {

                    # sampling true tsf
                    tsf <- sample(trueTSF, ss)
                    tsf[tsf == 0] <- 0.1
                    # transforming uncensored tsf sample into 'Surv' object
                    tsfUncensored <- Surv(tsf, rep(1, length(tsf)))
                    # applying censoring function
                    tsf <- censFnc(tsf, 100, 300)
                    ### estimating FC from censored samples
                    cox <- coxFitFnc(tsf)
                    weib <- weibFitFnc(tsf)
                    exp <- expFitFnc(tsf)
                    # estimating FC from uncensored samples
                    coxUncensored <- coxFitFnc(tsfUncensored)
                    weibUncensored <- weibFitFnc(tsfUncensored)
                    expUncensored <- expFitFnc(tsfUncensored)
                    #
                    tmp <- data.frame(cox, weib, exp, coxUncensored, weibUncensored, expUncensored)
                    return(round(tmp, 1))

                }
#                 #### testing methods
#                 ## gotta load survivalBootstrap first
#                 trueFC <- as.numeric(unique(survivalBootstrap %>%
#                     filter(fireCycle == fc,
#                            treatment == as.character(treat),
#                            replicate == r) %>%
#                     select( trueFC300)))
#
#                 plot(density(tmp$cox), col = "seagreen4")#, breaks = 20)
#                 lines(density(tmp$coxUncensored), col = "seagreen4", lty=3)#, breaks = 20)
#                 #lines(density(tmp$weib), col = "goldenrod2")
#                 #lines(density(tmp$exp), col = "indianred4")
#                 abline(v = trueFC)#, breaks = 20)
#                 abline(v = mean(tmp$cox), col = "seagreen4")#, breaks = 20)
#                 abline(v = mean(tmp$coxUncensored), col = "seagreen4", lty = 3)#, breaks = 20)
#
                tmp <- data.frame(fireCycle = as.numeric(fc),
                                    treatment = treat,
                                    replicate = as.numeric(r),
                                    sampleSize = ss,
                                   tmp)
                tmp <- melt(tmp, id.vars = c("fireCycle", "treatment", "replicate", "sampleSize"),
                     meas.vars = c("cox", "weib", "exp", "coxUncensored", "weibUncensored", "expUncensored"),
                     variable.name = "method",
                     value.name = "estimate")
                return(tmp)
            }
        }
    }
}
stopCluster(cl)
##
save(survivalEstimates, file = "survivalEstimatesFullDF.RData")

