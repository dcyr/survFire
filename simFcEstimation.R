####################################################################
####################################################################
####################################################################
rm(list=ls())
####################################################################
####################################################################
#setwd("/media/dcyr/Windows7_OS/Travail/SCF/fcEstimationExp")
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
# the following design took 2h40min to run with 3 cores on my machine
sampleSize <- c(10, 25, 50, 75, 94, 150, 250, 500)
nTrials <- 400
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
t1 <- Sys.time()
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
                    cox <- coxFitFnc(tsf)$cycle
                    weib <- weibFitFnc(tsf)$cycle
                    exp <- expFitFnc(tsf)$cycle
                    # estimating FC from uncensored samples
                    coxUncensored <- coxFitFnc(tsfUncensored)$cycle
                    weibUncensored <- weibFitFnc(tsfUncensored)$cycle
                    expUncensored <- expFitFnc(tsfUncensored)$cycle
                    #
                    tmp <- data.frame(cox, weib, exp, coxUncensored, weibUncensored, expUncensored)
                    return(round(tmp, 1))

                }
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
t2 <- Sys.time()
##
print(t2-t1)
save(survivalEstimates, file = "survivalEstimates.RData")

