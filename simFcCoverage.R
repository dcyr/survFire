####################################################################
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


# the following design took XhXXmin to run with 3 cores on my machine
sampleSize <- c(25, 50, 75, 94, 150, 250, 500)
replicates <- unique(tsfFinal$replicate)
resamplingEffort <- c(0.5, 1, 2)
nBootstrap <- 1000

## shrinking table to a collection of sample of maximum sample size
tsfSample <- tsfFinal %>%
    filter(replicate %in% replicates)
## creating ID variable
tsfSample <- mutate(tsfSample, ID = as.numeric(as.factor(paste(fireCycle, treatment, replicate))))

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

###
registerDoSNOW(cl)
t1 <- Sys.time()
###
survivalBootstrap <- foreach(i = unique(tsfSample$ID), .combine="rbind") %do% {
    ###
    df <- tsfSample %>%
        filter(ID == i)
    fc <- unique(df$fireCycle)
    treat <- unique(df$treatment)
    r <- unique(df$replicate)
    #
    tmp <- foreach(ss = sampleSize, .combine="rbind") %do% {
        # sampling true tsf
        tsf <- sample_n(df, ss)
        tsf <- tsf$tsfFinal
        tsf[tsf == 0] <- 0.1
        # transforming uncensored tsf sample into 'Surv' object
        tsfUncensored <- Surv(tsf, rep(1, length(tsf)))
        # applying censoring function
        tsf <- censFnc(tsf, 100, 300)
        #
        tmp <- foreach(rs = resamplingEffort, .combine="rbind") %do% {#

            tmp <- foreach(r = 1:nBootstrap, .combine="rbind",
                           .packages = c("survival")) %dopar% {#
                ### Bootstrap sampling
                tsfIndex <- 1:nrow(tsf)
                tsfIndex <- sample(tsfIndex, rs*max(tsfIndex), replace = TRUE)
                tsfBoot <- tsf[tsfIndex]
                tsfBootUncensored <- tsfUncensored[tsfIndex]
                ### estimating FC from censored samples
                cox <- coxFitFnc(tsfBoot)$cycle
                weib <- weibFitFnc(tsfBoot)$cycle
                exp <- expFitFnc(tsfBoot)$cycle
                ### estimating FC from uncensored samples
                coxUncensored <- coxFitFnc(tsfBootUncensored)$cycle
                weibUncensored <- weibFitFnc(tsfBootUncensored)$cycle
                expUncensored <- expFitFnc(tsfBootUncensored)$cycle
                tmp <- data.frame(cox, weib, exp, coxUncensored, weibUncensored, expUncensored)
                return(round(tmp, 1))
            }

            tmp <- data.frame(fireCycle = as.numeric(fc),
                              treatment = treat,
                              replicate = as.numeric(r),
                              sampleSize = ss,
                              resamplingEffort = rs,
                              tmp)
            tmp <- melt(tmp, id.vars = c("fireCycle", "treatment", "replicate", "sampleSize", "resamplingEffort"),
                        meas.vars = c("cox", "weib", "exp", "coxUncensored", "weibUncensored", "expUncensored"),
                        variable.name = "method",
                        value.name = "estimate")
            return(tmp)
        }
        return(tmp)
    }
    return(tmp)
}



stopCluster(cl)
t2 <- Sys.time()
##
print(t2-t1)
save(survivalBootstrap, file = "survivalBootstrap.RData")

