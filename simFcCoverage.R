####################################################################
####################################################################
####################################################################
rm(list=ls())
####################################################################
###################################################################
# setwd("/media/dcyr/Windows7_OS/Travail/SCF/fcEstimationExp")
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
require(boot)
source("../scripts/censFnc.R")
source("../scripts/fcEstSurvFnc.R")
############################################
tsfFinal <- get(load(paste(outputFolder, "tsfFinal.RData", sep="/")))


# the following design took XhXXmin to run with 3 cores on my machine
sampleSize <- c(25, 50, 75, 94, 150, 250, 500)
replicates <- unique(tsfFinal$replicate)
bootMethod <- c("norm", "perc","bca", "basic")
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
###
survivalBootstrap <- foreach(i = unique(tsfSample$ID), .combine="rbind",
                             .packages = c("survival", "boot", "dplyr", "foreach")) %dopar% {
    ###
    df <- tsfSample %>%
        filter(ID == i)
    fc <- unique(df$fireCycle)
    treat <- unique(df$treatment)
    r <- unique(df$replicate)
    #ssh -i "dcyr-AWS.pem" root@52.87.168.224

    tmp <- foreach(ss = sampleSize, .combine="rbind") %do% {
        # sampling true tsf
        tsf <- sample_n(df, ss)
        tsf <- tsf$tsfFinal
        tsf[tsf == 0] <- 0.1
        # applying censoring function
        tsf <- censFnc(tsf, 100, 300)
        ##############################################
        ### bootstrop estimation of FC from censored samples
        coxBoot <- boot(as.matrix(tsf), statistic = coxFitFnc, R = nBootstrap, sim = "ordinary")
        weibBoot <- boot(as.matrix(tsf), statistic = weibFitFnc, R = nBootstrap, sim = "ordinary")
        expBoot <- boot(as.matrix(tsf), statistic = expFitFnc, R = nBootstrap, sim = "ordinary")

        tmp <- foreach(m = bootMethod, .combine="rbind") %do% {
            ### computing 95%CI (sometimes fails)
            try(coxCI <- boot.ci(coxBoot, type = m))
            try(weibCI <- boot.ci(weibBoot, type = m))
            try(expCI <- boot.ci(expBoot, type = m))

            ###
            cox <- data.frame(method = "cox", estimate = round(coxBoot$t0,1), ll = NA, ul = NA)
            weib <- data.frame(method = "weib", estimate = round(weibBoot$t0,1), ll = NA, ul = NA)
            exp <- data.frame(method = "exp", estimate = round(expBoot$t0,1), ll = NA, ul = NA)

            try(cox[,c("ll", "ul")] <- round(t(sapply(coxCI[-(1:3)],function(x) tail(c(x),2))),1))
            try(weib[,c("ll", "ul")] <- round(t(sapply(weibCI[-(1:3)],function(x) tail(c(x),2))),1))
            try(exp[,c("ll", "ul")] <- round(t(sapply(expCI[-(1:3)],function(x) tail(c(x),2))),1))

            tmp <- rbind(cox, weib, exp)

            tmp <- data.frame(fireCycle = as.numeric(fc),
                              treatment = treat,
                              replicate = as.numeric(r),
                              bootMethod = m,
                              sampleSize = ss,
                              tmp)
            return(tmp)

        }
        return(tmp)
    }
    return(tmp)
}
stopCluster(cl)
## a little cleaning up
rownames(survivalBootstrap) <- 1:nrow(survivalBootstrap)
##
save(survivalBootstrap, file = "survivalBootstrapFullDF.RData")
