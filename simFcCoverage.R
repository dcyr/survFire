####################################################################
####################################################################
####################################################################
rm(list=ls())
####################################################################
####################################################################
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
bootMethod <- c("basic", "bca")
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
survivalBootstrap <- foreach(i = unique(tsfSample$ID), .combine="rbind",
                             .packages = c("survival", "boot", "dplyr", "foreach")) %dopar% {
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
        # applying censoring function
        tsf <- censFnc(tsf, 100, 300)
        tmp <- foreach(m = bootMethod, .combine="rbind") %do% {
            ##############################################
            ### bootstrop estimation of FC from censored samples
            coxBoot <- boot(as.matrix(tsf), statistic = coxFitFnc, R = nBootstrap, sim = "ordinary")
            weibBoot <- boot(as.matrix(tsf), statistic = weibFitFnc, R = nBootstrap, sim = "ordinary")
            expBoot <- boot(as.matrix(tsf), statistic = expFitFnc, R = nBootstrap, sim = "ordinary")


            ### computing 95%CI (sometimes fails)
            try(coxBoot <- boot.ci(coxBoot, type = m))
            try(weibBoot <- boot.ci(weibBoot, type = m))
            try(expBoot <- boot.ci(expBoot, type = m))

            ###
            cox <- c(method = "cox", estimate = NA, ll = NA, ul = NA)
            weib <- c(method = "weib", estimate = NA, ll = NA, ul = NA)
            exp <- c(method = "exp", estimate = NA, ll = NA, ul = NA)

            try(cox <- data.frame(method = "cox", estimate = coxBoot$t0,
                              ll = ifelse(class(coxBoot) =="bootci", coxBoot[[m]][4], NA),
                              ul = ifelse(class(coxBoot) =="bootci", coxBoot[[m]][5], NA)))
            try(weib <- data.frame(method = "weib", estimate = weibBoot$t0,
                               ll = ifelse(class(weibBoot) =="bootci", weibBoot[[m]][4], NA),
                               ul = ifelse(class(weibBoot) =="bootci", weibBoot[[m]][5], NA)))
            try(exp <- data.frame(method = "exp", estimate = expBoot$t0,
                              ll = ifelse(class(expBoot) =="bootci", expBoot[[m]][4], NA),
                              ul = ifelse(class(expBoot) =="bootci", expBoot[[m]][5], NA)))

            tmp <- rbind(cox, weib, exp)

            tmp <-data.frame(fireCycle = as.numeric(fc),
                             treatment = treat,
                             replicate = as.numeric(r),
                             bootMethod = m,
                             sampleSize = ss,
                             round(tmp, 1))
            return(tmp)

        }
        return(tmp)
    }
    return(tmp)
}
t2 <- Sys.time()
stopCluster(cl)
rownames(survivalBootstrap) <- 1:nrow(survivalBootstrap)
##
print(t2-t1)
save(survivalBootstrap, file = "survivalBootstrapFullDF.RData")
