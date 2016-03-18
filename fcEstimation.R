####################################################################
####################################################################
####################################################################
rm(list=ls())
####################################################################
###################################################################
setwd("/media/dcyr/Windows7_OS/Travail/SCF/fcEstimationExp")
wwd <- paste(getwd(), Sys.Date(), sep = "/")
dir.create(wwd)
setwd(wwd)
rm(wwd)
###
tsfTable <- read.csv("../data/tsf.csv")
### create 'surv' object
tsfTable <- tsfTable[complete.cases(tsfTable),]

##################################
######## bootstrapping
require(survival)
require(boot)
source("../scripts/fcEstSurvFnc.R")
##################################
tsf <- Surv(tsfTable$samplingDate - tsfTable$fireDate, tsfTable$status)
weight <- tsfTable$weight

bootMethod <- c("basic", "bca", "norm", "perc")
nBootstrap <- 25000

coxBoot <- boot(as.matrix(tsf), statistic = coxFitFnc, R = nBootstrap, weights = weight)
weibBoot <- boot(as.matrix(tsf), statistic = weibFitFnc, R = nBootstrap, weights = weight)
expBoot <- boot(as.matrix(tsf), statistic = expFitFnc, R = nBootstrap, weights = weight)

for (m in bootMethod) {
    coxCI <- boot.ci(coxBoot, type = m)
    weibCI <- boot.ci(weibBoot, type = m)
    expCI <- boot.ci(expBoot, type = m)

    cox <- data.frame(method = "Cox", bootMethod = m, estimate = round(coxBoot$t0,1), ll = NA, ul = NA)
    weib <- data.frame(method = "Weibull", bootMethod = m, estimate = round(weibBoot$t0,1), ll = NA, ul = NA)
    exp <- data.frame(method = "Exponential", bootMethod = m, estimate = round(expBoot$t0,1), ll = NA, ul = NA)

    cox[,c("ll", "ul")] <- round(t(sapply(coxCI[-(1:3)],function(x) tail(c(x),2))),1)
    weib[,c("ll", "ul")] <- round(t(sapply(weibCI[-(1:3)],function(x) tail(c(x),2))),1)
    exp[,c("ll", "ul")] <- round(t(sapply(expCI[-(1:3)],function(x) tail(c(x),2))),1)

    if (m == bootMethod[1]) {
        fcEst <- rbind(cox, weib, exp)
    } else {
        fcEst <- rbind(fcEst, cox, weib, exp)
    }
}
rownames(fcEst) <- 1:nrow(fcEst)
save(fcEst, file = "fcEst.RData")
##################################


##################################
get(load("../compiledOutputs/fcEst.RData"))
##################################
fcEst$bootMethod <- factor(fcEst$bootMethod, levels = c("basic", "bca", "norm", "perc"))
require(ggplot2)
png(filename = paste0("fcEstimation.png"),
    width = 6, height = 4, units = "in", res = 600, pointsize = 10,
    bg = "white")

    ggplot(fcEst, aes(method, estimate, bootMethod = method, color=method, linetype = bootMethod )) +
        #geom_point(size=3, position = "dodge") + # point for group mean
        geom_pointrange(aes(ymax=ul, ymin=ll), width=0.8, size = 0.2,
                        position= position_dodge(width = 0.3)) +
        scale_linetype_manual(values = 1:4, guide = guide_legend(title = "Type of CI", reverse = T)) +
        # geom_errorbar(aes(ymax=ul, ymin=ll), width=0.2, position= position_dodge(width = 0.3)) +
        scale_colour_manual(values = c("seagreen4", "goldenrod2", "indianred4"), guide = F) +
        coord_flip() +
        ylim(90, 525) +
        geom_text(aes(method, ll, group = bootMethod, color=method, label = round(ll)),
                  position= position_dodge(width = 0.35),
                  vjust = 0.5, hjust = 1.2, size = 2) +
        geom_text(aes(method, ul, group = bootMethod, color=method, label = round(ul)),
                  position= position_dodge(width = 0.35),
                  vjust = 0.5, hjust = -.2, size = 2) +
        geom_text(aes(method, estimate, label = round(estimate)),
                  #position = position_dodge(width = 0.3),
                  vjust = 2.5, hjust = .5, size = 3) +
        labs(title = "Fire cycle estimates and 95% confidence intervals for\n the case study area, Côte-Nord, Eastern Canada",
             y = "\nLength of fire cycle (years)",
             x = "") +
        theme_minimal()

dev.off()

