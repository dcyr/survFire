################################################################################
################################################################################
rm(list=ls())
require(dplyr)
require(raster)
require(reshape2)
################################################################################
################################################################################
setwd("/media/dcyr/Windows7_OS/Travail/SCF/fcEstimationExp")
outputFolder <- paste(getwd(), "compiledOutputs", sep="/")
wwd <- paste(getwd(), Sys.Date(), sep="/")
dir.create(wwd)
setwd(wwd)
rm(wwd)

### fetching total area from an output raster stack
refRaster <- get(load("../outputs/simOutput_125_0.RData"))
refRaster <- refRaster[[1]]$tsf[[1]]
### removing 10-km edge (should check if that's coherent with upstream output compiling)
e <- extent(refRaster, 11, nrow(refRaster)-10, 11, ncol(refRaster)-10)
refRaster <- crop(refRaster, e)
totalArea <- ncell(refRaster)
totalArea <- ncell(refRaster) * prod(res(refRaster))/10000 ### total area in hectare


################################################################################
################################################################################
######
######      300-yrs simulations with replicates
######
require(dplyr)
require(zoo)
################################################################################
################################################################################

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

################################################################################
################################################################################
######
######      Confidence intervals coverage
######
require(ggplot2)
#require(quantreg)
require(RColorBrewer)
################################################################################
################################################################################

####################################################################
######  bootstrap estimates
survivalBootstrap <- get(load(paste(outputFolder, "survivalBootstrap.RData", sep="/")))

################################################################################
################################################################################


##### Filtering results that were deemed uninteresting through a trial and error process
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
t1 <- Sys.time()
survivalBootstrap <- merge(survivalBootstrap, trueFC)
t2 <- Sys.time()


fc  <- 250

# plot distribution
df <- survivalBootstrap %>%
    filter(resamplingEffort == 1) %>%
    #filter(treatment == 0) %>%
    filter(fireCycle == fc)

#head(df)


ciDensity <- ggplot(df, aes(x = estimate, colour = method)) +
    stat_density(geom="line",position="identity") +
    xlim(0, 3*fc) +
    facet_grid(sampleSize ~ treatment, scales = "free_y") +
    geom_vline(xintercept = fc, colour="grey", linetype = 3, size = 0.75) +
    labs(title = paste0("Bootstrap estimates"),
             y = "Density\n",
             x = "\nYears")


png(filename = paste0("CiDensity_", fc, ".png"),
    width = 10, height = 10, units = "in", res = 600, pointsize=10)

    print(ciDensity +
              theme_bw() +
              theme(legend.position="top", legend.direction="horizontal",
                    axis.text.x = element_text(angle = 45, hjust = 1)))

dev.off()


# merge with trueFC

# compute coverage

# plot coverage



# ################################################################################
# ################################################################################


# ################################################################################
# ################################################################################
# ######
# ######      Plotting 300-yrs simulations results
# ######
# require(ggplot2)
# require(quantreg)
# require(RColorBrewer)
# require(dichromat) ### for 'rollmean'
# ################################################################################
# ################################################################################
#
# ###
# trueFCSummary <- trueFC %>%
#     group_by(fireCycle, treatment) %>%
#     summarise(meanTSF = round(mean(meanTSF),1),
#               mean50 = round(mean(trueFC50), 1),
#               mean150 = round(mean(trueFC150), 1),
#               mean300 = round(mean(trueFC300),1))
# trueFCSummary <- trueFCSummary %>%
#     filter(fireCycle <= 500)
# trueFCSummary <- droplevels(trueFCSummary)
#
# ###
# df <- melt(trueFC, id.vars = c("fireCycle", "treatment", "replicate"),
#             measure.vars = c("trueFC50", "trueFC150", "trueFC300", "meanTSF"))
# df <- df %>%
#     filter(variable != "trueFC150") %>%
#     filter(fireCycle <= 500)
# df <- droplevels(df)
#
# ### renaming / reordering levels for nicer plotting
# # fireCycle
# fcLevels <- unique(paste0(df$fireCycle, "-yrs. FC"))
# df$fireCycle <- factor(paste0(df$fireCycle, "-yrs. FC"), levels = fcLevels)
# trueFCSummary$fireCycle <- factor(paste0(trueFCSummary$fireCycle, "-yrs. FC"), levels = fcLevels)
# # treatment
# tLevels <- c("Decreasing fire activity",
#              "Constant fire activity",
#              "Increasing fire activity")
# levels(df$treatment) <- tLevels
# df$treatment <- factor(df$treatment, levels = rev(tLevels))
# levels(trueFCSummary$treatment) <- tLevels
#
#
# ###
# df$variable <- factor(df$variable, levels = c("meanTSF", "trueFC300", "trueFC50"))
# df$variable <- factor(as.numeric(df$variable))
# levels(df$variable) <- c("Final mean TSF  ",
#                           "True FC (entire simulation)  ",
#                           "True FC (last 50 years)  ")
#
#
# #################
#
# nRep <- max(trueFC$replicate)
#
# hist_sim300 <- ggplot(df, aes(x = value, fill = variable)) +
#     geom_histogram(position="dodge", binwidth = 0.075) +
#     facet_grid(treatment ~ fireCycle, scales = "free_x") +
#     scale_x_log10(breaks = c(30, 62.5, 125, 250, 500, 1000, 2000, 4000)) +
#     geom_vline(data = trueFCSummary,  aes(xintercept = meanTSF),
#                colour="black", linetype = 3, size = 0.5) +
#     geom_vline(data = trueFCSummary,  aes(xintercept = mean300),
#                colour="darkolivegreen", linetype = 3, size = 0.75) +
#     geom_vline(data = trueFCSummary,  aes(xintercept = mean50),
#                colour="red3", linetype = 3, size = 0.75) +
#     scale_fill_manual("", values = c("black", "darkolivegreen", "red3")) +
#     labs(title = paste0("Realized fire activity for all 300-yrs simulations\n(",
#                         nRep, " replicates per treatment)"),
#      y = "Number of replicates\n",
#      x = "\nYears (log scale)")
#
#
# png(filename="realizedFC300.png",
#     width = 10, height = 7, units = "in", res = 600, pointsize=10)
#
#     print(hist_sim300 +
#               theme_bw() +
#               theme(legend.position="top", legend.direction="horizontal",
#                     axis.text.x = element_text(angle = 45, hjust = 1)))
#
# dev.off()
#
# ################################################################################
# ################################################################################

################################################################################
################################################################################
######
######      Plotting 300-yrs simulations residuals
######
##### loading FC estimation obtained from simulated field sampling experiment
survivalEstimates <- get(load(paste(outputFolder, "survivalEstimates.RData", sep ="/")))
survivalEstimates$sampleSize <- as.numeric(survivalEstimates$sampleSize)
################################################################################
################################################################################


##### Filtering results that were deemed uninteresting through a trial and error process
##### (Usually because it creates CI so wide they were useless)
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


# ################################################################################
# ##### Boxplot - Residuals (Residuals by SampleSize, for each fire cycles and treatments)
# ################################################################################
# ### Computing residuals to plot
# residualsDF <- survivalEstimates %>%
#     mutate(residual300 =  estimate - trueFC300,
#            residual150 = estimate - trueFC150,
#            residual50 = estimate - trueFC50) %>%
#     group_by(fireCycle, treatment, sampleSize, method) %>%
#     summarise(meanResidual300 = round(mean(residual300), 1),
#               meanResidual150 = round(mean(residual150), 1),
#               meanResidual50 = round(mean(residual50), 1),
#               p005Residual300 = round(quantile(residual300, 0.005), 1),
#               p025Residual300 = round(quantile(residual300, 0.025), 1),
#               p05Residual300 = round(quantile(residual300, 0.05), 1),
#               p50Residual300 = round(quantile(residual300, 0.5), 1),
#               p95Residual300 = round(quantile(residual300, 0.95), 1),
#               p975Residual300 = round(quantile(residual300, 0.975), 1),
#               p995Residual300 = round(quantile(residual300, 0.995), 1))
#
# ### renaming factors for nice plotting
# fcLevels <- unique(paste0(residualsDF$fireCycle, "-yrs. FC"))
# residualsDF$fireCycle <- factor(paste0(residualsDF$fireCycle, "-yrs. FC"), levels = fcLevels)
# residualsDF$treatment <- factor(residualsDF$treatment)
# levels(residualsDF$treatment) <- c("Decreasing fire activity",
#                                    "Constant fire activity",
#                                    "Increasing fire activity")
# residualsDF$treatment <- factor(residualsDF$treatment, levels = rev(levels(residualsDF$treatment)))
# residualsDF$fireCycle <- factor(residualsDF$fireCycle)
#
# ### The 'data.frame' to plot
# df <- residualsDF %>%
#     filter(method %in% c("cox", "weib", "exp"))
# df <- droplevels(df)
# ### renaming 'method'
# levels(df$method) <- c("Cox", "Weibull", "Exponential")
#
# ###############################################################
# ### Boxplot
#
# boxPlotGraph <- ggplot(df, aes(x = factor(sampleSize))) +#,color = method)) +
#     geom_hline(yintercept = 0, size = 0.25, color = "grey") +
#     geom_boxplot(aes(ymin = p025Residual300,
#                      lower = p05Residual300,
#                      middle = meanResidual300,
#                      upper = p95Residual300,
#                      ymax = p975Residual300,
#                      #group = sampleSize,
#                      fill = method),
#                  position = position_dodge(width = 0.5),
#                  width = 0.5,
#                  stat = "identity",
#                  size = 0.3,
#                  alpha = 1) + #
#     facet_grid(fireCycle ~ treatment, scales = "free") +
#     labs(title = "Fire cycle estimation error from simulated forest fire history\nreconstruction affected by censoring\n",
#          y = "Estimate - True fire cycle\n(residuals)",
#          x = "\nSample size")
#
# ###
#
# png(filename="residualsCens.png",
#     width = 10, height = 6, units = "in", res = 600, pointsize=10)
#
#     print(boxPlotGraph +
#               scale_fill_brewer(type = "qual") +#, palette = 1, direction = 1)
#
#
#               theme_bw())
# # print(g4 + theme_grey(base_size=72) +
# #           scale_x_continuous(breaks = breaks) +
# #           theme(axis.text.x = element_text(size=44, angle = 45, hjust = 1),
# #                 axis.text.y = element_text(size=48),
# #                 strip.text.x = element_text(size=40),
# #                 strip.text.y = element_text(size=60)))
# dev.off()
# ###############################################################
# ###############################################################



# ################################################################################
# ################################################################################
# ######
# ######      Fire-size distribution figure
# ######
# library(MASS)
# ################################################################################
# ################################################################################
#
# ### empirical reference
# fireObs <- read.csv("../data/fireObs.csv", header=TRUE)
# fireSizeObs <- fireObs[, "SIZE"]
# fireSizeObs <- fireSizeObs[order(fireSizeObs)]
# ### fitted fire size distribution
# fireSizeFit <- fitdistr(fireSizeObs, "lognormal")
# #########
# xlim <- c(0, 100000)
# options(scipen=7) ###force fixed notation (no scientific notation unless > 10E7 )
#
# png(filename = paste0("fireSizeCDF.png"),
#     width = 4, height = 3.5, units = "in", res = 300, pointsize = 8,
#     bg = "white")
#
#     ### fitted fire size cdf
#     plot(plnorm(1:max(fireSizeObs), mean = fireSizeFit$estimate[1],
#                  sdlog = fireSizeFit$estimate[2]),
#          type = "l", xlim = xlim,
#          main = "Fire size cumulative probability distribution",
#          ylab = "Cumulative probability",
#          xlab = "Fire size (ha)", lty = 3)
#     ### empirical fire size cdf
#     n <- length(fireSizeObs)
#     lines(fireSizeObs, (1:n)/n, type = 's')
#     ###
#     legend(x = 0.7*xlim[2], y = 0.9,
#            c("Empirical CDF", "Fitted CDF", "Individual fires"),
#            lty = c(1, 3, 1), lwd = 1, col = c("black", "black", "indianred"),
#            cex=0.75,
#            bg="white")
#     ###
#     rug(fireSizeObs, ticksize = 0.03, side = 1, lwd = 1, col = "indianred")
#     abline(h = c(0,1), lwd = 0.2, lty = 2, col = "grey")
#
# dev.off()
# ######
# ################################################################################
# ################################################################################
#
#
#
# ################################################################################
# ################################################################################
# ######
# ######      2000-yrs constant fc simulations
# ######
# require(dplyr)
# require(zoo) ### for 'rollmean'
# ################################################################################
# ################################################################################
#
# ####################################################################
# ######  true fire cycle statistics computation
# output <- get(load(paste(outputFolder, "simOutputCompiled2000.RData", sep="/")))
# ######
#
# output <- output %>%
#     filter(replicate == 1) %>%
#     mutate(id = paste0(fireCycle, replicate, treatment))
# output$treatment <- as.factor(output$treatment)
#
# output <- output %>%
#     arrange(year, fireCycle)
#
# trueFCNames <- c("trueFC50", "trueFC150", "trueFC300")
# output[, trueFCNames] <- NA
#
# trueFC <- list()
# for (fc in unique(output$fireCycle)) {
#     trueFC[[fc]] <- output %>%
#         filter(fireCycle == fc)
#     year <- trueFC[[fc]]$year
#     aab <- trueFC[[fc]]$areaBurned_ha
#     ###
#     trueFC[[fc]][,"trueFC50"] <- totalArea/rollmean(aab, 50, align = "right", fill = NA)
#     trueFC[[fc]][,"trueFC150"] <- totalArea/rollmean(aab, 150, align = "right", fill = NA)
#     trueFC[[fc]][,"trueFC300"] <- totalArea/rollmean(aab, 300, align = "right", fill = NA)
# }
#
#
# ####################################################################
# ######  plotting 2000-yrs simulation examples
#
# yLimAAB <- c(0, max(as.numeric(lapply(trueFC, function(x) max(x$areaBurned_ha)))))
#
# for (fc in c(62.5, 125, 250, 500, 1000)) {
#     df <- trueFC[[fc]]
#     x <- df$year
#     meanTSF <- df$meanTSF
#     FC50 <- df$trueFC50
#     FC150 <- df$trueFC150
#     FC300 <- df$trueFC300
#     aab <- df$areaBurned_ha
#     ###
#     realizedFC <- totalArea/mean(aab)
#
#
#     png(filename = paste0("simTrueTSF", fc, ".png"),
#         width = 8, height = 5, units = "in", res = 300, pointsize = 8,
#         bg = "white")
#
#         #########
#         par(mfrow=c(2,1), mar=c(0,4,3,1))
#         options(scipen=7) ###force fixed notation (no scientific notation unless > 10E7 )
#         #
#         plot(x = x, y = aab, type="l", ylab="Annual area burned (ha)", xaxt="n", lwd=0.5, xlim=c(0,2000), ylim=yLimAAB)
#         abline(v = seq(from=0, to=2000, by= 100), col = "gray", lty=3, lwd=.5)
#         grid(col="gray", lty=3, lwd=.8)
#         #
#         text(max(x), y = yLimAAB[2]*.98, paste("Realized fire cycle (entire simulation) :", round(realizedFC, 1), "years"), adj = c(1,1), cex=1)
#         text(max(x), y = yLimAAB[2]*.92, paste("Mean percent annual area burned :", round(100*1/realizedFC, 2), "%"), adj = c(1,1), cex=1)
#         ###
#         par(mar=c(5,4,0,1))
#         ylim <- c(0,(realizedFC*3.5))
#         #ylimits <- c(0, 300)
#         plot(meanTSF, type="l", lwd=4, lty=1, col="black", ylab="Mean time since fire / True fire cycle", xlab="Simulated years", ylim=ylim)
#         abline(v = seq(from=0, to=2000, by= 100), col = "gray", lty=3, lwd=.5)
#         grid(col="gray", lty=3, lwd=.8)
#         abline(h=realizedFC, lty=1, lwd=1, col="black")
#         lines(FC50, lty=1, lwd=2, col="red3")
#         lines(FC150, lty=1, lwd=2, col="dodgerblue2")
#         lines(FC300, lty=1, lwd=2, col="darkolivegreen")
#         legend(0, ylim[2],
#                c("True mean TSF", "True FC (past 50 years)", "True FC (past 150 years)", "True FC (past 300 years)"),
#                lty = 1, lwd = c(4,2,2,2), col = c("black", "red3", "dodgerblue2", "darkolivegreen"),
#                cex=1,
#                bg="white")
#
#     dev.off()
# }
# ######
# ################################################################################
# ################################################################################