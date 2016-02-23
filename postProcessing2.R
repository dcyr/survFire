################################################################################
################################################################################
rm(list=ls())
require(dplyr)
require(raster)
require(reshape2)
require(zoo)
################################################################################
################################################################################
setwd("/media/dcyr/Windows7_OS/Travail/SCF/fcEstimationExp")
outputFolder <- paste(getwd(), "compiledOutputs", sep="/")
wwd <- paste(getwd(), Sys.Date(), sep="/")
dir.create(wwd)
setwd(wwd)
rm(wwd)

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
### choosing between censored and uncensored data sets
#FcEstimatedMethod <- c("cox", "weib", "exp")
FcEstimatedMethod <- c("cox", "weib", "exp")

####################################################################
######  loading bootstrap estimates
survivalEstimates <- get(load(paste(outputFolder, "survivalEstimates.RData", sep ="/")))
survivalEstimates <- survivalEstimates %>%
    filter(method %in% FcEstimatedMethod) %>%
    filter(sampleSize >= 25)
survivalEstimates <- droplevels(survivalEstimates)
################################################################################

################################################################################
######  loading bootstrap estimates
survivalBootstrap <- get(load(paste(outputFolder, "survivalBootstrap.RData", sep="/")))
survivalBootstrap <- filter(survivalBootstrap, method %in% FcEstimatedMethod)
survivalBootstrap <- droplevels(survivalBootstrap)


### defining colors, renaming / reordering levels for nicer plotting
mColors <- c("seagreen4", "goldenrod2", "indianred4")
# fireCycle
fcLevels <- unique(survivalEstimates$fireCycle)
fcLevels <- fcLevels[order(fcLevels)]
fcLevels <- paste0(fcLevels, "-yrs. FC")
survivalEstimates$fireCycle <- factor(paste0(survivalEstimates$fireCycle, "-yrs. FC"), levels = fcLevels)
survivalBootstrap$fireCycle <- factor(paste0(survivalBootstrap$fireCycle, "-yrs. FC"), levels = fcLevels)
# treatment
tLevels <- c("Decreasing fire activity",
             "Constant fire activity",
             "Increasing fire activity")
survivalEstimates$treatment <- as.factor(survivalEstimates$treatment)
survivalBootstrap$treatment <- as.factor(survivalBootstrap$treatment)
levels(survivalEstimates$treatment) <- levels(survivalBootstrap$treatment) <- tLevels
survivalEstimates$treatment <- factor(survivalEstimates$treatment, levels = rev(tLevels))
survivalBootstrap$treatment <- factor(survivalBootstrap$treatment, levels = rev(tLevels))

# sampleSize
ssLevels <- unique(survivalEstimates$sampleSize)
ssLevels <- ssLevels[order(ssLevels)]
ssLevels <- paste("sample size:", ssLevels)
survivalEstimates$sampleSize <- factor(paste("sample size:", survivalEstimates$sampleSize), levels = ssLevels)
survivalBootstrap$sampleSize <- factor(paste("sample size:", survivalBootstrap$sampleSize), levels = ssLevels)


# method
levels(survivalEstimates$method) <-
    levels(survivalBootstrap$method) <- c("Cox", "Weibull", "Exponential")

### plot distribution (suppl. material)
for (fc in unique(survivalBootstrap$fireCycle)) {
    fcNum <- as.numeric(gsub("-yrs. FC", "", fc))
    # plot distribution
    df <- survivalBootstrap %>%
        #filter(treatment == 0) %>%
        filter(fireCycle == fc) %>%
        mutate(residual = estimate - trueFC300)

    df <- droplevels(df)

    ### residual summary (to plot mean values)
    residualsSummary <- df %>%
        group_by(fireCycle, treatment, sampleSize, method) %>%
        summarise(meanResidual = round(mean(residual),1))

    ### Plot
    ciDensity <- ggplot(df, aes(x = residual, colour = method)) +
        geom_vline(xintercept = 0, colour="grey", linetype = 3, size = 0.75) +
        stat_density(geom="line", position="identity", size = 1) +
        xlim(-fcNum, fcNum) +
        facet_grid(sampleSize ~ treatment, scales = "free_y") +

        geom_vline(data = residualsSummary,  aes(xintercept = meanResidual, colour = method),
                   linetype = 3, size = 1) +
        scale_colour_manual(values = mColors) +
        #scale_colour_brewer(type = "qual") +
        labs(title = paste0("Distribution of bootstrap resampled fire cycle estimates\n(residuals; ", fcNum, "; n = 10000)"),
             y = "Density\n",
             x = "Residuals\nEstimated FC - True value (years)")

    ### Printing plot
    png(filename = paste0("coverageDensity_", fcNum, "UnCensored.png"),
        width = 10, height = 10, units = "in", res = 600, pointsize=10)

        print(ciDensity +
                  theme_bw() +
                  theme(legend.position="top", legend.direction="horizontal",
                        axis.text.x = element_text(angle = 45, hjust = 1),
                        strip.text.y = element_text(size = 8))) #, palette = 1, direction = 1)

    dev.off()
}


### computing bootstrap 95CI
coverageDf <- survivalBootstrap %>%
    mutate(ll = ll-estimate,
           ul = ul-estimate) %>%
    group_by(fireCycle, treatment, sampleSize, method) %>%
    summarise(meanll = mean(ll),
              meanul = mean(ul),
              medianll = median(ll),
              medianul = median(ul))



coverageDf[,c("meanll", "meanul")]
hist(coverageDf$meanul, xlim=c(-500, 500))

summary(coverageDf$ll - coverageDf$medianll)

### computing coverage rate
coverageDf <- merge(survivalEstimates, coverageDf)

#
coverageSummary <- coverageDf %>%
    mutate(residual = estimate - trueFC300) %>%
    group_by(fireCycle, treatment, sampleSize, method) %>%
    summarise(coverage = sum(ll < residual & ul > residual)/n())

# converting sampleSize back to numerical values
coverageSummary$sampleSize <- as.numeric(gsub("sample size: ", "", coverageSummary$sampleSize))

### plotting
coverageSummaryPlot <- ggplot(coverageSummary, aes(x = sampleSize, y = coverage, colour = method)) +
    facet_grid(fireCycle ~ treatment) +
    geom_hline(yintercept = 0.95, linetype = "dashed", size = 0.5, col = "grey") +
    geom_line() +
    ylim(0.75, 1) +
    scale_colour_manual(values = c("seagreen4", "goldenrod2", "indianred4")) +
    #?scale_colour_brewer(type = "qual") +
    scale_linetype_manual(values=c("dotted", "solid", "twodash")) +
    labs(title = paste0("Coverage rate of bootstrap 95% confidence intervals\n(10000 resampling)\n"),
         y = "Coverage rate\n",
         x = "\nSample size")


### printing plot
png(filename = paste0("coverageUncensored.png"),
    width = 10, height = 10, units = "in", res = 600, pointsize=10)

    print(coverageSummaryPlot +
              theme_bw() +
              theme(axis.text.x = element_text(angle = 45, hjust = 1),
                    #legend.position="right", legend.direction="horizontal",
                    strip.text.y = element_text(size = 10))) #, palette = 1, direction = 1)

dev.off()


################################################################################
################################################################################


################################################################################
################################################################################
######
######      Plotting 300-yrs simulations results
######
require(ggplot2)
require(quantreg)
require(RColorBrewer)
require(dichromat) ### for 'rollmean'
################################################################################
################################################################################

###
trueFC <- distinct(survivalEstimates[, c("fireCycle", "treatment", "replicate", "meanTSF",
                                      "trueFC50", "trueFC150", "trueFC300", "meanAAB")])
trueFCSummary <- trueFC %>%
    group_by(fireCycle, treatment) %>%
    summarise(meanTSF = round(mean(meanTSF),1),
              mean50 = round(mean(trueFC50), 1),
              mean150 = round(mean(trueFC150), 1),
              mean300 = round(mean(trueFC300),1))

###
df <- melt(trueFC, id.vars = c("fireCycle", "treatment", "replicate"),
            measure.vars = c("trueFC50", "trueFC150", "trueFC300", "meanTSF"))
df <- df %>%
    filter(variable != "trueFC150")
df <- droplevels(df)

### renaming / reordering levels for nicer plotting

###
df$variable <- factor(df$variable, levels = c("meanTSF", "trueFC300", "trueFC50"))
df$variable <- factor(as.numeric(df$variable))
levels(df$variable) <- c("Final mean TSF  ",
                          "True FC (entire simulation)  ",
                          "True FC (last 150 years)  ")


#################

nRep <- max(trueFC$replicate)

hist_sim300 <- ggplot(df, aes(x = value, fill = variable)) +
    geom_histogram(position="dodge", binwidth = 0.075) +
    facet_grid(treatment ~ fireCycle, scales = "free_x") +
    scale_x_log10(breaks = c(30, 62.5, 125, 250, 500, 1000, 2000, 4000)) +
    geom_vline(data = trueFCSummary,  aes(xintercept = meanTSF),
               colour="black", linetype = 3, size = 0.5) +
    geom_vline(data = trueFCSummary,  aes(xintercept = mean300),
               colour="darkolivegreen", linetype = 3, size = 0.75) +
    geom_vline(data = trueFCSummary,  aes(xintercept = mean50),
               colour="red3", linetype = 3, size = 0.75) +
    scale_fill_manual("", values = c("black", "darkolivegreen", "red3")) +
    labs(title = paste0("Realized fire activity for all 300-yrs simulations\n(",
                        nRep, " replicates per treatment)"),
     y = "Number of replicates\n",
     x = "\nYears (log scale)")


png(filename="realizedFC300.png",
    width = 10, height = 7, units = "in", res = 600, pointsize=10)

    print(hist_sim300 +
              theme_bw() +
              theme(legend.position="top", legend.direction="horizontal",
                    axis.text.x = element_text(angle = 45, hjust = 1)))

dev.off()

################################################################################
################################################################################

################################################################################
################################################################################
######
######      Plotting 300-yrs simulations residuals
######
################################################################################
###############################################################################

################################################################################
##### Boxplot - Residuals (Residuals by SampleSize, for each fire cycles and treatments)
################################################################################
### Computing residuals to plot
residualsDF <- survivalEstimates %>%
    mutate(residual300 =  estimate - trueFC300,
           residual150 = estimate - trueFC150,
           residual50 = estimate - trueFC50) %>%
    group_by(fireCycle, treatment, sampleSize, method) %>%
    summarise(meanResidual300 = round(mean(residual300), 1),
              meanResidual150 = round(mean(residual150), 1),
              meanResidual50 = round(mean(residual50), 1),
              p005Residual300 = round(quantile(residual300, 0.005), 1),
              p025Residual300 = round(quantile(residual300, 0.025), 1),
              p05Residual300 = round(quantile(residual300, 0.05), 1),
              p50Residual300 = round(quantile(residual300, 0.5), 1),
              p95Residual300 = round(quantile(residual300, 0.95), 1),
              p975Residual300 = round(quantile(residual300, 0.975), 1),
              p995Residual300 = round(quantile(residual300, 0.995), 1))

### renaming sample size
ss <- as.factor(gsub("sample size: ", "", residualsDF$sampleSize))
ss <- factor(ss, levels = levels(ss)[order(as.numeric(levels(ss)))])
residualsDF$sampleSize <- ss

###############################################################
### Boxplot

boxPlotGraph <- ggplot(residualsDF, aes(x = factor(sampleSize))) +#,color = method)) +
    geom_hline(yintercept = 0, size = 0.25, color = "grey") +
    geom_boxplot(aes(ymin = p025Residual300,
                     lower = p05Residual300,
                     middle = meanResidual300,
                     upper = p95Residual300,
                     ymax = p975Residual300,
                     #group = sampleSize,
                     fill = method),
                 position = position_dodge(width = 0.5),
                 width = 0.5,
                 stat = "identity",
                 size = 0.3,
                 alpha = 1) + #
    facet_grid(fireCycle ~ treatment, scales = "free") +
    labs(title = "Fire cycle estimation error from simulated forest fire history\nreconstruction affected by censoring\n",
         y = "Estimate - True fire cycle\n(residuals)",
         x = "\nSample size")

###
png(filename="residualsCens.png",
    width = 10, height = 6, units = "in", res = 600, pointsize=10)

    print(boxPlotGraph +
              scale_fill_manual(values = mColors) +
              #scale_fill_brewer(type = "qual") +#, palette = 1, direction = 1)


              theme_bw())
dev.off()
###############################################################
###############################################################



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