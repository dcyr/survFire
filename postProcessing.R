################################################################################
################################################################################
rm(list=ls())
require(raster)
################################################################################
################################################################################
setwd("/media/dcyr/Windows7_OS/Travail/SCF/fcEstimationExp")
outputFolder <- paste(getwd(), "compiledOutputs", sep="/")
wwd <- paste(getwd(), Sys.Date(), sep="/")
dir.create(wwd)
setwd(wwd)
rm(wwd)

### fetching total area
refRaster <- get(load("../outputs/simOutput_125_0.RData"))
refRaster <- refRaster[[1]]$tsf[[1]]
### removing 10-km edge (should check if that's coherent with upstream output compiling)
e <- extent(refRaster, 11, nrow(refRaster)-10, 11, ncol(refRaster)-10)
refRaster <- crop(refRaster, e)
totalArea <- ncell(refRaster)
totalArea <- ncell(refRaster) * prod(res(refRaster))/10000 ### total area in hectare

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
#
#
#

################################################################################
################################################################################
######
######      300-yrs simulations with replicates
######
require(dplyr)
require(zoo) ### for 'rollmean'
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


####################################################################
######  Plotting 300-yrs simulations results
require(ggplot2)
############################################

### illustrate the following graphically (with geom_histogram?)
trueFC %>%
    group_by(fireCycle, treatment) %>%
    summarise(meanTSF = mean(meanTSF),
              mean50 = mean(trueFC50),
              mean150 = mean(trueFC150),
              mean300 = mean(trueFC300))
###
rm(output)
###


############################################
##### loading FC estimation obtained from simulated field sampling experiment

survivalEstimates <- get(load(paste(outputFolder, "survivalEstimates.RData", sep ="/")))
survivalEstimates$sampleSize <- as.numeric(survivalEstimates$sampleSize)

#colnames(survivalEstimates)
#colnames(trueFC)

survivalEstimates <- merge(survivalEstimates, trueFC)
#head(survivalEstimates, 10)

df <- survivalEstimates
df <- survivalEstimates %>%
    filter(method == "cox")
# head(df)
residualsDF <- survivalEstimates %>%
    mutate(residual300 =  estimate - trueFC300,
           residual150 = estimate - trueFC150,
           residual50 = estimate - trueFC50) %>%
    group_by(fireCycle, treatment, sampleSize, method) %>%
    summarise(meanResidual300 = mean(residual300),
              meanResidual150 = mean(residual150),
              meanResidual50 = mean(residual50))


ylim = c(-200, 200)

boxPlotGraph <- ggplot(df, aes(x = sampleSize, y = estimate - trueFC300,
                                 group = sampleSize, fill = method)) +
    #colour = Species, order=as.numeric(Species))) +
    geom_boxplot() +
    #stat_sum_df("mean_cl_boot", geom = "smooth", mapping = aes(group = sampleSize)) + ## "mean_cl_boot' being a
    #stat_summary(fun.y="sum", geom="line", size=3) +
    # ylim(ylim) +
    #xlim(xlim) +
    ylim(ylim) +
    facet_grid(fireCycle ~ treatment, scales = "free") +
    #scale_colour_manual(values=cols) +
    #guides(fill = guide_legend(reverse = TRUE)) +
    labs(title = "Put title here",
         y = "Estimated - True fire cycle\n(residuals)",
         x = "Sample size")

png(filename=paste("biomassSpp", i, "_", a, "_SubReg", sr, ".png", sep=""),
    width = 2400, height = 2400, units = "in", pointsize=72)
print(g4 + theme_grey(base_size=72) +
          scale_x_continuous(breaks = breaks) +
          theme(axis.text.x = element_text(size=44, angle = 45, hjust = 1),
                axis.text.y = element_text(size=48),
                strip.text.x = element_text(size=40),
                strip.text.y = element_text(size=60)))
dev.off()

df[which(is.infinite(df$estimate)), "estimate"] <- 10000

summary(survivalEstimates)

ylim <- c(-200, 200)
xlim <- c(0, max(df$sampleSize)+10)
g4 <- ggplot(df, aes(x = sampleSize, y = estimate - trueFC300,
                     group = sampleSize)) +
                                    #colour = Species, order=as.numeric(Species))) +
    geom_boxplot(aes(fill = method)) +
    #stat_sum_df("mean_cl_boot", geom = "smooth", mapping = aes(group = sampleSize)) + ## "mean_cl_boot' being a
    #stat_summary(fun.y="sum", geom="line", size=3) +
   # ylim(ylim) +
    #xlim(xlim) +
    facet_grid(fireCycle ~ treatment, scales = "free") +
    #scale_colour_manual(values=cols) +
    #guides(fill = guide_legend(reverse = TRUE)) +
    labs(title = "Put title here",
         y = "Estimated - True fire cycle\n(residuals)",
         x = "Sample size")

png(filename=paste("biomassSpp", i, "_", a, "_SubReg", sr, ".png", sep=""),
    width = 2400, height = 2400, units = "in", pointsize=72)
    print(g4 + theme_grey(base_size=72) +
              scale_x_continuous(breaks = breaks) +
              theme(axis.text.x = element_text(size=44, angle = 45, hjust = 1),
                    axis.text.y = element_text(size=48),
                    strip.text.x = element_text(size=40),
                    strip.text.y = element_text(size=60)))
dev.off()


}
