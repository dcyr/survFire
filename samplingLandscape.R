####################################################################
####################################################################
####################################################################
rm(list=ls())
####################################################################
####################################################################
setwd("/media/dcyr/Windows7_OS/Travail/SCF/fcEstimationExp")
outputFolder <- paste(getwd(), "outputs", sep="/")
wwd <- paste(getwd(), Sys.Date(), sep="/")
dir.create(wwd)
setwd(wwd)
rm(wwd)

####################################################################
####################################################################
######
######      sampling landscape
######
####################################################################
####################################################################

x <- list.files(outputFolder)
simInfo <- gsub(".RData", "", x)
simInfo <- strsplit(simInfo, "_")
fc <- as.numeric(lapply(simInfo, function(x) x[2]))
corr <- as.numeric(lapply(simInfo, function(x) x[3]))

require(doSNOW)
clusterN <- 3  ### choose number of nodes to add to cluster.

cl = makeCluster(clusterN)
registerDoSNOW(cl)

###########################################
fireCycle <- meanTSF <- year <- replicate <- areaBurned <- treatment <- list()
###########################################
outputList <- foreach(i = seq_along(x)) %dopar%  {
    #for (i in seq_along(x)) {
    require(raster)
    ##
    output <- get(load(paste(outputFolder, x[i], sep="/")))
    fireCycle <- meanTSF <- year <- replicate <- areaBurned <- treatment <- list()
    ##
    for (r in seq_along(output)) {

        tsfStack <- output[[r]][["tsf"]]
        convFactor <- prod(res(tsfStack))/10000### to convert to hectares
        ## remove edge effet
        e <- extent(tsfStack, 11, nrow(tsfStack)-10, 11, ncol(tsfStack)-10)
        tsfStack <- crop(tsfStack, e)
        w <- tsfStack[[1]]
        w[] <- 1
        ## compiling statistics
        meanTSF[[r]] <- round(as.numeric(zonal(tsfStack, w, mean))[-1], 1)
        areaBurned[[r]] <- as.numeric(zonal(tsfStack==0, w, sum))[-1] * convFactor
        year[[r]] <- 1:nlayers(tsfStack)
        replicate[[r]] <- rep(r, nlayers(tsfStack))
        fireCycle[[r]] <- rep(fc[i], nlayers(tsfStack))
        treatment[[r]] <- rep(corr[i], nlayers(tsfStack))
        rm(tsfStack)
        print(r)
    }
    ##
    fireCycle <- as.numeric(do.call("cbind", fireCycle))
    treatment <- as.numeric(do.call("cbind", treatment))
    replicate <- as.numeric(do.call("cbind", replicate))
    year <- as.numeric(do.call("cbind", year))
    areaBurned_ha <- as.numeric(do.call("cbind", areaBurned))
    meanTSF <- as.numeric(do.call("cbind", meanTSF))
    ##
    output <- data.frame(fireCycle, treatment, replicate, year, areaBurned_ha, meanTSF)
    return(output)
    #print(paste(i, "of", length(x)))
}

stopCluster(cl)
output <- do.call("rbind", outputList)
save(output, file = "simOutputCompiled.RData")



#
#
# for (i in )
#
# meanTSF <- as.numeric(do.call('cbind', meanTSF))
# areaBurned <- as.numeric(do.call('cbind', areaBurned))
# year <- as.numeric(do.call('cbind', year))
# replicate <- as.numeric(do.call('cbind', replicate))
# fireCycle <- as.numeric(do.call('cbind', fireCycle))
#
#
# outputSummary
#
#
#
# require(ggplot2)
# require(dplyr)
# df <- as.data.frame(outputSummary) %>%
#     #filter(replicate == 1) %>%
#     mutate(areaBurned = areaBurned * 100)
#
# mean(df$meanTSF)
#
# g <- ggplot(df, aes(x=year, y=areaBurned,
#                     group=replicate)) +
#     #stat_summary(fun.data = "median_hilow", geom = "smooth", size=5) +
#     geom_line(aes(group=replicate),
#               size=1,
#               type=1,
#               alpha=0.5)
#
#     +
#     scale_colour_brewer(palette="Set2") +
#     scale_fill_brewer(palette="Set2") +
#     facet_grid(Scenario ~ Treatment) +
#     guides(colour = guide_legend(reverse = TRUE),
#            fill = guide_legend(reverse = TRUE)) +
#     xlim(c(2000, 2100)) +
#     labs(title = paste("Fire activity - ", ecoName, "\n", sep=""),
#          y=ifelse(i=="PAAB", "Annual area Burned (%)\n", "Mean fire return interval (years)\n"),
#          x="Year")


#samplingEffort <- 100
#binWidth <- c(5,25,50,75,100,125,150,175,200,250,300)
#censInt <-c(100,300)
#spinOffDuration <- 2 * max(fireCycle)
