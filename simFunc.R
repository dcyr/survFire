####################################################################
####################################################################
###### simulations function
####################################################################
####################################################################

sim <- function(initialLandscape, simDuration, fireCycle, fireSizeFit,
                distribType = NULL, outputTSF = TRUE, outputFire = FALSE) {

    require(raster)
    source("../fireSpreadFunc.R")

    scaleFactor <- prod(res(initialLandscape))/10000 # converting number of pixels to hectares

    ##############################################################
    ##############################################################
    ##### fire regime
    ##############################################################
    ##############################################################
    # average fraction of the "burnable" area burned annually
    fAAB <- 1 / fireCycle
    # fire size distrib
    estimates <- names(fireSizeFit$estimate)
    if (estimates == "rate") {
        distribType <- "exp"
        fireSizeMean <- 1/fireSizeFit$estimate["rate"]
    }
    if (estimates[1] == "meanlog" & estimates[2] == "sdlog") {
        distribType <- "lognorm"
        fireSizeMean <- exp(fireSizeFit$estimate["meanlog"] + 0.5*fireSizeFit$estimate["sdlog"]^2)
    }
    # fire sequence (annually)
    nFiresMean <- fAAB * sum(values(initialLandscape>0), na.rm = T) * scaleFactor / fireSizeMean
    nFireSequence <- rpois(lambda = nFiresMean, n = 1:simDuration)

    ##############################################################
    ##############################################################
    ##### aging and burning
    ##############################################################
    ##############################################################
    fires <- list()
    for (y in seq_along(nFireSequence)) {

        ### aging landscape
        if (y == 1) {
            tsf <- initialLandscape
            timeSinceFire <- list()
        } else {
            tsf <- tsf + 1
        }
        ### burning
        nFires <- nFireSequence[y]
        ### skip if no fire events
        if (nFires > 0) {
            eligible <- tsf > 0

            ## burning
            if (distribType == "exp") {
                fSize <- round(rexp(nFires, rate = fireSizeFit$estimate))
            }
            if (distribType == "lognorm") {
                fSize <- round(rlnorm(nFires, meanlog = fireSizeFit$estimate[1],
                                sdlog = fireSizeFit$estimate[2]))
            }

            f <- stack(fireSpread(eligible = eligible, fireSize = fSize))
            if (nlayers(f) == 1) {
                tsf[f > 0] <- 0
            } else {
                tsf[calc(f, sum, na.rm = T) > 0] <- 0
            }
            fires[[y]] <- f
        }

        timeSinceFire[[y]] <- tsf

        t2 <- Sys.time()
        print(paste0("year ", y, " (", nFires, " fires)", round(t2-t1, 2)))
    }
    ## preparing outputs
    output <- list()
    if (outputFire){
        output[["fires"]] <- fires
    }
    if (outputTSF){
        output[["tsf"]] <- timeSinceFire
    }
    return(output)
}

