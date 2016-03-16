####################################################################
####################################################################
###### Simulation function
######
###### Runs one simulation based of inputs provided upstream
###### Outputs TSF maps and/or individual fire map
######
###### Dominic Cyr
####################################################################
####################################################################

sim <- function(initialLandscape, simDuration, fireCycle, fireSizeFit,
                distribType = NULL, outputTSF = TRUE, outputFire = FALSE,
                corr = 0) {


    if(abs(corr)>0.75) stop("correlation should be between -0.75 and 0.75")

    require(raster)
    source("../scripts/fireSpreadFnc.R")

    #### converting number of pixels to hectares
    scaleFactor <- prod(res(initialLandscape))/10000

    ##############################################################
    ##############################################################
    ####  fire regime
    require(MASS)
    ##############################################################
    ##############################################################

    ### average fraction of the "burnable" area burned annually
    fAAB <- 1 / fireCycle

    ########
    ### fire size distribution
    estimates <- names(fireSizeFit$estimate)
    if (estimates == "rate") {
        distribType <- "exp"
        fireSizeMean <- 1/fireSizeFit$estimate["rate"]
    }
    if (estimates[1] == "meanlog" & estimates[2] == "sdlog") {
        distribType <- "lognorm"
        fireSizeMean <- exp(fireSizeFit$estimate["meanlog"] +
                                0.5*fireSizeFit$estimate["sdlog"]^2)
    }

    ########
    ### Fire sequence (annually)
    nFiresMean <- fAAB * sum(is.na(values(initialLandscape))==FALSE) * scaleFactor / fireSizeMean
    ### Generate a yearly sequence of number of fires
    nFireSequence <- rpois(lambda = nFiresMean, n =simDuration)


    if (!(corr == 0)) { # Imposing correlation between number of fires and sim year
        # Ordering sequence, then permutating until desired correlation is achieved
        repeat {
            nFireSequence <- nFireSequence[order(nFireSequence,
                                                 decreasing = ifelse(corr > 0, F, T))]
            corrTmp <- cor(nFireSequence, 1:simDuration)

            iter <- 1
            while(abs(corrTmp - corr)>0.01) {
                x <- sample(1:simDuration, 2)
                nFireSequence[x] <- nFireSequence[rev(x)]
                #nFireSequence[i[2]] <- nFireSequence[i[1]]
                corrTmp <- cor(nFireSequence, 1:simDuration)
                iter <- iter+1
                if (iter > 1000) break # in case algo misses target correlation
            }
            if  (abs(corrTmp - corr) < 0.01) break  # This acceptable difference may need to be larger when simDuration is small
        }
    }

    ##############################################################
    ##############################################################
    ##### aging and burning
    ##############################################################
    ##############################################################
    fires <- list()
    for (y in seq_along(nFireSequence)) {

        #### aging landscape
        if (y == 1) {
            tsf <- initialLandscape
            timeSinceFire <- list()
        } else {
            tsf <- tsf + 1
        }
        #### burning
        nFires <- nFireSequence[y]
        #### skip if no fire events
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
            if (sum(values(eligible>0)) > 0) { ## skip if there's nothing eligible to burn
                f <- stack(fireSpread(eligible = eligible, fireSize = fSize))
                if (nlayers(f) == 1) {
                    tsf[f > 0] <- 0
                } else {
                    tsf[calc(f, sum, na.rm = T) > 0] <- 0
                }
                fires[[y]] <- f
            }

        }
        timeSinceFire[[y]] <- tsf

        print(paste0("year ", y, " (", nFires, " fires)"))
    }
    ## preparing outputs
    output <- list()
    if (outputFire){ ### these outputs record all individual fires (larger)
        output[["fires"]] <- fires
    }
    if (outputTSF){ ### these outputs record TSF maps
        output[["tsf"]] <- stack(timeSinceFire)
    }
    return(output)
}

