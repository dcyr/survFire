################################################################################
################################################################################
####
#### A simple cellular automaton to simulate fire spread
#### Dominic Cyr
####
################################################################################
################################################################################

fireSpread <- function(eligible, w = NULL,
                       fireSize, probSpread =  0.32) {

    require(raster)
    scaleFactor <- prod(res(eligible) / 100) # convert pixels into hectares

    ##default propagation kernel ("queen's case")
    if(is.null(w))  {
        w <- matrix(1, nrow = 3, ncol = 3)
    }

    f <- eligible

    fArea <- list()
    for (i in seq_along(fireSize)) {
        f[] <- 0
        fs <- fireSize[i]
        ## igniting fire
        ignition <- sample(which(values(eligible)>0), 1)
        f[ignition] <- 1
        fireFront <- f
        #
        fSize <- 1*scaleFactor

        while(fSize < fs) {
            remaining <- (fs - fSize)/scaleFactor
            ### looking up neighbours
            # all cells adjacent to a burned cell
            fireNeighbour <- focal(fireFront, w = w, pad = T, padValue = 0) > 0
            # substraction those that already burned
            fireNeighbour <- fireNeighbour - (f>0)
            indices <- which(values(fireNeighbour) == 1)
            ### remove pixels that already burned
            indices <- setdiff(indices, which(values(eligible) == 0))
            ### fire spread
            # burning only a proportion of neighbours
            rVec <- runif(length(indices))
            indices <- indices[union(which.min(rVec), which(rVec < probSpread))] ## burn at least one pixel until target fire size is achieved
            # trim so that final fire size does not exceed target fire size
            if (length(indices) > remaining) {
                indices <- sample(indices, size = remaining)
            }


            if (length(indices) == 0) break  ### fire stops burning

            f[indices] <- 1 ##
            # resetting fire front
            fireFront[] <- 0
            # updating fire front for next iteration
            fireFront[indices] <- 1 ## area burned in the current time step
            ## updating current fire size
            fSize <- fSize + length(indices)*scaleFactor
            # record iteration on raster
            f[f>0] <- f[f>0] + 1
        }
        f[f==0] <- NA
        # storing individual fire events
        fArea[[i]] <- f-1 ## spatial
        eligible[f>0] <- 0
    }
    return(stack(fArea))
}

