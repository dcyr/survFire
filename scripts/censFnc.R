#########
### censoring function
#########
censFnc <- function(x, min, max) {
    time <-runif(length(x), min = min, max = max)
    cens <- as.numeric(x<time)
    x[cens == 0] <- round(time[cens == 0]) ## 0=> censored (alive), 1 => not censored (death)
    return(Surv(time = x, event = cens))
}