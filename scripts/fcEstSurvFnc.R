#######################################
######## fire cycle estimation survival function
#######################################

expFitFn <- function(x) {        ### negative exponential fitting
    surv.exp <- survreg(x ~ 1, dist="exponential")
    b <- exp(surv.exp$coefficients)
    list(cycle=b, freq=1/b)
}

weibFitFn <- function(x) {        ### weibull exponential fitting
    surv.weib <- survreg(x ~ 1, dist="weibull")
    b <- exp(surv.weib$coefficients)
    c <- 1/(surv.weib$scale)
    fc <- b*gamma((1/c)+1)
    list(cycle=fc, freq=1/fc)
}

coxFitFn <- function(x) {
    cox_reg <- coxph(x ~ 1)
    base.cox <- basehaz(cox_reg)
    ##fc <- max(base.cox[,2])/max(base.cox[,1])
    fc <- ifelse(length(base.cox[,1])==1, max(base.cox[,2])/max(base.cox[,1]), mean(c(base.cox[length(base.cox[,1])-1,2],base.cox[length(base.cox[,1]),2])) / mean(c(base.cox[length(base.cox[,1])-1,1],base.cox[length(base.cox[,1]),1])))
    list(cycle=fc, freq=1/fc)
}
