#######################################
######## fire cycle estimation survival function
#######################################

expFitFnc <- function(x, ind = NULL) {        ### negative exponential fitting
    if(!is.null(ind)) { ### a condition that is met when using 'boot' function
        x <- Surv(x[ind,"time"], x[ind,"status"])
    }
    surv.exp <- survreg(x ~ 1, dist="exponential")
    b <- exp(surv.exp$coefficients)
    return(b)
}

weibFitFnc <- function(x, ind = NULL) {        ### weibull exponential fitting
    if(!is.null(ind)) { ### a condition that is met when using 'boot' function
        x <- Surv(x[ind,"time"], x[ind,"status"])
    }
    surv.weib <- survreg(x ~ 1, dist="weibull")
    b <- exp(surv.weib$coefficients)
    c <- 1/(surv.weib$scale)
    fc <- b*gamma((1/c)+1)
    return(fc)
}

coxFitFnc <- function(x, ind = NULL) {
    if(!is.null(ind)) { ### a condition that is met when using 'boot' function
        x <- Surv(x[ind,"time"], x[ind,"status"])
    }
    cox_reg <- coxph(x ~ 1)
    base.cox <- basehaz(cox_reg)
    maxHaz <- max(base.cox[,1])
    index <- which(base.cox[,1] == maxHaz)
    fc <- mean(base.cox[index,"time"])/maxHaz
    #     base.cox <- distinct(base.cox, hazard)
    #     fc <- max(base.cox[,2])/max(base.cox[,1])
    # #     fc <- ifelse(length(base.cox[,1])==1, max(base.cox[,2]) / max(base.cox[,1]),
    #                  mean(c(base.cox[length(base.cox[,1])-1,2], base.cox[length(base.cox[,1]),2])) / mean(c(base.cox[length(base.cox[,1])-1,1],base.cox[length(base.cox[,1]),1])))
    return(fc)
    #return(list(cycle=fc, freq=1/fc))
}
