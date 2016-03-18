# Fire cycle estimation functions (and other scripts)
Dominic Cyr  

Updated on Mar 18 2016


-----------


The present repository contains a script that can be sourced to provide with functions that implement three different methods for estimating fire cycle from time since fire data: [_fcEstSurvFnc.R_][6].

Other scripts can be found as well and are not described here. In-file comments should provide a good enough description. The [simulation and data processing pipeline description][3] provides additionnal information about how those scripts interacts with others, as well as examples of how to use them, what inputs they need, etc.

All three functions contained in [_fcEstSurvFnc.R_][6], `expFitFnc`, `weibFitFnc` and `coxFitFnc`, are very simple to use, all work the same way and return a single numerical estimation of the length of the fire cycle. Using the case study presented in our [paper][1], here are examples of how to use them.


```r
require(survival)
source("./fcEstSurvFnc.R")
# loading time since fire (tsf) data
tsf <- read.csv("../data/tsf.csv")
tsf <- tsf[complete.cases(tsf),]
head(tsf)
```

```
##   ID      lat      long samplingDate fireDate status weight
## 1  1 49.21298 -68.44260         2004     1759      0   1.60
## 2  2 49.33618 -68.76995         2004     1855      0   1.60
## 3  3 49.72644 -68.48300         2004     1923      1   1.60
## 4  4 49.94262 -68.74809         2004     1645      0   0.31
## 5  5 49.95407 -68.76481         2004     1645      1   0.20
## 6  6 50.06534 -68.76459         2004     1735      1   0.80
```

```r
# First create a 'Surv' object
x <- Surv(tsf$samplingDate - tsf$fireDate, ### time since fire in years
     tsf$status) ## 0 (censored) / 1 (uncensored)

# functions directly return fire cycle
fc <- coxFitFnc(x)
fc
```

```
## [1] 229.0309
```

It is also possible (and recommended) to produce confidence intervals and to weigh observations in cases of unbalanced samples. Here we used the bias-corrected and accelerated (_BCa_) bootstrap, by Efron ([1987][7]), which adjusts for both bias and skewness in the bootstrap distribution. 


```r
# weight variable
w <- tsf$weight
# bootstrap resampling effort
n <- 10000
require(boot)
# create 'boot' object, note the 'weight' variable
coxBoot <- boot(as.matrix(x), statistic = coxFitFnc, R = n, weights = w, sim = "ordinary")
# weighted and bias corrected bootstrap
coxCI <- boot.ci(coxBoot, type = "bca") 
# 95% confidence interval
coxCI
```

```
## BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
## Based on 10000 bootstrap replicates
## 
## CALL : 
## boot.ci(boot.out = coxBoot, type = "bca")
## 
## Intervals : 
## Level       BCa          
## 95%   (160.6, 386.3 )  
## Calculations and Intervals on Original Scale
## Some BCa intervals may be unstable
```

Here you can note that the returned 'bootci' object mention that "Some BCa intervals may be unstable". That indicates that the size of the original sample is relatively small and induces a relativaly large variability of extreme quantiles. That may indicates that decreasing the confidence level from the default of 95% to 90% may provide with more reliable coverage performance.


```r
# weighted and bias corrected bootstrap
coxCI <- boot.ci(coxBoot, type = "bca", conf = 0.9) 
# 90% confidence interval
coxCI
```

```
## BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
## Based on 10000 bootstrap replicates
## 
## CALL : 
## boot.ci(boot.out = coxBoot, conf = 0.9, type = "bca")
## 
## Intervals : 
## Level       BCa          
## 90%   (174.3, 374.9 )  
## Calculations and Intervals on Original Scale
```

```r
# Extracting numerical values
coxEstimate <- data.frame(estimate = round(coxBoot$t0, 1), ll = NA, ul = NA)
coxEstimate[,c("ll", "ul")] <- round(t(sapply(coxCI[-(1:3)],function(x) tail(c(x),2))),1)
coxEstimate
```

```
##   estimate    ll    ul
## 1      229 174.3 374.9
```



Don't hesitate to contact me should you have questions, comments or suggestions for improvement or additionnal functionalities.

[Dominic Cyr][5]

[2]: https://github.com/dcyr/survFire/scripts
[3]: https://github.com/dcyr/survFire/blob/master/pipeline.md
[5]: http://dominiccyr.ca
[6]: https://github.com/dcyr/survFire/blob/master/scripts/fcEstSurvFnc.R
[7]: http://www.tandfonline.com/doi/abs/10.1080/01621459.1987.10478410
