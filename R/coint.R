## Table 1 (and Online Appendix Table 3)
## Cointegration regressions, persistence of residuals, tests for rank of cointegration

library(dynlm)
library(urca)
library(xtable)
source("R/data_fns.R")
source("R/util_fns.R")
source("R/ur_tests.R")

data <- loadData()
data$istar <- data$pistar.ptr + data$rstar.realtime

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
    ## if no command line arguments: analyze ten-year yield
    ## Table 1: analysis of 10y yield
    filename <- "tables/coint.tex"
    yield.var <- "y10"
} else if (length(args)==1 && args[1]=="level") {
    ## else if "level" is command line argument: analyze level of yield curve
    ## Online Appendix Table 2: analysis of level (PC1)
    filename <- "tables/coint_level.tex"
    yield.var <- "level"
} else {
    print(args)
    stop("invalid command line argument")
}

rstar.nms <- c("rstar.filt", "rstar.realtime", "rr.ewma")
rstar.desc <- c("filtered", "real-time", "mov.~avg.")

nms <- c("constant", "",
         "$\\pi_t^\\ast$", "",
         "$r_t^\\ast$", "",
         "$i_t^\\ast$", "",
         "$R^2$",
         "Memo: $r^\\ast$",
         "SD", "$\\hat{\\rho}$", "Half-life", "ADF", "PP", "LFST",
         "Johansen $r=0$", "Johansen $r=1$",
         "ECM $\\hat{\\alpha}$", "")
tbl <- data.frame(matrix(NA, length(nms), 4+length(rstar.nms)))
tbl[,1] <- sprintf("%-15s", nms)
tbl[, 2] <- ""
mod <- lm(y10 ~ 1, data)
v <- drop(sandwich::NeweyWest(mod, lag=6, prewhite=FALSE))
tbl[1, 2] <- sprintf("%4.2f", mod$coef)
tbl[2, 2] <- sprintf("(%4.2f)", sqrt(v))
tbl[11:16,2] <- persistence_stats(data$y10)

## DOLS
p <- 4 # both DOLS and Johansen test
## create first differences for DOLS
for (nm in c("pistar.ptr", rstar.nms, "istar"))
    data[paste0("d.", nm)] <- c(NA, diff(data[[nm]]))
col <- 3
for (j in 0:(length(rstar.nms)+1)) {
    if (j==0) {
        ## pi-star only
        fmla <- get(yield.var) ~ pistar.ptr + L(d.pistar.ptr, -p:p)
        regrows <- c(1,3)
        tbl[c(5,6,7,8,10), col] <- ""
        ind <- 1:2 # which coefficients
        X <- data[, c(yield.var, "pistar.ptr")]
        ecmod.fmla <- dy ~ L(resid, 1) + L(dy, 1:p) + L(d.pistar.ptr, 1:p)
    } else if (j>0 & j<=length(rstar.nms)) {
        ## with r-star
        regrows <- c(1,3,5)
        rvar <- rstar.nms[j]
        ## fmla <- get(yield.var) ~ pistar.ptr + get(rvar) # OLS
        fmla <- formula(paste0(yield.var, " ~ pistar.ptr + ", rvar, " + L(d.pistar.ptr, -p:p) + L(d.", rvar, ", -p:p)"))
        tbl[c(7,8), col] <- ""
        tbl[10, col] <- rstar.desc[j]
        ind <- 1:3
        X <- data[, c(yield.var, "pistar.ptr", rvar)]
        ecmod.fmla <- formula(paste0("dy ~ L(resid, 1) + L(dy, 1:p) + L(d.pistar.ptr, 1:p) + L(d.", rvar, ", 1:p)"))
    } else {
        ## with i-star
        regrows <- c(1,7)
        fmla <- get(yield.var) ~ istar + L(d.istar, -p:p)
        tbl[10, col] <- "real-time"
        tbl[c(3,4,5,6),col] <- ""
        ind <- 1:2
        X <- data[, c(yield.var, "istar")]
        ecmod.fmla <- dy ~ L(resid, 1) + L(dy, 1:p) + L(d.istar, 1:p)
    }
    ## DOLS
    mod <- dynlm(fmla, data=ts(data))
    V <- sandwich::NeweyWest(mod, lag=6, prewhite=FALSE)
    tbl[regrows, col] <- sprintf("%4.2f", mod$coef[ind])
    tbl[regrows+1, col] <- sprintf("(%4.2f)", sqrt(diag(V))[ind])
    tbl[9, col] <- sprintf("%4.2f", summary(mod)$r.squared)
    ## persistence statistics
    if (j == 0) {
        data$resid <- data[[yield.var]] - mod$coef[1] - mod$coef[2]*data$pistar.ptr
    } else if (j>0 & j<=length(rstar.nms)) {
        data$resid <- data[[yield.var]] - mod$coef[1] - mod$coef[2]*data$pistar.ptr - mod$coef[3]*data[[rvar]]
    } else {
        data$resid <- data[[yield.var]] - mod$coef[1] - mod$coef[2]*data$istar
    }
    tbl[11:16, col] <- persistence_stats(data$resid, nvar=ifelse(j==0, 2, 3))
    ## Johansen test for cointegration rank
    ## rval <- vars::VARselect(X)
    ## cat("Selection criterion for Johansen VAR:\n")
    ## print(rval$selection)
    ## cat("using", p, "lags\n")
    vecmod <- ca.jo(X, type="trace", ecdet="const", K=p)
    for (r in 0:1) {
        ind <- length(vecmod@teststat) - r
        result <- sprintf("%4.2f", vecmod@teststat[ind])
        sig <- vecmod@teststat[ind] > vecmod@cval[ind,]
        if (any(sig))
            result <- paste(c(result, rep("*", sum(sig))), collapse="")
        tbl[17+r, col] <- result
    }
    ## Error-correction model
    data$dy <- c(NA, diff(data[[yield.var]]))
    ecmod <- dynlm::dynlm(ecmod.fmla, data=ts(data))
    b <- ecmod$coef
    V <- sandwich::vcovHC(ecmod, "HC0") # robust SEs
    SEs <- sqrt(diag(V))
    tstats <- b/SEs
    pvals <- 2*pt(abs(tstats), df=ecmod$df.residual, lower.tail=FALSE)
    tbl[19, col] <- sprintf("%4.2f", b[2])
    tbl[20, col] <- sprintf("(%4.2f)", SEs[2])
    col <- col+1
}
print(tbl)

sink(filename)
print(xtable(tbl),
      include.rownames=FALSE, include.colnames=FALSE, only.contents=TRUE,
      sanitize.text.function=function(x){x}, hline.after=c(10, 16, 18))
sink()
