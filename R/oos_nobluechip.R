## Out-of-sample forecasts with DTSM (with observed i*)
##  if Blue Chip data is not available
## Table 5, top panel only

library(KFAS)
library(xtable)
source("R/dtsm_fns.R")
source("R/data_fns.R")
source("R/util_fns.R")
set.seed(616)

## load data
data <- loadData()
data$istar <- (data$pistar.ptr + data$rstar.realtime)/100
yield.cols <- attr(data, "yield.cols")
mats <- as.numeric(substr(yield.cols, 2, 5))
dt <- 1/4   # quarterly frequency
T <- nrow(data)
t0 <- 20  ## with five years of data
h <- c(4,10,20,30,40) # forecast horizons
H <- max(h)
t1 <- max(which(data$yyyymm<200800)) ## last date to estimate models
cat("# First forecast date:", data$yyyymm[t0], "\n")
cat("# Last estimation date:", data$yyyymm[t1], "\n")
cat("# Last forecast date (quarterly, models only):", data$yyyymm[T-H], "\n")
cat("# Number of quarterly forecast dates:", T-H-t0+1, "\n")
N <- 3 # number of factors
J <- length(yield.cols)
yieldvars <- "y10"

## forecasts + errors: dates x horizons x variables
yhat_jsz <- array(NA, c(T-H-t0+1, length(h), length(yieldvars)))
yhat_se <- array(NA, c(T-H-t0+1, length(h), length(yieldvars)))
e_rw <- array(NA, c(T-H-t0+1, length(h), length(yieldvars)))
e_jsz <- array(NA, c(T-H-t0+1, length(h), length(yieldvars)))
e_se <- array(NA, c(T-H-t0+1, length(h), length(yieldvars)))

t2 <- T-H
cat("# Last forecast date (quarterly):", data$yyyymm[t2], "\n")

llk.jsz.prev <- -Inf
llk.se.prev <- -Inf
start <- Sys.time()
for (t in t0:t2) {
    cat("##################################################\n")
    cat("## Forecast date:", data$yyyymm[t], "\n")
    cat("Duration so far:", format(Sys.time()-start), "\n")

    if (t <= T-H) {
        cat("## Quarterly forecasts\n")
        ## Random walk forecasts
        for (j in seq_along(yieldvars)) {
            y <- data[[yieldvars[j]]]/100
            e_rw[t-t0+1, , j] <- y[t+h] - y[t]
        }
    }

    Y <- as.matrix(data[1:t, yield.cols])/100
    W <- eigen(cov(Y))$vectors
    WN <- t(W[,1:N])   # N x J
    WN <- WN * sign(WN[,1])
    cP <- Y %*% t(WN)  # nobs x N
    istar <- data$istar[1:t]

    ## JSZ model
    if (data$yyyymm[t] < 200800) {
        cat("## 1) Estimate JSZ model -- stationary, observed cP\n")
        ## starting values
        if (t == t0) {
            lm <- ar.ols(cP, aic=FALSE, order.max=1, intercept=TRUE, demean=FALSE)
            L <- t(chol(lm$var.pred))
            pars.start <- getStartingValuesForMLE(L, Y, WN, mats, dt)
            theta.start <- pars2theta_jsz(pars.start)
        } else {
            theta.start <- pars2theta_jsz(mod.jsz)
        }
        cat("LLK at starting values:", -obj_jsz(theta.start, Y, WN, mats, dt), "\n")
        theta <- get_optim(theta.start, obj_jsz, Y, WN, mats, dt)
        pars <- theta2pars_jsz(theta, N)
        llk.jsz <- -obj_jsz(theta, Y, WN, mats, dt)
        cat("LLK at optimal values:", llk.jsz, "\n")
        if (llk.jsz < llk.jsz.prev)
            cat("# Note: likelihood not increased with additional data !!\n")
        llk.jsz.prev <- llk.jsz
        res.llk <- jsz.llk(Y, WN, K1Q.X=diag(pars$lamQ-1), Sigma.cP=pars$Omega.cP, mats=mats, dt=dt)
        mod.jsz <- within(c(pars, res.llk), {
            mu <- K0P.cP
            Phi <- K1P.cP + diag(N)
            PhiQ <- K1Q.cP + diag(N)
            Yhat <- rep(1, nrow(cP)) %*% AcP + cP %*% BcP
        })
        persistence(mod.jsz)
    }
    ## calculate forecasts
    EcP <- cP[t,] # E_t cP_t+j   j in 1:H
    for (j in 1:H) {
        EcP <- drop(mod.jsz$mu + mod.jsz$Phi %*% EcP)
        EY <- mod.jsz$AcP + EcP %*% mod.jsz$BcP
        if (t <= T-H && any(h==j)) {
            yhat_jsz[t-t0+1, which(h==j), ] <- EY[match(yieldvars, yield.cols)]
            e_jsz[t-t0+1, which(h==j), ] <- as.matrix(data[t+j, yieldvars]/100) - EY[match(yieldvars, yield.cols)]
        }
    }

    ## DTSM with shifting endpoint
    if (data$yyyymm[t] < 200800) {
        cat("## 2) Estimate DTSM with observed shifting endpoint\n")
        cat("## a) optimize with fresh starting values\n")
        pars.start <- startval_ose(Y, istar, mod.jsz, WN, mats, dt)
        theta.start <- pars2theta_ose(pars.start)
        theta <- get_optim(theta.start, obj_ose, istar, Y, WN, mats, dt)
        cat("LLK =", llk.se <- -obj_ose(theta, istar, Y, WN, mats, dt), "\n")
        if (t>t0) {
            cat("## b) also optimize with values from previous iteration\n")
            theta.start <- pars2theta_ose(mod.se)
            theta2 <- get_optim(theta.start, obj_ose, istar, Y, WN, mats, dt)
            cat("LLK =", llk.se.2 <- -obj_ose(theta2, istar, Y, WN, mats, dt), "\n")
            if (llk.se.2 > llk.se) {
                theta <- theta2
                llk.se <- llk.se.2
            }
        }
        pars <- theta2pars_ose(theta, WN, mats, dt)
        cat("LLK at optimal values:", llk.se, "\n")
        if (llk.se < llk.se.prev)
            cat("# Note: likelihood not increased with additional data !!\n")
        llk.se.prev <- llk.se
        rval <- llk_ose(Y, istar, WN, pars$kinfQ, pars$lamQ, pars$p, pars$a, pars$Sigma, pars$sigma.tau, mats, dt)
        stopifnot(all.equal(sum(rval$llk), llk.se))
        mod.se <- within(c(pars, rval), {
            istar <- data$istar[1:t]
            cP <- cP
            PhiQ <- K1Q.cP + diag(N)
            Yhat <- rep(1, t) %*% AcP + cP %*% BcP
            ystar <- AcP + (Pbar + istar[t] * gamma) %*% BcP # Y^\ast_t
        })
        persistence(mod.se)
    }
    ## calculate forecasts
    EPtilde <- cP[t,] - mod.se$Pbar - data$istar[t] * mod.se$gamma
    for (j in 1:H) {
        EPtilde <- drop(mod.se$Phi %*% EPtilde)
        EcP <- mod.se$Pbar + data$istar[t] * mod.se$gamma + EPtilde
        EY <- mod.se$AcP + EcP %*% mod.se$BcP
        if (t <= T-H && any(h==j)) {
            yhat_se[t-t0+1, which(h==j), ] <- EY[match(yieldvars, yield.cols)]
            e_se[t-t0+1, which(h==j), ] <- as.matrix(data[t+j, yieldvars]/100) - EY[match(yieldvars, yield.cols)]
        }
    }
}
cat("This took:", format(Sys.time() - start), "\n")

rmse <- function(e)
    apply(e, 2, function(x) sqrt(mean(x^2)))
mae <- function(e)
    apply(e, 2, function(x) mean(abs(x)))

for (j in seq_along(yieldvars)) {
    cat("##", yieldvars[j], "\n")
    tbl_rmse <- 100*rbind(rmse(e_rw[,,j]), rmse(e_jsz[,,j]), rmse(e_se[,,j]), NA, NA)
    tbl_mae <- 100*rbind(mae(e_rw[,,j]), mae(e_jsz[,,j]), mae(e_se[,,j]), NA, NA)
    for (i in seq_along(h)) {
        ## SE vs RW
        tbl_rmse[4, i] <- myDMtest(e_se[,i,j], e_rw[,i,j], alternative="less", h=h[i], power=2)$p.value
        tbl_mae[4, i] <- myDMtest(e_se[,i,j], e_rw[,i,j], alternative="less", h=h[i], power=1)$p.value
        ## SE vs JSZ
        tbl_rmse[5, i] <- myDMtest(e_se[,i,j], e_jsz[,i,j], alternative="less", h=h[i], power=2)$p.value
        tbl_mae[5, i] <- myDMtest(e_se[,i,j], e_jsz[,i,j], alternative="less", h=h[i], power=1)$p.value
    }
    colnames(tbl_rmse) <- h; colnames(tbl_mae) <- h
    rownames(tbl_rmse) <- c("Random walk (\\emph{RW})", "Fixed endpoint (\\emph{FE})", "Observed shifting endpoint (\\emph{OSE})", "$p$-value: \\emph{OSE} $\\geq$ \\emph{RW}", "$p$-value: \\emph{OSE} $\\geq$ \\emph{FE}")
    rownames(tbl_mae) <- rownames(tbl_rmse)
    cat("RMSE:\n")
    print(round(tbl_rmse, 2))
    ## cat("MAE:\n")
    ## print(round(tbl_mae, 2))
}

