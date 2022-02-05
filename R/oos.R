## Out-of-sample forecasts with DTSM (with observed i*)
## Table 5

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
yieldvars <- c("y10")

## forecasts + errors: dates x horizons x variables
yhat_jsz <- array(NA, c(T-H-t0+1, length(h), length(yieldvars)))
yhat_se <- array(NA, c(T-H-t0+1, length(h), length(yieldvars)))
e_rw <- array(NA, c(T-H-t0+1, length(h), length(yieldvars)))
e_jsz <- array(NA, c(T-H-t0+1, length(h), length(yieldvars)))
e_se <- array(NA, c(T-H-t0+1, length(h), length(yieldvars)))

## Bluechip forecasts of 10y yield
annAvg <- function(x) {
    ## for quarterly forecasts or future yields, calculate annual averages
    stopifnot(length(x)==20)
    vapply(1:5, function(j) mean(x[((j-1)*4+1):(j*4)]), numeric(1))  ## could do this more elegantly with rollmean()
}

bc <- loadBlueChip()
bc <- bc[bc$year5<=2017,]
## forecast errors: dates x horizons x models
e_bc <- array(NA, c(nrow(bc), 5, 4))
rownames(e_bc) <- bc$yyyymm
data$yyyyq <- data$year*10 + (data$month-1) %/% 3 + 1 # useful to match blue chip forecast quarters
cat("# First Bluechip forecast date:", bc$yyyymm[1], "\n")
cat("# Last Bluechip forecast date:", tail(bc$yyyymm, 1), "\n")
cat("# Number of Bluechip forecast dates:", nrow(bc), "\n")
data$bc_date <- FALSE
data$bc_flag <- FALSE
for (i in 1:nrow(bc)) {
    t <- max(which(data$yyyymm < bc$yyyymm[i]))
    data$bc_flag[t] <- TRUE
    data$bc_date[t] <- bc$yyyymm[i]
}
t2 <- max(T-H, max(which(data$bc_flag==TRUE))) ## last forecast date
cat("# Last forecast date (Blue Chip or quarterly):", data$yyyymm[t2], "\n")

llk.jsz.prev <- -Inf
llk.se.prev <- -Inf
start <- Sys.time()
for (t in t0:t2) {
    cat("##################################################\n")
    cat("## Forecast date:", data$yyyymm[t], "\n")
    cat("Duration so far:", format(Sys.time()-start), "\n")

    ## Blue Chip
    if (data$bc_flag[t]) {
            cat("## Bluechip forecast date:", data$bc_date[t], "\n")
            ## determine forecast horizons
            ## - 20 horizons, from first quarter in year1 to last quarter in year5
            ## - these will be averaged over years
            i <- which(bc$yyyymm == data$bc_date[t])
            h1 <- min(which(data$year == bc$year1[i])) - t
            h2 <- h1+20-1
            cat("Quarterly forecast horizons from:", data$yyyymm[t+h1], "to", data$yyyymm[t+h2], "\n")
            cat("Bluechip forecasts from:", bc$year1[i], "to", bc$year5[i], "\n")
            ## future yields
            yfut <- annAvg(data$y10[t+(h1:h2)])
            ## Blue Chip
            e_bc[i, ,1] <- yfut - unlist(bc[i, paste0("f", 1:5)])
            ## Random walk
            e_bc[i, ,2] <- yfut - data$y10[t]
    }

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
    yhat_bc <- numeric(20)
    for (j in 1:H) {
        EcP <- drop(mod.jsz$mu + mod.jsz$Phi %*% EcP)
        EY <- mod.jsz$AcP + EcP %*% mod.jsz$BcP
        if (t <= T-H && any(h==j)) {
            yhat_jsz[t-t0+1, which(h==j), ] <- EY[match(yieldvars, yield.cols)]
            e_jsz[t-t0+1, which(h==j), ] <- as.matrix(data[t+j, yieldvars]/100) - EY[match(yieldvars, yield.cols)]
        }
        if (data$bc_flag[t] && (j>=h1 & j<=h2))
            yhat_bc[j-h1+1] <- 100*EY[mats==10]
    }
    if (data$bc_flag[t])
        e_bc[i, ,3] <- yfut - annAvg(yhat_bc)

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
    yhat_bc <- numeric(20)
    for (j in 1:H) {
        EPtilde <- drop(mod.se$Phi %*% EPtilde)
        EcP <- mod.se$Pbar + data$istar[t] * mod.se$gamma + EPtilde
        EY <- mod.se$AcP + EcP %*% mod.se$BcP
        if (t <= T-H && any(h==j)) {
            yhat_se[t-t0+1, which(h==j), ] <- EY[match(yieldvars, yield.cols)]
            e_se[t-t0+1, which(h==j), ] <- as.matrix(data[t+j, yieldvars]/100) - EY[match(yieldvars, yield.cols)]
        }
        if (data$bc_flag[t] && (j>=h1 & j<=h2))
            yhat_bc[j-h1+1] <- 100*EY[mats==10]
    }
    if (data$bc_flag[t])
        e_bc[i, ,4] <- yfut - annAvg(yhat_bc)

}
cat("This took:", format(Sys.time() - start), "\n")

rmse <- function(e)
    apply(e, 2, function(x) sqrt(mean(x^2)))
mae <- function(e)
    apply(e, 2, function(x) mean(abs(x)))
for (j in seq_along(yieldvars)) {
    cat("##", yieldvars[j], "- quarterly\n")
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
    if (yieldvars[j] == "y10")
        print(xtable(tbl_rmse, digits=2),
              include.rownames=TRUE, include.colnames=FALSE, only.contents=TRUE,
              sanitize.text.function=function(x){x}, hline.after=3,
              file = "tables/oos_10y_rmse.tex")
}
## Bluechip
cat("##", yieldvars[j], "- Blue Chip\n")
tbl_rmse <- rbind(apply(e_bc, c(3,2), function(x) sqrt(mean(x^2))), NA, NA, NA)
rownames(tbl_rmse) <- c("Blue Chip (\\emph{BC})", "Random walk (\\emph{RW})", "Fixed endpoint (\\emph{FE})", "Observed shifting endpoint (\\emph{OSE})", "$p$-value: \\emph{OSE} $\\geq$ \\emph{BC}", "$p$-value: \\emph{OSE} $\\geq$ \\emph{RW}", "$p$-value: \\emph{OSE} $\\geq$ \\emph{FE}")
tbl_mae <- rbind(apply(e_bc, c(3,2), function(x) mean(abs(x))), NA, NA, NA)
rownames(tbl_mae) <- rownames(tbl_rmse)
colnames(tbl_rmse) <- 1:5; colnames(tbl_mae) <- 1:5
## Diebold-Mariano
for (j in 1:5) {
    ## OSE vs BC
    tbl_rmse[5, j] <- myDMtest(e_bc[,j,4], e_bc[,j,1], alternative="less", h=j*4, power=2)$p.value
    tbl_mae[5, j] <- myDMtest(e_bc[,j,4], e_bc[,j,1], alternative="less", h=j*4, power=1)$p.value
    ## OSE vs RW
    tbl_rmse[6, j] <- myDMtest(e_bc[,j,4], e_bc[,j,2], alternative="less", h=j*4, power=2)$p.value
    tbl_mae[6, j] <- myDMtest(e_bc[,j,4], e_bc[,j,2], alternative="less", h=j*4, power=1)$p.value
    ## OSE vs FE
    tbl_rmse[7, j] <- myDMtest(e_bc[,j,4], e_bc[,j,3], alternative="less", h=j*4, power=2)$p.value
    tbl_mae[7, j] <- myDMtest(e_bc[,j,4], e_bc[,j,3], alternative="less", h=j*4, power=1)$p.value
}
cat("RMSE:\n")
print(round(tbl_rmse, 2))
## cat("MAE:\n")
## print(round(tbl_mae, 2))
print(xtable(tbl_rmse, digits=2),
      include.rownames=TRUE, include.colnames=FALSE, only.contents=TRUE,
      sanitize.text.function=function(x){x}, hline.after=4,
      file = "tables/oos_10y_bluechip_rmse.tex")

