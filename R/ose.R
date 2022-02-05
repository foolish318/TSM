## Estimate DTSM with observed shifting endpoint (OSE)
## save estimation results in results/ose.RData

library(KFAS)
source("R/dtsm_fns.R")
source("R/data_fns.R")
source("R/util_fns.R")

filename <- "results/ose.RData"

cat("# Estimating model\n")
cat("#", filename, "\n")
## load data
data <- loadData()
data$istar <- (data$pistar.ptr + data$rstar.realtime)/100
yield.cols <- attr(data, "yield.cols")
mats <- as.numeric(substr(yield.cols, 2, 5))
dt <- 1/4   # quarterly frequency
## ## post 1985
## data <- data[data$yyyymm>198500,]
## yields
Y <- as.matrix(data[yield.cols])/100
W <- eigen(cov(Y))$vectors
N <- 3 # number of factors
WN <- t(W[,1:N])   # N x J
J <- ncol(Y)
nobs <- nrow(Y)
cP <- Y %*% t(WN)  # nobs x N
est.sample <- data$yyyymm<200800 # estimation sample - leave out ZLB period
set.seed(616)

## FE model for starting values
## stationary, observed cP -- to get Q-dynamics for starting values
cat("# 1) Estimate JSZ (FE) model to get good starting values\n")
mod.jsz <- estimate_jsz(Y, est.sample, WN, mats, dt)

## OSE model
cat("# 2) Estimate OSE model\n")
cat("# Get starting values\n")
pars.start <- startval_ose(Y[est.sample,], data$istar[est.sample], mod.jsz, WN, mats, dt)
cat("P-eigenvalues:\n")
print(abs(eigen(pars.start$Phi)$values))
cat("Q-eigenvalues:\n")
print(pars.start$lamQ)
cat("gamma:\n")
print(pars.start$gamma)
cat("Pbar:\n")
print(pars.start$Pbar)
theta.start <- pars2theta_ose(pars.start)
pars <- theta2pars_ose(theta.start, WN, mats, dt)
theta.check <- pars2theta_ose(pars)
stopifnot(all.equal(theta.start, theta.check))
cat("# Numerical optimization of OSE LLK\n")
llk.start <- -obj_ose(theta.start, data$istar[est.sample], Y[est.sample,], WN, mats, dt)
cat("LLK at starting values:", llk.start, "\n")
res <- llk_ose(Y[est.sample,], data$istar[est.sample], WN, pars$kinfQ, pars$lamQ, pars$p, pars$a, pars$Sigma, pars$sigma.tau, mats, dt)
stopifnot(all.equal(llk.start, sum(res$llk)))
theta <- get_optim(theta.start, obj_ose, data$istar[est.sample], Y[est.sample,], WN, mats, dt, trace=1)
pars <- theta2pars_ose(theta, WN, mats, dt)
optimal.llk <- -obj_ose(theta, data$istar[est.sample], Y[est.sample,], WN, mats, dt)
cat("# LLK at optimal values:", optimal.llk, "\n")
rval <- llk_ose(Y[est.sample,], data$istar[est.sample], WN, pars$kinfQ, pars$lamQ, pars$p, pars$a, pars$Sigma, pars$sigma.tau, mats, dt)
stopifnot(all.equal(sum(rval$llk), optimal.llk))
## full-sample results (including ZLB)
rval <- llk_ose(Y, data$istar, WN, pars$kinfQ, pars$lamQ, pars$p, pars$a, pars$Sigma, pars$sigma.tau, mats, dt)
cat("Full sample:", sum(rval$llk), "\n")
mod <- c(pars, rval)

mod <- within(mod, {
    muQ <- K0Q.cP
    PhiQ <- K1Q.cP + diag(N)
    cPstar <- rep(1, nobs) %*% rbind(Pbar) + istar %*% rbind(gamma)
    Ystar <- rep(1, nobs) %*% AcP + cPstar %*% BcP
    Yhat <- Yhat <- rep(1, nrow(cP)) %*% AcP + cP %*% BcP
    Ytilde <- Y - Ystar
    Z <- cbind(istar, cP)
    mu.Z <- cbind(c(0, (diag(N) - Phi)%*%Pbar))
    Phi.Z <- cbind(c(1, (diag(N) - Phi)%*%gamma), rbind(0, Phi))
    Sigma.Z <- matrix(NA, N+1, N+1)
    Sigma.Z[1,1] <- sigma.tau^2
    Sigma.Z[-1,1] <- gamma * sigma.tau^2
    Sigma.Z[1,-1] <- gamma * sigma.tau^2
    Sigma.Z[-1,-1] <- Sigma
})

stopifnot(check_restrictions(mod$Pbar, mod$gamma, mod$rho0.cP, mod$rho1.cP))
persistence(mod)

save(data, mod, mod.jsz, est.sample, Y, WN, N, J, nobs, mats, yield.cols, dt, file = filename, compress = TRUE)

