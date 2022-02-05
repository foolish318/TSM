## ESE model: analyze and combine results of separate MCMC chains

graphics.off()
source("R/dtsm_fns.R") # jsz.loadings()
source("R/util_fns.R") # makeCovMat()

filenames <- paste0("results/", c("ese_616.RData", "ese_1.RData", "ese_2.RData", "ese_3.RData", "ese_4.RData"))
filename_all <- "results/ese.RData" # for combined results

cat("Loading and combining MCMC chains...\n")
cat(filenames[1], "\n")
load(filenames[1])
M <- length(mcmc.ind)
N <- ncol(lamQ)
gamma.ind <- setdiff(1:3, which(apply(gamma, 2, sd)==0))[1:2]
Pbar.ind <- gamma.ind

## L parallel chains
L <- length(filenames)
parnames <- c("100*kinfQ", "lamQ_1", "gamma_1", "gamma_2", "Pbar_1", "Pbar_2", "100*sig.tau", "100*sig.e", "eig(Phi)", paste("100*Sigma", 1:3))
theta <- array(NA, c(L, M, 12))
## colnames(theta) <- parnames
theta[1,,] <- cbind(100*kinfQ, lamQ[,1], gamma[,gamma.ind], Pbar[,Pbar.ind], 100*sqrt(sigtau2), 100*sqrt(sige2), eigenval, 100*sqrt(Sigma[,c(1,4,6)]))
for (j in 2:L) {
    cat(filenames[j], "\n")
    load(filenames[j])
    theta[j,,] <- cbind(100*kinfQ, lamQ[,1], gamma[,gamma.ind], Pbar[,Pbar.ind], 100*sqrt(sigtau2), 100*sqrt(sige2), eigenval, 100*sqrt(Sigma[,c(1,4,6)]))
}

## plot all cumulative ergodic means
dev.new()
par(mfrow = c(2,ceiling(length(parnames)/2)),
    mar=c(3,2,1,1))
for (i in seq_along(parnames)) {
    em <- matrix(NA, L, M)
    for (j in 1:L)
        em[j,] <- cumsum(theta[j,,i])/1:M
    plot(em[1,], type="l", ylim=range(em), main=parnames[i])
    for (j in 2:L)
        lines(em[j,])
}

## flatten theta - multiple chains to one combined chain
theta <- matrix(theta, L*M, 12)

## summary statistics of MCMC posterior distribution
tbl <- matrix(NA, length(parnames), 5)
rownames(tbl) <- parnames
colnames(tbl) <- c("MCMC mean", "median", "2.5%", "97.5%", "SD")
tbl[,1] <- colMeans(theta, na.rm=TRUE)   # posterior means
tbl[,2:4] <- t(apply(theta, 2, quantile, c(.5, .025, .975), na.rm=TRUE))
tbl[,5] <- apply(theta, 2, sd)
print(round(tbl,4))

## flatten chains for all parameters
nms <- c("AcP", "BcP", "gamma", "Pbar", "Phi", "tau", "cP", "sigtau2", "sige2", "kinfQ", "lamQ", "Sigma", "Sigma.tilde", "eigenval", "alpha")
nms <- unique(nms)
load(filenames[1])
for (i in seq_along(nms))
    assign(paste0(nms[i], "_all"), get(nms[i]))
for (j in 2:L) {
    load(filenames[j])
    for (i in seq_along(nms))
        assign(paste0(nms[i], "_all"), abind::abind(get(paste0(nms[i], "_all")), get(nms[i]), along=1))
}
for (i in seq_along(nms))
    assign(nms[i], get(paste0(nms[i], "_all")))
rm(list = paste0(nms, "_all"))
M <- nrow(lamQ)

cat("Mean acceptance probabilities:\n")
print(colMeans(alpha))

cat("Writing", filename_all, "...\n")

## save combined results
save(data, zlb.ind, Y, WN, mats, yield.cols, dt, prior, kinfQ, lamQ, gamma, sigtau2, sige2, eigenval, Pbar, Phi, Sigma, Sigma.tilde, AcP, BcP, alpha, tau, cP, pars2theta, makegamma, makePbar, file = filename_all, compress=TRUE)
