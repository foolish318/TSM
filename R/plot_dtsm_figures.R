## Plot figures based on DTSM estimates
## Figures 4, 5, 7

graphics.off()
library(dplyr)
source("R/util_fns.R")
source("R/dtsm_fns.R")

## OSE model
load("results/ose.RData")
T <- nrow(data)
J <- ncol(Y)
## TP
mod <- within(mod, {
    Yhat <- rep(1, T) %*% AcP + cP %*% BcP
    loads.rn <- gaussian.loadings(round(mats/dt), mu.Z, Phi.Z - diag(N+1), Sigma.Z, rho0.cP*dt, c(0, rho1.cP)*dt, dt)
    Yrn <- rep(1, T) %*% loads.rn$A + Z %*% loads.rn$B
    Ytp <- Yhat - Yrn
    frn <- 2*Yrn[,mats==10] - Yrn[,mats==5]
    ftp <- 2*Ytp[,mats==10] - Ytp[,mats==5]
})

## ESE model
load("results/ese.RData")
data$istar.ese <- 100*apply(tau, 2, mean)
data$istar.lb <- 100*apply(tau, 2, quantile, .025)
data$istar.ub <- 100*apply(tau, 2, quantile, .975)
M <- nrow(AcP)
Yhat <- array(NA, c(T, J, M))
Yrn <- array(NA, c(T, J, M))
Ystar <- array(NA, c(T, J, M))
for (j in 1:M) {
    Yhat[,,j] <- rep(1, T) %*% AcP[j,,drop=FALSE] + cP[j,,] %*% BcP[j,,]
    cPstar.j <- rep(1, T) %*% Pbar[j,,drop=FALSE] + tau[j,] %*% gamma[j,,drop=FALSE]
    Ystar[,,j] <- rep(1, T) %*% AcP[j,,drop=FALSE] + cPstar.j %*% BcP[j,,]
    loads <- jsz.loadings(WN, diag(lamQ[j,]-1), kinfQ[j], makeCovMat(Sigma[j,]), mats, dt)
    Ptilde.j <- cP[j,,] - rep(1, T) %*% Pbar[j,,drop=FALSE] - tau[j,] %*% gamma[j,,drop=FALSE]
    Yrn[,,j] <- getYrn(tau[j,], Ptilde.j, matrix(Phi[j,], N, N), sigtau2[j], makeCovMat(Sigma.tilde[j,]), loads$rho0.cP, loads$rho1.cP, mats, dt)
}
data$yhat.ese <- 100*rowMeans(Yhat[,mats==10,])
data$fhat.ese <- 2*data$yhat.ese - data$yhat.ese
ystar <- 100*Ystar[,mats==10,]
data$ystar.ese <- rowMeans(ystar)
data$ystar.lb <- apply(ystar, 1, quantile, 0.025)
data$ystar.ub <- apply(ystar, 1, quantile, 0.975)
data$ytilde.ese <- 100*Y[,mats==10] - data$ystar.ese
frn <- 100*(2*Yrn[,mats==10,] - Yrn[,mats==5,])
data$frn.ese <- rowMeans(frn)
data$frn.lb <- apply(frn, 1, quantile, 0.025)
data$frn.ub <- apply(frn, 1, quantile, 0.975)
data$ftp.ese <- data$fhat.ese - data$frn.ese

## adding OSE results to data.frame
data$f <- 100*(2*Y[,mats==10] - Y[,mats==5])
data$ystar.ose <- 100*mod$Ystar[,mats==10]
data$yhat.ose <- 100*mod$Yhat[,mats==10]
data$ytilde.ose <- 100*mod$Ytilde[,mats==10]
data$frn.ose <- 100*mod$frn
data$ftp.ose <- 100*mod$ftp
data$frn.jsz <- 100*mod.jsz$frn
data$ftp.jsz <- 100*mod.jsz$ftp

##################################################
## Figure 4 - i*
pdf("figures/dtsm_istar.pdf", width=7, height=5.5, pointsize=10)
par(mar=c(3,4,1,1), mgp=c(2,.6,0))
cols <- c("black", "steelblue", "red", "green", "blue")
yrange <- range(data$y10, data$istar.lb, data$istar.ub, 0)
plot(data$date, data$y10, ylim=yrange, col=cols[1], type="l", lwd=2, xlab="", ylab="Percent", xaxs = "i", yaxs = "i")
polygon(x = c(data$date, rev(data$date)), y = c(data$istar.lb, rev(data$istar.ub)), density=NA, col="lightblue")
lines(data$date, data$y10, col=cols[1], lwd=2)
lines(data$date, data$istar.ese, col=cols[2], lwd=2)
lines(data$date, data$istar, col=cols[3], lwd=2)
legend("topright", c("Ten-year yield", "ESE model estimate of i*", "Real-time proxy estimate of i*"), col=cols, lwd=2)
box()
dev.off()

##################################################
## Figure 5 - trend and cycle of ten-year yield
pdf("figures/dtsm_y10.pdf", width=7, height=5.5, pointsize=10)
par(mar=c(3,4,1,1), mgp=c(2,.6,0))
cols <- c("gray", "black", "steelblue", "red", "steelblue", "red")
lwds <- c(2,1,2,2,1,1)
ltys <- c(1,2,1,1,1,1)
yrange <- range(0, data$y10, data$ystar.ub, data$ytilde.ese, data$ytilde.ose)
plot(data$date, data$y10, type="l", ylim=yrange, main="", col=cols[1], lwd=lwds[1], lty=ltys[1], ylab="Percent", xlab="", xaxs="i", yaxs="i")
polygon(x = c(data$date, rev(data$date)), y = c(data$ystar.lb, rev(data$ystar.ub)), density=NA, col="lightblue")
lines(data$date, data$y10, col=cols[1], lwd=lwds[1],lty=ltys[1])
lines(data$date, data$yhat.ose, col=cols[2], lwd=lwds[2], lty=ltys[2])
lines(data$date, data$ystar.ese, col=cols[3], lwd=lwds[3], lty=ltys[3])
lines(data$date, data$ystar.ose, col=cols[4], lwd=lwds[4], lty=ltys[4])
lines(data$date, data$ytilde.ese, col=cols[5], lwd=lwds[5], lty=ltys[5])
lines(data$date, data$ytilde.ose, col=cols[6], lwd=lwds[6], lty=ltys[6])
abline(h=0, lty=2)
legend("topright", c("Observed yield", "Fitted yield", "Trend component (ESE model)", "Trend component (OSE model)", "Cycle component (ESE model)", "Cycle component (OSE model)"), lwd=lwds, lty=ltys, col=cols, bg="white")
box()
dev.off()

##################################################
## Figure 7 - term premium
pdf("figures/dtsm_tp.pdf", width=7, height=5, pointsize=10)
par(mfrow=c(1,2), mar = c(2,3,2,0.5), mgp=c(2,.6,0))
yrange <- range(0, data$f, data$frn.jsz, data$ftp.jsz, data$frn.ose, data$ftp.ose, data$frn.ese, data$ftp.ese)
## left panel
cols <- c("black", "steelblue", "red")
plot(data$date, data$f, type="l", ylim=yrange, col=cols[1], lwd=2, ylab="Percent", xlab="", xaxs="i", yaxs="i")
title("Fixed-endpoint model", font.main=1, cex.main=1.1)
lines(data$date, data$frn.jsz, col=cols[2], lwd=2)
lines(data$date, data$ftp.jsz, col=cols[3], lwd=2)
legend("topright", c("Forward rate", "Risk-neutral rate", "Term premium"), col=cols, lwd=2)
## right panel
cols <- c("black", "steelblue", "forestgreen", "red", "orange")
lwds <- c(2,2,2,2,2)
par(mar = c(2,2,2,1.5))
plot(data$date, data$f, ylim=yrange, col=cols[1], type="l", lwd=lwds[1], xlab="", ylab="", xaxs = "i", yaxs = "i")
title("Shifting-endpoint models", font.main=1, cex.main=1.1)
lines(data$date, data$frn.ese, col=cols[2], lwd=lwds[2])
lines(data$date, data$frn.ose, col=cols[3], lwd=lwds[3])
lines(data$date, data$ftp.ese, col=cols[4], lwd=lwds[4])
lines(data$date, data$ftp.ose, col=cols[5], lwd=lwds[5])
box()
legend("topright", c("Forward rate", "Risk-neutral rate (ESE)", "Risk-neutral rate (OSE)", "Term premium (ESE)", "Term premium (OSE)"), col=cols, lty=1, lwd=2, bg="white")
dev.off()

###########################################################
## summary statistics for secular decline in interest rates
years <- c(1980:1982, 2015:2017)
yearagg <- data[data$year %in% years,] %>%
    select(year, f, frn.jsz, ftp.jsz, frn.ose, ftp.ose, frn.ese, ftp.ese) %>%
    group_by(year) %>%
    summarise_all(funs(mean))
print(round(yearagg,3))
tmp <- yearagg %>%
    group_by(year > 2000) %>%
    summarise_all(funs(mean)) %>%
    select(year, f, ftp.jsz, ftp.ose, ftp.ese)
tmp <- as.matrix(rbind(tmp, tmp[2,]-tmp[1,]))
rownames(tmp) <- tmp[,1]
tmp <- tmp[,-1]
print(round(tmp,1))
cat("Percent of decline accounted for by term premium:\n")
tmp2 <- tmp[3,2:4]/tmp[3,1]*100
print(round(tmp2, 0))
cat("Percent of decline accounted for by expectations:\n")
print(round(100-tmp2))
