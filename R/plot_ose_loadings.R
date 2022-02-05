## Plot loadings for OSE model on i*
## - Figure 6
## - Online Appendix Figure 3

graphics.off()
source("R/dtsm_fns.R")
load("results/ose.RData")

## Online Appendix Figure 3
## loadings of forward rates and yields on i* for very long maturities
H <- 400
floads <- numeric(H)  # loadings of f(n-1)
PhiQn <- diag(N)
for (n in 1:H) {
    floads[n] <- t(mod$rho1.cP) %*% PhiQn %*% mod$gamma
    PhiQn <- PhiQn %*% mod$PhiQ
}
yloads <- cumsum(floads) / 1:H
pdf("figures/istar_loadings_humpshape.pdf", width=6, height=4, pointsize=10)
par(mar=c(4,4,1,1), mgp=c(2,.6,0))
plot(floads, type="l", lwd=2, xlab="Maturity in quarters", ylab="Coefficient")
lines(yloads, lwd=2, col="steelblue")
legend("topright", c("Forward rate loadings", "Yield loadings"), col=c("black", "steelblue"), lty=1, lwd=2)
abline(h=1, lty=2)
abline(v=60, lty=2)
dev.off()

## Maturity (in quarters) when forward rate loadings drop below one
print(min(which(floads<1)))
## Peak of forward rate loadings
print(i <- which(floads==max(floads)))
print(floads[i])
## Peak of yield loadings
print(i <- which(yloads==max(yloads)))
print(yloads[i])

##################################################
## Figure 6: loadings of yields on i*
loads_istar <- drop(crossprod(mod$BcP, mod$gamma))
## regression in the data using i*-proxy
p <- 4
library(dynlm)
data$distar <- c(NA, diff(data$istar))
coef_istar <- numeric(J)
se_istar <- numeric(J)
for (j in 1:J) {
    data$yield <- Y[,j]
    ## reg <- lm(yield ~ istar, data)
    reg <- dynlm(yield ~ istar + L(distar, -p:p), data=ts(data))
    V <- sandwich::NeweyWest(reg, lag=6, prewhite=FALSE)
    coef_istar[j] <- reg$coef[2]
    se_istar[j] <- sqrt(diag(V))[2]
}
cat("Regression coefficient of y10 on i*:", coef_istar[mats==10], "SE =", se_istar[mats==10], "\n")
## simulate many small samples
cat("simulating 5000 small samples...\n")
B <- 5000
coef_shortsim <- matrix(NA, B, J)
for (i in 1:B) {
    if (i %% 500 == 0) cat(i, "\n")
    sim <- sim_unsp_model(T=nobs, mod)
    ## DOLS
    simdata <- data.frame(istar = sim$istar)
    simdata$distar <- c(NA, diff(simdata$istar))
    for (j in 1:J) {
        simdata$yield <- sim$Y[,j]
        reg <- dynlm(yield ~ istar + L(distar, -p:p), data=ts(simdata))
        coef_shortsim[i,j] <- reg$coef[2]
    }
}
coef_sim_q <- apply(coef_shortsim, 2, quantile, c(0.5, 0.025, 0.975))
coef_sim_mean <- colMeans(coef_shortsim)

pdf("figures/dtsm_loadings.pdf", width=7, height=5.5, pointsize=10)
par(mar=c(3,4,1,1), mgp=c(2,.6,0))
cols <- c("black", "steelblue", "forestgreen")
yrange <- range(0, max(coef_sim_q), max(coef_istar+1.96*se_istar))
plot(mats, coef_istar, type="l", lwd=2, ylim=yrange, ylab="Coefficient", xlab="Maturity (years)", xaxs = "i", yaxs = "i", col=cols[1])
polygon(c(mats, rev(mats)), c(coef_istar+1.96*se_istar, rev(coef_istar-1.96*se_istar)), density=NA, col="lightgray")
lines(mats, coef_istar, lwd=2, col=cols[1])
lines(mats, loads_istar, lwd=2, col=cols[2])
lines(mats, coef_sim_q[2,], lwd=2, lty=2, col=cols[2])
lines(mats, coef_sim_q[3,], lwd=2, lty=2, col=cols[2])
abline(h=1, lty=2)
legend("bottomright", c("Regression coefficients", "Model loadings: population", "Model loadings: small samples, 95% CI"), col=cols[c(1,2,2)], lwd=2, lty=c(1,1,2))
box()
dev.off()
