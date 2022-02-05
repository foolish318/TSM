## Simulate small samples of yields and trends and predict excess returns
## Table 4

source("R/util_fns.R")
source("R/dtsm_fns.R")
set.seed(616)

B <- 5000

## ESE model
load("results/ese.RData")
N <- nrow(WN)
T <- nrow(data)
J <- ncol(Y)

## OSE model
load("results/ose.RData")
OSE <- mod

## main results
R2.OSE <- matrix(NA, B, 3)
R2.ESE <- matrix(NA, B, 3)
R2.FE <- matrix(NA, B, 3)
M <- nrow(lamQ)
for (i in 1:B) {
    ## OSE model
    sim.OSE <- sim_unsp_model(nobs, OSE)
    res <- predictReturns(sim.OSE$Y, sim.OSE$istar, h=1)
    R2.OSE[i,] <- c(res, diff(res))
    ## ESE model
    ## draw parameters from posterior distribution
    j <- sample(M, 1)
    ESE <- list(AcP = AcP[j,,drop=FALSE], BcP = BcP[j,,],
                gamma = gamma[j,], Pbar = Pbar[j,],
                Phi = matrix(Phi[j,], N, N),
                sigma.tau = sqrt(sigtau2[j]),
                sigma.e = sqrt(sige2[j]),
                Sigma.tilde = makeCovMat(Sigma.tilde[j,]))
    sim.ESE <- sim_unsp_model(nobs, ESE)
    res <- predictReturns(sim.ESE$Y, sim.ESE$istar, h=1)
    R2.ESE[i,] <- c(res, diff(res))
    ## under the null hypothesis -- fixed endpoint
    sim.FE <- sim_jsz_model(nobs, mod.jsz)
    res <- predictReturns(sim.FE$Y, sim.OSE$istar, h=1) # using (unrelated) i* from OSE model
    R2.FE[i,] <- c(res, diff(res))
}
tbl <- data.frame(matrix(NA, 7, 4))
tbl[,1] <- c("Data", "\\emph{FE} model", "", "\\emph{OSE} model", "", "\\emph{ESE} model", "")
## actual data
res <- predictReturns(Y, data$istar, h=1)
tbl[1,2:4] <- sprintf("%4.2f", c(res, res[2]-res[1]))
## FE model simulated data
tbl[2,2:4] <- sprintf("%4.2f", colMeans(R2.FE))
for (i in 1:3)
    tbl[3, 1+i] <- sprintf("[%4.2f, %4.2f]", quantile(R2.FE[,i], 0.025), quantile(R2.FE[,i], 0.975))
## OSE model simulated data
tbl[4,2:4] <- sprintf("%4.2f", colMeans(R2.OSE))
for (i in 1:3)
    tbl[5, 1+i] <- sprintf("[%4.2f, %4.2f]", quantile(R2.OSE[,i], 0.025), quantile(R2.OSE[,i], 0.975))
## ESE model simulated data
tbl[6,2:4] <- sprintf("%4.2f", colMeans(R2.ESE))
for (i in 1:3)
    tbl[7, 1+i] <- sprintf("[%4.2f, %4.2f]", quantile(R2.ESE[,i], 0.025), quantile(R2.ESE[,i], 0.975))
cat("Quarterly returns:\n")
print(tbl)

sink("tables/dtsm_returns.tex")
print(xtable::xtable(tbl, digi=3),
      include.rownames=FALSE, include.colnames=FALSE, only.contents=TRUE,
      sanitize.text.function=function(x){x}, hline.after=NULL)
sink()

