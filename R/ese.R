## Estimation of DTSM with unspanned endpoint -- ESE model
## MCMC sampler

library(KFAS)
library(MCMCpack)
library(mvtnorm)
source("R/dtsm_fns.R")
source("R/data_fns.R")
source("R/util_fns.R")

## normalization for ESE is different than for OSE
## (these are also defined in dtsm_fns.R but overwritten here)
makegamma <- function(a, dummy) c(1, a)
makePbar <- function(p, dummy, dummy2) c(0, p)

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==1) {
    SEED <- as.numeric(args[1])
}
if (!exists("SEED")) SEED <- 616

filename <- paste0("results/ese_", SEED, ".RData")

SIGMA <- 0.06
ALPHA <- 100
LAMBDA1 <- 0.2
LAMBDA2 <- 1.0

N <- 3 # number of factors
ltri <- lower.tri(matrix(0, N, N), diag=TRUE)

cat("# Estimating model\n")
cat("#", filename, "\n")
## load data
data <- loadData()
data$istar <- data$pistar.ptr + data$rstar.realtime
mats <- c(0.25, 0.5, 1:15)
yield.cols <- paste0("y", mats)
dt <- 1/4   # quarterly frequency

## yields
Y <- as.matrix(data[yield.cols])/100
J <- ncol(Y)
W <- 1*cbind(mats==0.25, mats==2, mats==10)
WN <- t(W)   # N x J
nobs <- nrow(data)

## MCMC chains
set.seed(SEED)
cat("Random seed:", SEED, "\n")
M <- 100000 # length of MCMC chain
cat("length of chain:", M, "\n")
scale.thetaQ <- c(100, 0.1, 0.1, 0.1)
scale.ap <- c(1, 1, 10, 10)
scale.sigtau2 <- 0.1
scale.sigma <- c(100, 100, 100, 100, 100, 100)
scale.theta <- c(scale.thetaQ, scale.ap, scale.sigtau2, scale.sigma)
nu <- 5            # for independence proposals, multivariate t draws

## ZLB
## 3m yield: Dec-2008 to Nov-2015 (based on the fed funds rate
## 6m yield: same
## 1y yield: Oct-2011 to Nov-2014 (beginning based on Swanson-Williams, end based on level of yield)
## 2y yield: Oct-2012 to May-2013 (same)
zlb.ind <- matrix(0, nobs, J)
zlb.ind[data$yyyymm >= 200812 & data$yyyymm < 201512, mats==0.25] <- 1 # 3m yield
zlb.ind[data$yyyymm >= 200812 & data$yyyymm < 201512, mats==0.5] <- 1  # 6m yield
zlb.ind[data$yyyymm >= 201110 & data$yyyymm < 201411, mats==1] <- 1    # 1y yield
zlb.ind[data$yyyymm >= 201210 & data$yyyymm < 201305, mats==2] <- 1    # 1y yield

## get starting values
est.sample <- data$yyyymm<200800 # estimation sample - leave out ZLB period
cat("# 1) standard JSZ model -- stationary, observed cP -- to get Q-dynamics\n")
mod.jsz <- estimate_jsz(Y, est.sample, WN, mats, dt)
cat("# 2) starting values for gamma, Pbar, Sigma.tilde\n")
cP <- Y %*% t(WN)
pars.start <- startval_ose(Y[est.sample,], data$istar[est.sample]/100, mod.jsz, WN, mats, dt)
cat("# after JSZ etc:\n")
cat("sigtau2 =", pars.start$sigtau2, "\n")
cat("P-eigenvalues:\n")
print(abs(eigen(pars.start$Phi)$values))
cat("Q-eigenvalues:\n")
print(pars.start$lamQ)
cat("gamma:\n")
print(pars.start$gamma)
cat("Pbar:\n")
print(pars.start$Pbar)

res <- lsfit(data$istar, cP)
psi.data <- diag(cov(res$residuals))
Psi.data <- diag(psi.data)
sigma <- sqrt(psi.data)
print(sigma*100)

prior <- within(list(), {
    ## sigma_tau^2 ~ IG(alpha_eta/2, delta_eta/2)
##################################################
    alpha_eta  <- ALPHA     # tightness. Del Negro et al use 100
    mean_sigtau2  <- SIGMA^2/400     # 0.1^2*400 - SD of change over 100 years is 0.1=10%
    ## interpretation/compare: Del Negro pi*: (2%)^2=Var(change over 400 quarters)=400*Var(change over one quarter)
##################################################
    delta_eta <- mean_sigtau2 * (alpha_eta - 2)
    ## VAR: independent Normal-Wishart, Minnessota
    ## Sigma_tilde ~ IW(kappa, Psi)
    kappa <- N + 2                       # low precision
    Psi <- Psi.data # pars.start$Sigma              # prior mean = Psi/(kappa-N-1) = Psi
    ## vec(Phi) ~ N(0; Omega)
    Phi <- matrix(0, N, N) # mean
    ## SDs -- Minnessota
    lambda_1 <- LAMBDA1 ## Del Negro et al used 0.2
    lambda_2 <- LAMBDA2
    ## SD(Phi_ij) = sigma_i*lambda_1*lambda_2/sigma_j if i \neq j
    ##            = lambda_1 if i == j
    SDs <- diag(3)*lambda_1
    for (i in 1:3)
        for (j in setdiff(1:3, i))
            SDs[i,j] <- lambda_1*lambda_2*sigma[i]/sigma[j]
    omega <- as.numeric(SDs^2)
    Omega <- diag(omega)
    Omega_inv <- diag(1 / omega)
    ## ap ~ N
    ap_mean <- c(rowSums(WN)[2:3],0,0)
    ap_cov <- diag(c(0.2, 0.2, 0.05^2, 0.05^2))  # before: 0.01^2, 0.02^2
    ap_cov_inv <- solve(ap_cov)
    ## sige2 = meas. error variance ~ IG(alpha_e/2, delta_e/2)
    alpha_e <- 4
    mean_sige2 <- 0.001^2  # ten basis points
    delta_e <- mean_sige2 * (alpha_e - 2)
})

## print(prior[c("alpha_eta", "mean_sigtau2", "lambda_1", "lambda_2", "ap_mean", "ap_cov")])
cat("# Priors:\n")
cat(100*sqrt(prior$mean_sigtau2*400), "percent per 100 years (prior mean sigtau2)\n")
cat("alpha =", prior$alpha_eta, "\n")
cat("lambda_1 =", prior$lambda_1, "\n")
cat("lambda_2 =", prior$lambda_2, "\n")

ssm <- function(Y, pars) {
    ## state vector: 5 state variables
    ##  (1) constant, (2) i*, (3-5) Ptilde(t)
    N <- length(pars$lamQ)
    J <- ncol(Y)
    ## measurement equation
    Z <- cbind(t(pars$AcP + crossprod(pars$Pbar, pars$BcP)), # constant
               t(crossprod(pars$gamma, pars$BcP)),         # tau
               t(pars$BcP))                                # Ptilde(t)
    H <- diag(J)*pars$sige2
    ## transition equation
    T <- diag(2+N)
    T[3:5,3:5] <- pars$Phi
    R <- matrix(0, 2+N, N+1); R[2:(N+2),] <- diag(N+1)
    Q <- matrix(0, 1+N, 1+N); Q[1,1] <- pars$sigtau2; Q[-1,-1] <- pars$Sigma.tilde
    ## initial conditions
    a1 <- c(1, 0.02, rep(0, 3))
    P1 <- diag(c(0, rep(1, 4)))
    ## model
    Y[zlb.ind==1] <- NA  #  yields that are constrained at ZLB treated as unobserved
    SSModel(Y ~ -1 + SSMcustom(Z=Z, T=T, R=R, Q=Q, a1=a1, P1=P1), H=H)
}

pars2theta <- function(pars) {
    rval <- scale.theta * c(pars$kinfQ, log(-pars$dlamQ/(1+pars$dlamQ)), pars$a, pars$p, log(pars$sigtau2), t(chol(pars$Sigma))[ltri])
    c(rval, as.numeric(pars$Phi), log(pars$sige2))
}

theta2pars <- function(theta) {
    theta[1:15] <- theta[1:15] / scale.theta
    stopifnot(length(theta)==15+9+1)
    pars <- list(kinfQ = theta[1],
                 dlamQ = -exp(theta[2:4])/(1+exp(theta[2:4])),
                 a = theta[5:6],
                 p = theta[7:8],
                 sigtau2 = exp(theta[9]))
    Sigmachol <- matrix(0, N, N)
    Sigmachol[ltri] <- theta[10:15]; theta <- tail(theta, -15)
    pars$Sigma <- Sigmachol %*% t(Sigmachol)
    pars$lamQ <- 1 + cumsum(pars$dlamQ)
    pars$Phi <- matrix(head(theta, N^2), N, N); theta <- tail(theta, -N^2)
    pars$sige2 <- exp(theta)
    pars
}

## log density of prior for Phi
logprior_Phi <- function(Phi) {
    phi <- as.numeric(Phi)
    n <- length(phi)
    -n/2*log(2*pi)-1/2*sum(log(diag(prior$Omega)))-1/2*sum(phi^2/diag(prior$Omega))
}

obj <- function(theta, prior, Y, WN, mats, dt) {
    pars <- theta2pars(theta)
    if (pars$sigtau2 > 0.1) return(1e6)
    if (pars$sige2 > 0.1) return(1e6)
    if (any(abs(eigen(pars$Phi)$values)>1)) return(1e6)
    if (any(is.na(pars$lamQ)) || any(pars$lamQ < 0)) return(1e6)
    loads <- jsz.loadings(WN, diag(pars$lamQ-1), pars$kinfQ, pars$Sigma, mats, dt)
    pars$gamma <- c(1, pars$a)
    pars$Pbar <- c(0, pars$p)
    pars$Sigma.tilde <- pars$Sigma - outer(pars$gamma, pars$gamma)*pars$sigtau2
    if (invalidCovMat(pars$Sigma.tilde)) return(1e6)
    pars <- c(pars, loads[c("AcP", "BcP")])
    mod <- ssm(Y, pars)
    llk <- KFS(mod, smoothing="none")$logLik
    lprior <- log(dinvgamma(pars$sigtau2, prior$alpha_eta/2, prior$delta_eta/2)) +
        logprior_Phi(pars$Phi) +
        log(diwish(pars$Sigma.tilde, prior$kappa, prior$Psi)) +
        log(dinvgamma(pars$sige2, prior$alpha_e/2, prior$delta_e/2)) +
        dmvnorm(c(pars$a, pars$p), prior$ap_mean, prior$ap_cov, log=TRUE)
    - (llk + lprior)
}

stopifnot(all.equal(logprior_Phi(pars.start$Phi),
                    dmvnorm(as.numeric(pars.start$Phi), mean = as.numeric(prior$Phi), sigma = prior$Omega, log=TRUE)))
stopifnot(all.equal(logprior_Phi(pars.start$Phi), sum(dnorm(as.numeric(pars.start$Phi), sd=as.numeric(prior$SDs), log=TRUE))))

cat("# 3) find mode of posterior\n")
theta.start <- pars2theta(pars.start)
theta.check <- pars2theta(theta2pars(theta.start))
stopifnot(all.equal(theta.start, theta.check))
lp.start <- -obj(theta.start, prior, Y, WN, mats, dt)
cat("log-posterior at starting values:", lp.start, "\n")

theta.mode <- get_optim(theta.start, obj, prior, Y, WN, mats, dt, trace=1)
cat("look at scaling:\n")
print(cbind(theta.mode[1:length(scale.theta)]/scale.theta, scale.theta, theta.mode[1:length(scale.theta)]))
pars <- theta2pars(theta.mode)
optimal.lp <- -obj(theta.mode, prior, Y, WN, mats, dt)
cat("# log-posterior at optimal values:", optimal.lp, "\n")
loads <- jsz.loadings(WN, diag(pars$lamQ-1), pars$kinfQ, pars$Sigma, mats, dt)
pars.mode <- c(pars, loads[c("AcP", "BcP", "BX")])
pars.mode$gamma <- c(1, pars.mode$a)
pars.mode$Pbar <- c(0, pars.mode$p)
pars.mode$Sigma.tilde <- pars.mode$Sigma - outer(pars.mode$gamma, pars.mode$gamma)*pars.mode$sigtau2

## random starting values for MCMC chain
H <- numDeriv::hessian(obj, theta.mode, method="Richardson", method.args=list(), prior, Y, WN, mats, dt)
theta.cov <- makePD(solve(H))
cat("Objective function at mode:", obj(theta.mode, prior, Y, WN, mats, dt), "\n")
valid <- FALSE
while (!valid) {
    theta <- drop(theta.mode + 2*t(chol(theta.cov)) %*% rnorm(length(theta.mode)))
    valid <- (obj(theta, prior, Y, WN, mats, dt) < 1e+6)
}
cat("Objective function at random start:", obj(theta, prior, Y, WN, mats, dt), "\n")
pars <- theta2pars(theta)
loads <- jsz.loadings(WN, diag(pars$lamQ-1), pars$kinfQ, pars$Sigma, mats, dt)
pars <- c(pars, loads[c("AcP", "BcP", "BX", "rho0.cP", "rho1.cP")])
pars$gamma <- makegamma(pars$a, loads$rho1.cP)
pars$Pbar <- makePbar(pars$p, loads$rho0.cP, loads$rho1.cP)
pars$Sigma.tilde <- pars$Sigma - outer(pars$gamma, pars$gamma)*pars$sigtau2

## double-check that first factor is short rate
loads <- jsz.loadings(WN, diag(pars$lamQ-1), pars$kinfQ, pars$Sigma, mats, dt)
stopifnot(all.equal(Y[,1], drop(loads$rho0.cP) + drop(Y %*% W %*% loads$rho1.cP), check.attributes=FALSE))

cat("# Posterior mode:\n")
cat("sigtau2 =", pars.mode$sigtau2, "\n")
cat("vs prior mean", prior$mean_sigtau2, "\n")
cat("P-eigenvalues:\n")
print(abs(eigen(pars.mode$Phi)$values))
cat("Q-eigenvalues:\n")
print(pars.mode$lamQ)
cat("gamma:\n")
print(pars.mode$gamma)
cat("Pbar:\n")
print(pars.mode$Pbar)

llk_Q <- function(Y, cP, AcP, BcP, sige2) {
    Yhat <- rep(1, nobs) %*% rbind(AcP) + cP %*% BcP
    Q.errors <- Y - Yhat
    Q.errors[zlb.ind==1] <- 0
    -.5*sum(Q.errors^2)/sige2
}

llk_P <- function(tau, cP, Pbar, gamma, Phi, Sigma.tilde, sigtau2) {
    ## calculate value of P-likelihood
    Ptilde <- cP - rep(1, nobs) %*% rbind(Pbar) - tau %*% rbind(gamma)
    e <- t(Ptilde[2:nobs,]) - Phi %*% t(Ptilde[1:(nobs-1),])  # N x (T-1)  # VAR errors (Ptilde only)
    llk_var <- -.5*(nobs-1)*log(det(Sigma.tilde))-.5*sum(e*(solve(Sigma.tilde, e)))
    llk_rw <- -.5*(nobs-1)*log(sigtau2) -.5*sum(diff(tau)^2/sigtau2)
    llk_var + llk_rw
}

prop_thetaQ <- function(thetaQ, Y, pars) {
    ## Globals: WN, mats, dt, nobs
    thetaQ <- thetaQ / scale.thetaQ
    prop <- list(kinfQ = thetaQ[1],
                 dlamQ = -exp(thetaQ[2:4])/(1+exp(thetaQ[2:4])))
    prop$lamQ <- 1 + cumsum(prop$dlamQ)
    if (any(is.na(prop$lamQ)) || any(prop$lamQ < 0) || all(prop$lamQ < 0.9)) return(NULL)
    ## Q-likelihood -- loadings and fitting errors change
    loads <- jsz.loadings(WN, diag(prop$lamQ-1), prop$kinfQ, pars$Sigma, mats, dt)
    prop <- c(prop, loads[c("AcP", "BcP", "BX")])
    prop$llkQ <- llk_Q(Y, pars$cP, prop$AcP, prop$BcP, pars$sige2)
    prop
}

obj_thetaQ <- function(thetaQ, Y, pars) {
    ## value of neg. log cond. posterior
    res <- prop_thetaQ(thetaQ, Y, pars)
    if (is.null(res)) return(1e6)
    rval <- -res$llkQ
}

prop_sigma <- function(sigma, Y, pars) {
    ## Globals: WN, mats, dt, nobs, prior
    Sigmachol <- matrix(0, N, N)
    Sigmachol[ltri] <- sigma/scale.sigma
    prop <- list(Sigma = Sigmachol %*% t(Sigmachol))
    if (invalidCovMat(prop$Sigma)) return(NULL)
    prop$Sigma.tilde <- prop$Sigma - outer(pars$gamma, pars$gamma)*pars$sigtau2
    if (invalidCovMat(prop$Sigma.tilde)) return(NULL)
    loads <- jsz.loadings.post(WN, pars$BX, diag(pars$lamQ-1), pars$kinfQ, prop$Sigma, mats, dt)
    prop <- c(prop, loads[c("AcP")])  # BcP and rho1.cP not affected by change in Sigma
    prop$llkQ <- llk_Q(Y, pars$cP, prop$AcP, pars$BcP, pars$sige2)
    prop$llkP <- llk_P(pars$tau, pars$cP, pars$Pbar, pars$gamma, pars$Phi, prop$Sigma.tilde, pars$sigtau2)
    prop$lprior <- log(diwish(prop$Sigma.tilde, prior$kappa, prior$Psi))
    if (!is.finite(prop$lprior)) return(NULL)
    prop
}

obj_sigma <- function(sigma, Y, pars) {
    res <- prop_sigma(sigma, Y, pars)
    if (is.null(res)) return(1e6)
    -(res$llkQ + res$llkP + res$lprior)
}

prop_sigtau2 <- function(sigtau2, pars) {
    ## Globals: WN, mats, dt, nobs, prior
    prop <- list(sigtau2 = sigtau2,
                 Sigma.tilde = pars$Sigma - outer(pars$gamma, pars$gamma)*sigtau2)
    if (invalidCovMat(prop$Sigma.tilde)) return(NULL)
    prop$llkP <- llk_P(pars$tau, pars$cP, pars$Pbar, pars$gamma, pars$Phi, prop$Sigma.tilde, prop$sigtau2)
    prop$lprior <- log(dinvgamma(prop$sigtau2, prior$alpha_eta/2, prior$delta_eta/2)) +
        log(diwish(prop$Sigma.tilde, prior$kappa, prior$Psi))
    if (!is.finite(prop$lprior)) return(NULL)
    prop
}

obj_sigtau2 <- function(lsigtau2, pars) {
    res <- prop_sigtau2(exp(lsigtau2), pars)
    if (is.null(res)) return(1e6)
    -(res$llkP + res$lprior)
}

prop_ap <- function(ap, pars) {
    ## Globals: WN, mats, dt, nobs
    prop <- list(a = ap[1:2], p = ap[3:4])
    prop$gamma <- c(1, prop$a)
    prop$Pbar <- c(0, prop$p)
    prop$Sigma.tilde <- pars$Sigma - outer(prop$gamma, prop$gamma)*pars$sigtau2
    if (invalidCovMat(prop$Sigma.tilde)) return(NULL)
    prop$llkP <- llk_P(pars$tau, pars$cP, prop$Pbar, prop$gamma, pars$Phi, prop$Sigma.tilde, pars$sigtau2)
    prop$lprior <- log(diwish(prop$Sigma.tilde, prior$kappa, prior$Psi)) +
        dmvnorm(ap, prior$ap_mean, prior$ap_cov, log=TRUE)
    if (!is.finite(prop$lprior)) return(NULL)
    prop
}

obj_ap <- function(ap, pars) {
    res <- prop_ap(ap, pars)
    if (is.null(res)) return(1e6)
    -(res$llkP + res$lprior)
}

cat("# 3) MCMC sampler...\n")
kinfQ <- numeric(M)*NA
lamQ <- matrix(NA, M, N)
gamma <- matrix(NA, M, N)
Pbar <- matrix(NA, M, N)
## apvar <- matrix(NA, M, 4)
## apvar_disp <- matrix(NA, M, 4)
sigtau2 <- numeric(M)*NA
sige2 <- numeric(M)*NA
Phi <- matrix(NA, M, N^2)
eigenval <- numeric(M)*NA
Sigma <- matrix(NA, M, N*(N+1)/2)
Sigma.tilde <- matrix(NA, M, N*(N+1)/2)
AcP <- matrix(NA, M, J)
BcP <- array(NA, c(M, N, J))
tau <- matrix(NA, M, nrow(data))
cP <- array(NA, c(M, nrow(data), N))
alpha <- matrix(0, M, 4) # Metropolis-Hastings acceptance probabilities
start <- Sys.time()
for (j in 1:M) {
    if (j %% 1000 == 0) {
        cat("Iteration", j, "\n")
        cat("Duration so far:", format(Sys.time()-start), "\n")
        cat("Mean of acceptance prob:", round(colMeans(head(alpha, j-1)), 3), "\n")
        cat("Cumul. mean eigenvalue:", mean(head(eigenval, j-1)), "\n")
    }

    ## ###########################################
    ## 1) draw states -- i-star and Ptilde -> cP
    sim <- simulateSSM(ssm(Y, pars), type="states")  ## (N+1) x T
    pars$tau <- sim[,2,1]
    pars$cP <- rep(1, nobs) %*% rbind(pars$Pbar) + pars$tau %*% rbind(pars$gamma) + sim[,3:5,1]
    pars$llkQ <- llk_Q(Y, pars$cP, pars$AcP, pars$BcP, pars$sige2)
    pars$llkP <- llk_P(pars$tau, pars$cP, pars$Pbar, pars$gamma, pars$Phi, pars$Sigma.tilde, pars$sigtau2)

    ## ###########################################
    ## 2) thetaQ - kinfQ, lamQ -- IMH
    thetaQ_current <- scale.thetaQ * c(pars$kinfQ, log(-pars$dlamQ/(1+pars$dlamQ)))
    ## res <- prop_thetaQ(thetaQ_current, Y, pars)
    ## stopifnot(isTRUE(all.equal(res$llkQ, pars$llkQ)))
    res <- optim(thetaQ_current, obj_thetaQ, gr=NULL, Y, pars, method="L-BFGS-B", hessian=TRUE)
    thetaQ_mean <- res$par
    thetaQ_mvtscale <- makePD(solve(res$hessian)) * (nu-2)/nu
    thetaQ_prop <- thetaQ_mean + rmvt(1, sigma = thetaQ_mvtscale, df = nu)
    prop <- prop_thetaQ(thetaQ_prop, Y, pars)
    if (!is.null(prop)) {
        ## ratio of log proposal density
        lr_prop <- dmvt(thetaQ_current - thetaQ_mean, sigma = thetaQ_mvtscale, df = nu, log = TRUE) -
            dmvt(thetaQ_prop - thetaQ_mean, sigma = thetaQ_mvtscale, df = nu, log = TRUE)
        (alpha[j, 1] <- min(exp(prop$llkQ - pars$llkQ + lr_prop), 1))
        if (runif(1) < alpha[j, 1]) {
            nms <- c("kinfQ", "dlamQ", "lamQ", "AcP", "BcP", "BX", "llkQ")
            ## stopifnot(setequal(names(prop), nms))
            pars[nms] <- prop[nms]
        }
    }

    ## ###########################################
    ## 3) Sigma -- Independence MH
    sigma_current <- scale.sigma * as.numeric(t(chol(pars$Sigma))[ltri])
    lprior_current <- log(diwish(pars$Sigma.tilde, prior$kappa, prior$Psi))
    ## res <- prop_sigma(sigma_current, Y, pars)
    ## stopifnot(isTRUE(all.equal(res$llkQ + res$llkP + res$lprior, pars$llkQ + pars$llkP + lprior_current)))
    res <- optim(sigma_current, obj_sigma, gr=NULL, Y, pars, method="L-BFGS-B", hessian=TRUE)
    sigma_mean <- res$par
    sigma_mvtscale <- makePD(solve(res$hessian)) * (nu-2)/nu
    sigma_prop <- sigma_mean + rmvt(1, sigma = sigma_mvtscale, df = nu)
    prop <- prop_sigma(sigma_prop, Y, pars)
    if (!is.null(prop)) {
        lr_prop <- dmvt(sigma_current - sigma_mean, sigma = sigma_mvtscale, df = nu, log = TRUE) -
            dmvt(sigma_prop - sigma_mean, sigma = sigma_mvtscale, df = nu, log = TRUE)
        lr_prior <-  prop$lprior - lprior_current; prop$lprior <- NULL
        (alpha[j, 2] <- min(exp(prop$llkQ - pars$llkQ + prop$llkP - pars$llkP + lr_prop + lr_prior), 1))
        if (runif(1) < alpha[j, 2]) {
            nms <- c("Sigma", "Sigma.tilde", "AcP", "llkQ", "llkP")
            ## stopifnot(setequal(names(prop), nms))
            pars[nms] <- prop[nms]
        }
    }

    ## ###########################################
    ## 4) ap (gamma/Pbar) -- Independence MH
    ap_current <- c(pars$a, pars$p)
    lprior_current <- log(diwish(pars$Sigma.tilde, prior$kappa, prior$Psi)) +
        dmvnorm(ap_current, prior$ap_mean, prior$ap_cov, log=TRUE)
    ## stopifnot(isTRUE(all.equal(-obj_ap(ap_current, pars), pars$llkP + lprior_current)))

    ## ## ## proposal based on conditional posterior that ignores Sigma_tilde implications
    ## Sigma_inv <- solve(pars$Sigma.tilde)
    ## y <- as.numeric(t(pars$cP[2:nobs,] - pars$cP[1:(nobs-1),] %*% t(pars$Phi)))
    ## D <- cbind(pars$tau[2:nobs], 1, -pars$tau[1:(nobs-1)], -1)
    ## H <- rbind(0, diag(4)[1:2,], 0, diag(4)[3:4,])
    ## h <- c(1,0,0,0,0,0)
    ## U <- rbind(diag(6), kronecker(diag(2), pars$Phi))
    ## z <- y - kronecker(D, diag(3)) %*% U %*% h
    ## tmp <- t(H) %*% t(U) %*% kronecker(crossprod(D), Sigma_inv) %*% U %*% H
    ## apvar[j,] <- diag(solve(prior$ap_cov_inv + tmp))
    ## apvar_disp[j,] <- diag(solve(tmp))
    ## (ap_mean <- drop(ap_cov %*% (prior$ap_cov_inv %*% prior$ap_mean + t(H) %*% t(U) %*% kronecker(t(D), Sigma_inv) %*% z)))
    ## (ap_scale <- t(chol(ap_cov)))
    ## gibbs_mean <- ap_mean

    ## proposal based on mode and Hessian of conditional posterior, found numerically
    res <- optim(ap_current, obj_ap, gr=NULL, pars, method="L-BFGS-B", hessian=TRUE)
    ap_mean <- res$par
    ap_cov <- makePD(solve(res$hessian))
    ## cbind(gibbs_mean, ap_mean)
    ## obj_ap(gibbs_mean, pars)
    ## obj_ap(ap_mean, pars)
    ## prop_ap(gibbs_mean, pars)
    ## prop_ap(ap_mean, pars)
    ## gamma <- c(1, gibbs_mean[1:2])
    ## Sigma.tilde <- pars$Sigma - outer(gamma, gamma)*pars$sigtau2
    ## eigen(Sigma.tilde)
    ## eigen(pars$Sigma.tilde)

    ap_prop <- drop(ap_mean + t(chol(ap_cov)) %*% rnorm(4))
    prop <- prop_ap(ap_prop, pars)
    if (!is.null(prop)) {
        lr_prop <- logrationorm(ap_current, ap_prop, ap_mean, ap_cov)
        lr_prior <- prop$lprior - lprior_current; prop$lprior <- NULL
        (alpha[j, 3] <- min(exp(prop$llkP - pars$llkP + lr_prop + lr_prior), 1))
    }
    if (runif(1) < alpha[j, 3]) {
        nms <- c("a", "gamma", "p", "Pbar", "Sigma.tilde", "llkP")
        ## stopifnot(setequal(names(prop), nms))
        pars[nms] <- prop[nms]
    }

    ## ###########################################
    ## 5) sigtau2 -- Independence MH
    lprior_current <- log(dinvgamma(pars$sigtau2, prior$alpha_eta/2, prior$delta_eta/2)) +
        log(diwish(pars$Sigma.tilde, prior$kappa, prior$Psi))
    ## res <- prop_sigtau2(pars$sigtau2, pars)
    ## stopifnot(isTRUE(all.equal(res$llkP + res$lprior, pars$llkP + lprior_current)))

    ## ## proposal based on conditional posterior that ignores Sigma_tilde implications
    ## alpha_post <- (prior$alpha_eta + length(pars$tau) - 1)
    ## delta_post <- (prior$delta_eta + sum(diff(pars$tau)^2))
    ## sigtau2_prop <- rinvgamma(1, alpha_post/2, delta_post/2)
    ## lr_prop <- log(dinvgamma(pars$sigtau2, alpha_post/2, delta_post/2)) - log(dinvgamma(sigtau2_prop, alpha_post/2, delta_post/2))

    ## proposal based on mode and Hessian of conditional posterior, found numerically
    res <- optim(log(pars$sigtau2), obj_sigtau2, gr=NULL, pars, method="L-BFGS-B", hessian=TRUE)
    lsigtau2_mean <- drop(res$par)
    lsigtau2_var <- 1/drop(res$hessian)
    lsigtau2_prop <- lsigtau2_mean + sqrt(lsigtau2_var)*rt(1, df=nu)
    lr_prop <- dt((log(pars$sigtau2) - lsigtau2_mean)/sqrt(lsigtau2_var), df=nu, log=TRUE) -
        dt((lsigtau2_prop - lsigtau2_mean)/sqrt(lsigtau2_var), df=nu, log=TRUE)
    sigtau2_prop <- exp(lsigtau2_prop)
    prop <- prop_sigtau2(sigtau2_prop, pars)
    if (!is.null(prop)) {
        lr_prior <- prop$lprior - lprior_current; prop$lprior <- NULL
        alpha[j, 4] <- min(exp(prop$llkP - pars$llkP + lr_prior + lr_prop), 1)
        if (runif(1) < alpha[j, 4]) {
            nms <- c("sigtau2", "Sigma.tilde", "llkP")
            ## stopifnot(setequal(names(prop), nms))
            pars[nms] <- prop[nms]
        }
    }

    ## ###########################################
    ## 6) Phi -- Gibbs
    Ptilde <- pars$cP - rep(1, nobs) %*% rbind(pars$Pbar) - pars$tau %*% rbind(pars$gamma)
    Ymat <- t(Ptilde[2:nobs,])
    Zmat <- t(Ptilde[1:(nobs-1),])
    Sigma_inv <- solve(pars$Sigma.tilde)
    phi_cov <- solve(prior$Omega_inv + kronecker(tcrossprod(Zmat), Sigma_inv))
    phi_mean <- phi_cov %*% (prior$Omega_inv %*% as.numeric(prior$Phi) + kronecker(Zmat, Sigma_inv) %*% as.numeric(Ymat))
    Phi_mean <- matrix(phi_mean, N, N)
    phi_draw <- phi_mean + t(chol(phi_cov)) %*% rnorm(N^2)
    Phi_draw <- matrix(phi_draw, N, N)
    ev <- max(abs(eigen(Phi_draw)$values))
    if (ev < 1) {
        pars$Phi <- Phi_draw
        pars$llkP <- llk_P(pars$tau, pars$cP, pars$Pbar, pars$gamma, pars$Phi, pars$Sigma.tilde, pars$sigtau2)
    }
    eigenval[j] <- ev

    ## ###########################################
    ## 7) sige2 - measurement error variance
    Yhat <- rep(1, nobs) %*% rbind(pars$AcP) + pars$cP %*% pars$BcP
    Q.errors <- Y - Yhat
    Q.errors[zlb.ind==1] <- 0  # ignore ZLB yields -- treated as missing!
    pars$sige2 <- rinvgamma(1, (prior$alpha_e + sum(zlb.ind==0))/2, (prior$delta_e + sum(Q.errors^2))/2)
    pars$llkQ <- -.5*sum(Q.errors^2)/pars$sige2

    ## save draws
    kinfQ[j] <- pars$kinfQ
    lamQ[j,] <- pars$lamQ
    gamma[j,] <- pars$gamma
    Pbar[j,] <- pars$Pbar
    sigtau2[j] <- pars$sigtau2
    sige2[j] <- pars$sige2
    Phi[j,] <- as.numeric(pars$Phi)
    Sigma[j,] <- pars$Sigma[ltri]
    Sigma.tilde[j,] <- pars$Sigma.tilde[ltri]
    AcP[j,] <- pars$AcP
    BcP[j,,] <- pars$BcP
    tau[j,] <- pars$tau
    cP[j,,] <- pars$cP
}

## drop burn-in sample
mcmc.ind <- seq(M/2+1, M, by=1)
kinfQ <- kinfQ[mcmc.ind]; lamQ <- lamQ[mcmc.ind,]; gamma <- gamma[mcmc.ind,]
sigtau2 <- sigtau2[mcmc.ind]; sige2 <- sige2[mcmc.ind]
Pbar <- Pbar[mcmc.ind,]; Phi <- Phi[mcmc.ind,]; eigenval <- eigenval[mcmc.ind]
Sigma.tilde <- Sigma.tilde[mcmc.ind,]; Sigma <- Sigma[mcmc.ind,]
AcP <- AcP[mcmc.ind,]; BcP <- BcP[mcmc.ind,,]
tau <- tau[mcmc.ind,]
cP <- cP[mcmc.ind,,]

save(data, zlb.ind, Y, WN, mats, yield.cols, dt, prior, kinfQ, lamQ, gamma, sigtau2, sige2, eigenval, Pbar, Phi, Sigma, Sigma.tilde, AcP, BcP, alpha, pars, pars.mode, tau, cP, mcmc.ind, pars2theta, makegamma, makePbar, file = filename, compress=TRUE)

