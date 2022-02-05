## functions for Dynamic Term Structure Models

## compiled functions for affine DTSM loadings
## this creates affineLoadingsCpp() and affineLoadingsCppBonly()
Rcpp::sourceCpp("R/affine.cpp") # this requires the package Rcpp and Rtools

C_PENAL <- 1e6 # a constant to penalize violations of constraints in likelihood estimation

##################################################
## JSZ = Joslin, Singleton, Zhu (2011, RFS)
## this is the restricted special case of our model with a Fixed Endpoint (FE)

obj_jsz <- function(theta, Y, WN, mats, dt) {
    ## objective function for numerical optimization
    rval <- C_PENAL
    try({
        N <- nrow(WN)
        pars <- theta2pars_jsz(theta, N)
        if (checkPars_jsz(pars)) # check validty of parameters
            res <- jsz.llk(Y, WN, K1Q.X=diag(pars$lamQ-1), Sigma.cP=pars$Omega.cP, mats=mats, dt=dt)
        rval <- sum(res$llk) # return sum of negative log likelihoods
    }, silent = TRUE)
    rval
}

checkPars_jsz <- function(pars) {
    ## check restrictions on parameters
    valid <- TRUE
    ## diagonal elements of Sigma positive and bounded away from zero
    if (any(diag(pars$L)<1e-7)) valid <- FALSE
    if (any(diag(pars$Omega.cP)<1e-10)) valid <- FALSE
    ## eigenvalues of Phi.Q not explosive
    if (any(pars$lamQ>1)) valid <- FALSE
    ## eigenvalues sorted in decreasing order
    if (any(pars$dlamQ>0)) valid <- FALSE
    ## P-eigenvalues
    if (!is.null(pars$Phi) && is.null(pars$Phi.LR)) {
        ## only need to check ev's if they weren't already checked for calculating Phi.LR
        maxev <- max(abs(eigen(pars$Phi)$values))
        stopifnot(all.equal(maxev, 1))
    }
    valid
}

theta2pars_jsz <- function(theta, N) {
    ## "unreparameterize"
    ## convert theta vector to list of individual parameters
    ## theta:  N + N*(N+1)/2    N=3  -> 3 + 6 = 9
    ## Q parameters
    pars <- list(dlamQ=theta[1:N])
    pars$lamQ=cumsum(pars$dlamQ+c(1, rep(0, N-1)))
    ## P-innovation covariance matrix
    pars$L <- matrix(0,N,N)
    pars$L[lower.tri(pars$L, diag=TRUE)] <- tail(theta, -N)
    pars$Omega.cP <- pars$L %*% t(pars$L)
    pars
}

pars2theta_jsz <- function(pars) {
    ## reparameterize to satisfy key parameter constraints
    ## convert list of individual parameters to theta vector
    dlamQ <- c(pars$lamQ[1]-1,diff(pars$lamQ));
    return(c(dlamQ, pars$L[lower.tri(pars$L, diag=TRUE)]))
}

getStartingValuesForMLE <- function(L, Y, WN, mats, dt, nSeeds = 500) {
    ## find good starting values for lamQ
    ## for MLE estimation of JSZ model with
    ##  (a) observable risk factors and
    ##  (b) unrestricted P-dynamics
    ## other parameters: L is given, kinfQ and P-dynamics concentrated from LLK
    ## Arguments: L - Cholesky of cov-mat, Y - yields, WN - loadings, mats - maturities in yields, dt - years per timestep (e.g., 1/12), nSeeds - how many random seeds to try
    ## Value: list with starting values
    Omega.cP <- L %*% t(L)
    ## use optimization
    N <- nrow(L)
    obj <- function(dlamQ) {
        lamQ <- cumsum(dlamQ + c(1, rep(0, N-1)))
        if (any(lamQ>1)) return(C_PENAL)
        if (any(dlamQ>0)) return(C_PENAL)
        if (any(lamQ<.5)) return(C_PENAL)
        res.llk <- jsz.llk(Y, WN, K1Q.X = diag(lamQ-1), Sigma.cP = Omega.cP, mats=mats, dt=dt)
        sum(res.llk$llk)
    }
    dlamQ <- c(-.001, -.1, -.2)
    cat("neg.llk at arbitrary dlamQ:", obj(dlamQ), "\n")
    rval <- optim(dlamQ, obj, method="L-BFGS-B")
    cat("neg.llk optimized:         ", rval$value, "\n")
    best.lamQ <- cumsum(rval$par + c(1, rep(0, N-1)))
    cat('LLK for best lamQ:', -rval$value, "\n")
    list(lamQ = best.lamQ, L = L)
}

estimate_jsz <- function(Y, est.sample, WN, mats, dt) {
    ## estimation of JSZ model
    ## 1) observable risk factors
    cP <- Y %*% t(WN)  # nobs x N
    N <- nrow(WN)
    ## 2) get some reasonable starting values for covariance matrix
    lm <- ar.ols(cP[est.sample,], aic=FALSE, order.max=1, intercept=TRUE, demean=FALSE)
    L <- t(chol(lm$var.pred))
    ## 3) get good starting values for lamQ
    pars.start <- getStartingValuesForMLE(L, Y[est.sample,], WN, mats, dt)
    theta.start <- pars2theta_jsz(pars.start)
    ## 4) numerical optimization  for cov mat and lamQ
    ## (rest of parameters are concentrated out)
    cat("Likelihood at starting values:", -obj_jsz(theta.start, Y[est.sample,], WN, mats, dt), "\n")
    theta <- get_optim(theta.start, obj_jsz, Y[est.sample,], WN, mats, dt, trace=0)
    pars <- theta2pars_jsz(theta, N)  ## lamQ, Omega
    cat("Likelihood at optimal values:", -obj_jsz(theta, Y[est.sample,], WN, mats, dt), "\n")
    ## 5) get all model parameters
    res.llk <- jsz.llk(Y[est.sample,], WN, K1Q.X=diag(pars$lamQ-1), Sigma.cP=pars$Omega.cP, mats=mats, dt=dt)
    mod <- within(c(pars, res.llk), {
        N <- N; dt <- dt
        cP <- Y %*% t(WN)  # nobs x N
        muQ <- K0Q.cP
        mu <- K0P.cP
        PhiQ <- K1Q.cP + diag(N)
        Phi <- K1P.cP + diag(N)
        loads.rn <- gaussian.loadings(round(mats/dt), K0P.cP, K1P.cP, Omega.cP, rho0.cP*dt, rho1.cP*dt, dt)
        Yhat <- rep(1, nrow(cP)) %*% AcP + cP %*% BcP
        Yrn <- rep(1, nrow(cP)) %*% loads.rn$A + cP %*% loads.rn$B
        Ytp <- Yhat - Yrn
        fhat <- 2*Yhat[,mats==10] - Yhat[,mats==5]
        frn <- 2*Yrn[,mats==10] - Yrn[,mats==5]
        ftp <- 2*Ytp[,mats==10] - Ytp[,mats==5]
    })
    cat("RMSE =", 10000*sqrt(mean((Y-mod$Yhat)^2)), "basis points\n")
    persistence(mod)
    mod
}

jsz.llk <- function (yields.o, W, K1Q.X, kinfQ=NA, K0P.cP=NA, K1P.cP=NA, Sigma.cP, mats, dt, sigma.e=NA) {
    ## Compute the likelihood for a Gaussian term structure.
    ## Source "A New Perspective on Gaussian Dynamic Term Structure Models" by Joslin, Singleton and Zhu
    ##
    ## INPUTS:
    ## yields.o   : (T+1)*J,  matrix of observed yields (first row are t=0 observations, which likelihood conditions on)
    ## mats       : 1*J,      maturities in years
    ## dt         : scalar,   length of period in years
    ## W          : N*J,      vector of portfolio weights to fit without error.
    ## K1Q,X      : N*N,      normalized latent-model matrix (does not have to be diagonal, see form below)
    ## Sigma.cP   : N*N,      positive definite matrix that is the covariance of innovations to cP
    ##
    ## OPTIONAL INPUTS -- concentrated out if not supplied:
    ## kinfQ      : scalar,   for stationary model, long run mean of (annualized) short rate under Q is -kinfQ/K1(m1,m1)
    ## K0P.cP     : N*1       intercept in VAR for cP
    ## K1P.cP     : N*N       mean reversion matrix in VAR for cP
    ## sigma.e    : scalar    standard error of yield observation errors
    ##
    ## OUTPUT:
    ## llk        : T*1       time series of -log likelihoods (includes 2-pi constants)
    ## AcP        : 1*J       yt = AcP' + BcP'*Xt  (yt is J*1 vector)
    ## BcP        : N*J       AcP, BcP satisfy internal consistency condition that AcP*W' = 0, BcP*W' = I_N
    ## AX         : 1*J       yt = AX' + BX'*Xt
    ## BX         : N*J       Xt is the 'jordan-normalized' latent state
    ## ...
    ##
    ## The model takes the form:
    ##   r(t) = rho0.cP + rho1.cP'*cPt
    ##        = 1'*Xt  (Xt is the 'jordan-normalized' state
    ##        = 1 period discount rate (annualized)
    ##
    ## Under Q:
    ##   X(t+1) - X(t)   = K0Q.X  + K1Q.X*X(t)  + eps_X(t+1),   cov(eps_X(t+1)) = Sigma.X
    ##   cP(t+1) - cP(t) = K0Q.cP + K1Q.cP*X(t) + eps_cP(t+1),  cov(eps_cP(t+1)) = Sigma.cP
    ##   where Sigma.X is chosen to match Sigma.cP
    ## and K0Q_X(m1) = kinfQ where m1 is the multiplicity of the highest eigenvalue (typically 1)

    ## Under P:
    ##   cP(t+1) - cP(t) = K0P.cP + K1P.cP*X(t) + eps_cP(t+1),  cov(eps_cP(t+1)) = Sigma.cP
    ##
    ## Model yields are given by:
    ##   yt^m = AcP' + BcP'*cPt  (J*1)
    ## And observed yields are given by:
    ##  yt^o = yt^m + epsilon.e(t)
    ## where V*epsilon.e~N(0,sigma.e^2 I_(J-N))
    ## and V is an (J-N)*J matrix which projects onto the span orthogonal to the
    ## row span of W.  This means errors are orthogonal to cPt and cPt^o = cPt^m.
    ##

    ## Setup
    T <- nrow(yields.o)-1
    J <- ncol(yields.o)
    N <- nrow(W)
    cP <- yields.o %*% t(W) # (T+1)*N, cP stands for math caligraphic P.

    ## COMPUTE THE Q-LIKELIHOOD:
    ## First find the loadings for the model:
    ## yt = AcP' + BcP'*cPt, AcP is 1*J, BcP is N*J
    if (is.na(kinfQ)) {
        ## concentrate out kinfQ
        ## AcP = alpha0_cP*kinf + alpha1_cP
        rho0.cP <- 0
        ## AcP0, AX0 will be the loadings with rho0_cP = 0, which won't be correct
        loads <- jsz.loadings.rho0cP(W, K1Q.X, rho0.cP, Sigma.cP, mats, dt)
        ## [BcP, AcP0, K0Q_cPx, K1Q_cP, rho0_cP, rho1_cP, K0Q_X, K1Q_X, AX0, BX, Sigma_X, alpha0_cP, alpha1_cP, alpha0_X, alpha1_X, m1]
        BcP <- loads$BcP; BX <- loads$BX
        alpha0.X <- loads$alpha0.X; alpha1.X <- loads$alpha1.X
        alpha0.cP <- loads$alpha0.cP; alpha1.cP <- loads$alpha1.cP
        K1Q.X <- loads$K1Q.X; m1 <- loads$m1;
        ## back out kinfQ that BEST fits average yields
        ## Note: V*e_t ~ N(0, sigma^2 I_{J-N})
        V <- t(MASS::Null(t(W)))
        kinfQ <- drop(t(colMeans(yields.o[2:(T+1),]) - t(alpha1.cP) - t(BcP)%*%colMeans(cP[2:(T+1),]))%*%(t(V)%*%V%*%t(alpha0.cP)) / (alpha0.cP%*%t(V)%*%V%*%t(alpha0.cP)))
        ## alternative to get least squares estimate   y = x*kinfQ + u
        ## y <- V %*% (colMeans(yields.o[2:(T+1),]) - t(alpha1.cP) - t(BcP)%*%colMeans(cP[2:(T+1),]))
        ## x <- V %*% t(alpha0.cP)
        ## crossprod(x, y)/crossprod(x, x)
        ## get correct loadings
        AX <- alpha0.X*kinfQ + alpha1.X
        AcP <- alpha0.cP*kinfQ + alpha1.cP
        ## get these to return to caller (not used in jsz.llk)
        K0Q.X <- matrix(0,N,1);
        K0Q.X[m1] <- kinfQ;
        params <- jsz.rotation(W, K1Q.X, K0Q.X, dt, BX, AX);
        K0Q.cP <- params$K0Q.cP; K1Q.cP <- params$K1Q.cP;
        rho0.cP <- params$rho0.cP; rho1.cP <- params$rho1.cP
    } else {
        loads <- jsz.loadings(W, K1Q.X, kinfQ, Sigma.cP, mats, dt)
        BcP <- loads$BcP; AcP <- loads$AcP; K0Q.cP <- loads$K0Q.cP; K1Q.cP <- loads$K1Q.cP;
        rho0.cP <- loads$rho0.cP; rho1.cP <- loads$rho1.cP; K0Q.X <- loads$K0Q.X; K1Q.X <- loads$K1Q.X;
        AX <- loads$AX; BX <- loads$BX;
    }
    yields.m <- rep(1,T+1)%*%AcP + cP %*% BcP # (T+1)*J, model-implied yields
    yield.errors <- yields.o[2:(T+1),] - yields.m[2:(T+1),]; # T*J
    square_orthogonal_yield.errors <- yield.errors^2; # T*J, but N-dimensional projection onto W is always 0, so effectively (J-N) dimensional
    ## Compute optimal sigma.e if it is not supplied
    if (is.na(sigma.e))
        sigma.e <- sqrt( sum(square_orthogonal_yield.errors)/(T*(J-N)) )
    ## Q-likelihood
    term1 <- .5*rowSums(square_orthogonal_yield.errors)/sigma.e^2
    term2 <- (J-N)*.5*log(2*pi)
    term3 <- .5*(J-N)*log(sigma.e^2)
    llkQ <- term1 + term2 + term3 # 1*T

    ## COMPUTE THE P-LIKELIHOOD:
    if (any(is.na(K0P.cP))||any(is.na(K1P.cP))) {
        ## Run OLS to obtain maximum likelihood estimates of K0P, K1P
        var1 <- ar.ols(cP, order=1, aic=FALSE, demean=FALSE, intercept=TRUE)
        Phi <- var1$ar[,,]
        mu <- var1$x.intercept
        res <- eigen(Phi)
        if (max(abs(res$values))>0.99) {
            Phi <- makeStationary(Phi)
            mu <- (diag(N) - Phi) %*% colMeans(cP)
        }
        K1P.cP <- Phi-diag(N)
        K0P.cP <- mu
    }
    innovations = t(cP[2:(T+1),]) - (K0P.cP%*%matrix(1,1,T) + (K1P.cP+diag(N))%*%t(cP[1:T,])) # N*T
    llkP = .5*N*log(2*pi) + .5*log(det(Sigma.cP)) + .5*colSums(innovations*solve(Sigma.cP, innovations)) # 1*T

    list(llk=t(llkQ + llkP), llkP=llkP, llkQ=llkQ, AcP=AcP, BcP=BcP, cP=cP, AX=AX, BX=BX, K0P.cP=K0P.cP, K1P.cP=K1P.cP, K0Q.cP=K0Q.cP, K1Q.cP=K1Q.cP, rho0.cP=rho0.cP, rho1.cP=rho1.cP, K0Q.X=K0Q.X, K1Q.X=K1Q.X, kinfQ=kinfQ, sigma.e=sigma.e, kinfQ=kinfQ)
}

jsz.loadings <- function(W, K1Q.X, kinfQ, Sigma.cP, mats, dt) {
    ## Inputs:
    ##   mats       : 1*J,      maturities in years
    ##   dt         : scalar,   length of period in years
    ##   W          : N*J,      vector of portfolio weights to fit without error.
    ##   K1Q.X      : N*N
    ##   kinfQ      : scalar,   determines long run mean
    ##   Sigma.cP   : N*N  covariance of innovations
    ##
    ## Returns:
    ##   AcP    : 1*J
    ##   BcP    : N*J
    ##   K0Q.cP : N*1
    ##   K1Q.cP : N*N
    ##   rho0.cP: scalar
    ##   rho1.cP: N*1
    ##   K0Q.X  : N*1
    ##   K1Q.X  : N*N
    ##   AX  : 1*J
    ##   BX  : N*J
    ##
    ## This function:
    ## 1. Compute the loadings for the normalized model:
    ##     X(t+1) - X(t) = K0Q.X + K1Q.X*X(t) + eps_X(t+1), cov(eps_X)=Sigma_X
    ##     and r(t) = 1.X(t)
    ##     where r(t) is the annualized short rate, (i.e. price of 1-period zero coupon bond at time t is exp(-r(t)*dt))
    ##    and K0Q_X(m1) = kinf, and K0Q_X is 0 in all other entries.
    ##      m1 is the multiplicity of the first eigenvalue.
    ##    Sigma.X is not provided -> solved for so that Sigma.cP (below) is matched.
    ##    yt = AX' + BX'*Xt
    ##
    ## 2. For cPt = W*yt and the model above for Xt, find AcP, BcP so that
    ##    yt = AcP' + BcP'*cPt
    ##
    ## 3. Computes the rotated model parameters K0Q.cP, K1Q.cP, rho0.cP, rho1.cP

    J <- length(mats)
    N <- nrow(K1Q.X)
    rho0d <- 0
    rho1d <- rep(1,N)
    mats.periods <- round(mats/dt)

    adjK1QX <- jszAdjustK1QX(K1Q.X)
    K1Q.X <- adjK1QX$K1Q.X
    m1 <- adjK1QX$m1
    ## m1 <- 1
    PhiQ <- K1Q.X + diag(N)

    ## we need to compute Sigma.X by first computing BX
    ##
    ## First compute the loadings ignoring the convexity term -- BX will be correct
    ## yt = AX' + BX'*Xt
    ## yt is J*1, AX is 1*J, BX is N*J, Xt is N*1, W is N*J
    ##
    ## cPt = W*yt  (cPt N*1, W is N*J)
    ##     = W*AX' + W*BX'*Xt
    ##     = WAXp + WBXp*Xt
    ##
    ## Substituting:
    ## yt = AX' + BX'*(WBXp\(cPt - WAXp))
    ##    = (I - BX'*(WBXp\WAXp))*AX' + BX'*WBXp\cPt
    ##    = AcP' + BcP'*cPt
    ## where AcP = AX*(I - BX'*(WBXp\WAXp))'
    ##       BcP = (WBXp)'\BX
    ##
    ## Sigma.cP = W*BX'*Sigma_X*(W*BX')'
    ## Sigma.X = (W*BX')\Sigma.cP/(W*BX')'

    ## loads.X.prelim <- gaussian.loadings(mats.periods, K0Q.X, K1Q.X, matrix(0, N, N), rho0d*dt, rho1d*dt, dt)
    loads.X.prelim <- affineLoadingsBonlyCpp(PhiQ, rho1d*dt, mats.periods, dt)
    BX <- loads.X.prelim$B
    WBXp <- W %*% t(BX)  # N*N
    WBXpinv <- solve(WBXp) # N*N
    WBXpinvp <- t(WBXpinv)
    Sigma.X <- WBXpinv %*% Sigma.cP %*% WBXpinvp # (W*BX')\Sigma.cP/(BX*W');

    ## Now with Sigma_X in hand, compute loadings for AX
    K0Q.X <- matrix(0,N,1)
    K0Q.X[m1] <- kinfQ
    ## loads.X <- gaussian.loadings(mats.periods, K0Q.X, K1Q.X, Sigma.X, rho0d*dt, rho1d*dt, dt)
    loads.X <- affineLoadingsCpp(K0Q.X, PhiQ, Sigma.X, rho0d*dt, rho1d*dt, mats.periods, dt)
    AX <- loads.X$A  # 1*J
    ## Rotate the model to obtain the AcP, BcP loadings.
    ## (See above for calculation)
    WAXp <- W %*% t(AX)  # N*1
    BcP <- WBXpinvp %*% BX # N*J
    AcP <- AX - AX %*% (t(W) %*% WBXpinvp %*% BX)
    ## compute rotated model parameters
    K1Q.cP <- WBXp %*% K1Q.X %*% WBXpinv
    K0Q.cP <- WBXp %*% K0Q.X - K1Q.cP %*% WAXp
    rho1.cP <- WBXpinvp %*% rep(1,N)
    rho0.cP <- -t(WAXp) %*% rho1.cP

    list(AX=AX, BX=BX, AcP=AcP, BcP=BcP, K0Q.cP=K0Q.cP, K1Q.cP=K1Q.cP, rho0.cP=rho0.cP, K0Q.X=K0Q.X, K1Q.X=K1Q.X, rho1.cP=rho1.cP, Sigma.X=Sigma.X)

}

jsz.loadings.prelim <- function(W, K1Q.X, mats, dt) {
    ## same as first part of jsz.loadings
    J <- length(mats)
    N <- nrow(K1Q.X)
    rho1d <- rep(1,N)
    mats.periods <- round(mats/dt)
    adjK1QX <- jszAdjustK1QX(K1Q.X)
    K1Q.X <- adjK1QX$K1Q.X
    m1 <- adjK1QX$m1
    ## m1 <- 1
    PhiQ <- K1Q.X + diag(N)
    loads.X.prelim <- affineLoadingsBonlyCpp(PhiQ, rho1d*dt, mats.periods, dt)
    BX <- loads.X.prelim$B
    WBXp <- W %*% t(BX)  # N*N
    WBXpinv <- solve(WBXp) # N*N
    WBXpinvp <- t(WBXpinv)
    BcP <- WBXpinvp %*% BX # N*J
    rho1.cP <- WBXpinvp %*% rep(1,N)
    list(BcP=BcP, rho1.cP=rho1.cP, BX=BX)
}

jsz.loadings.post <- function(W, BX, K1Q.X, kinfQ, Sigma.cP, mats, dt) {
    ## loadings if the preliminary step (BX etc) already done
    J <- length(mats)
    N <- nrow(K1Q.X)
    rho0d <- 0
    rho1d <- rep(1,N)
    mats.periods <- round(mats/dt)
    adjK1QX <- jszAdjustK1QX(K1Q.X)
    K1Q.X <- adjK1QX$K1Q.X
    m1 <- adjK1QX$m1
    PhiQ <- K1Q.X + diag(N)
    WBXp <- W %*% t(BX)  # N*N
    WBXpinv <- solve(WBXp) # N*N
    WBXpinvp <- t(WBXpinv)
    Sigma.X <- WBXpinv %*% Sigma.cP %*% WBXpinvp # (W*BX')\Sigma.cP/(BX*W');
    K0Q.X <- matrix(0,N,1)
    K0Q.X[m1] <- kinfQ
    loads.X <- affineLoadingsCpp(K0Q.X, PhiQ, Sigma.X, rho0d*dt, rho1d*dt, mats.periods, dt)
    AX <- loads.X$A  # 1*J
    WAXp <- W %*% t(AX)  # N*1
    AcP <- AX - AX %*% (t(W) %*% WBXpinvp %*% BX)
    rho1.cP <- WBXpinvp %*% rep(1,N)
    rho0.cP <- -t(WAXp) %*% rho1.cP
    list(AX=AX, AcP=AcP, rho0.cP=rho0.cP)
}

jsz.rotation <- function(W, K1Q.X, K0Q.X, dt, BX, AX) {
    ## Inputs:
    ##   W          : N*J,      vector of portfolio weights to fit without error.
    ##   K1Q.X      : N*N
    ##   K0Q.X      : N*1
    ##   dt         : scalar,   length of period in years
    ##   BX         : N*J  (BX, AX) are optional (saves time)
    ##   AX         : 1*J
    ##
    ## Returns:  [K0Q.cP, K1Q.cP, rho0.cP, rho1.cP]
    ##   K0Q.cP : N*1
    ##   K1Q.cP : N*N
    ##   rho0.cP : scalar
    ##   rho1.cP : N*1
    ##
    ## r(t) = rho0.cP + rho1.cP'*cPt
    ##      = 1'*Xt
    ##      = 1 period discount rate (annualized)
    ##
    ## Under Q:
    ##   X(t+1) - X(t)   = K0Q.X  + K1Q.X*X(t)  + eps_X(t+1),   cov(eps_X(t+1)) = Sigma_X
    ##   cP(t+1) - cP(t) = K0Q.cP + K1Q.cP*X(t) + eps_cP(t+1),  cov(eps_cP(t+1)) = Sigma.cP
    ## Where Sigma_X is chosen to match Sigma.cP
    ##
    ## cPt = W*yt  (cPt N*1, W is N*J)
    ##     = W*AX' + W*BX'*Xt
    ##     = WAXp + WBXp*Xt
    ##
    ## Delta(cP) = WBXp*Delta(Xt)
    ##           = WBXp*(K1Q.X*Xt + sqrt(Sigma_X)*eps(t+1))
    ##           = WBXp*(K1Q.X)*(WBXp\(cPt - WAXp)) + sqrt(Sigma.cP)*eps(t+1)
    ##           = WBXp*(K1Q.X)/WBXp*cPt - WBXp*(K1Q.X)/WBXp*WAXp] + sqrt(Sigma.cP)*eps(t+1)
    ##
    ## rt = 1'*Xt  [annualized 1-period rate]
    ##    = 1'*(WBXp\(cPt - WAXp))
    ##    = [- 1'*(WBXp\WAXp)] + ((WBXp)'1)'*cPt

    N <- nrow(K1Q.X)
    WBXp <- W %*% t(BX)
    WAXp <- W %*% t(AX)
    WBXpinv <- solve(WBXp)

    K1Q.cP <- WBXp %*% K1Q.X %*% WBXpinv
    K0Q.cP <- WBXp %*% K0Q.X - K1Q.cP %*% WAXp

    rho1.cP = t(WBXpinv) %*% rep(1,N)
    rho0.cP = - t(WAXp) %*% rho1.cP

    list(K0Q.cP=K0Q.cP, K1Q.cP=K1Q.cP, rho0.cP=rho0.cP, rho1.cP=rho1.cP)
}

gaussian.loadings <- function(maturities, K0d, K1d, H0d, rho0d, rho1d, timestep=1) {
    ## maturities: M*1
    ## K0d      : N*1
    ## K1d      : N*1
    ## H0d      : N*N
    ## rho0d    : scalar
    ## rho1d    : N*1
    ## timestep : optional argument.
    ##
    ## By : N*M
    ## Ay : 1*M  (faster to not compute with only one output argument)
    ##
    ## r(t)   = rho0d + rho1d'Xt
    ##        = 1 period discount rate
    ## P(t)   =  price of  t-period zero coupon bond
    ##        = EQ0[exp(-r0 - r1 - ... - r(t-1)]
    ##        = exp(A+B'X0)
    ## yields = Ay + By'*X0
    ##   yield is express on a per period basis unless timestep is provided.
    ##   --For example, if the price of a two-year zero is exp(-2*.06)=exp(-24*.005),
    ##   --and we have a monthly model, the function will return Ay+By*X0=.005
    ##   --unless timestep=1/12 is provided in which case it returns Ay+By*X0=.06
    ##
    ## Where under Q:
    ##   X(t+1) - X(t) = K0d + K1d*X(t) + eps(t+1),  cov(eps(t+1)) = H0d
    ##
    ## A1 = -rho0d
    ## B1 = -rho1d
    ## At = A(t-1) + K0d'*B(t-1) .5*B(t-1)'*H0d*B(t-1) - rho0d
    ## Bt = B(t-1) + K1d'*B(t-1) - rho1d
    ##
    ## maturities: 1*M # of periods

    M = length(maturities)
    N = length(K0d)
    Atemp = 0
    Btemp = matrix(0,N,1)
    Ay = matrix(NA,1,M)
    By = matrix(NA,N,M)

    curr_mat = 1
    K0dp <- t(K0d)
    K1dp <- t(K1d)
    for (i in 1:maturities[M]) {
        Atemp <- Atemp + K0dp%*%Btemp +.5%*%t(Btemp)%*%H0d%*%Btemp - rho0d
        Btemp <- Btemp + K1dp%*%Btemp - rho1d

        if (i==maturities[curr_mat]) {
            Ay[1,curr_mat] <- -Atemp/maturities[curr_mat]
            By[,curr_mat] <- -Btemp/maturities[curr_mat]
            curr_mat <- curr_mat + 1
        }
    }

    gaussian.loadings <- list(A = Ay/timestep, B = By/timestep)
}

jszAdjustK1QX <- function(K1Q.X, eps1=1e-3) {
    ## function [K1Q_X, isTypicalDiagonal, m1] = jszAdjustK1QX(K1Q_X, eps1);
    ##
    ## This function adjusts diagonal K1Q_X to give a non-diagonal but more
    ## computationally tractable K1Q_X.
    ##
    ## K1Q_X can fall into a few cases:
    ##   0. diagonal
    ##   1. not diagonal
    ##   2. zero eigenvalue
    ##   3. near repeated roots
    ## In cases 1-3, the diagonal closed form solver doesn't work, so compute differently.
    ## In case 1-2, we will use the recursive solver, though there are more efficient methods.
    ## In case 3, we will add a positive number above the diagonal.  this gives a different member of the set of observationally equivalent models.
    ##   So for example:
    ##      [lambda1, 0; 0, lambda2] is replaced by [lambda1, f(lambda2-lambda1); 0, lambda2] when abs(lambda1 - lambda2)<eps0
    ##   By making f not just 0/1, it will help by making the likelihood
    ##   continuous if we parameterize by kinf. (not an issue in some cases.)
    ##
    ## We also order the diagonal of diagonal K1Q.


    ## Cutoff function sets the super diagonal entry to something between 0 and
    ## 1, depending on how close the eigenvalues are.
    cutoff.fun <- function(x, eps1) {
        eps1 = 1e-3;
        eps0 = 1e-5;
        ##    xc <- 1*(x<eps0) + (1 - (x - eps0)/(eps1 - eps0))*(x>=eps0 && x<eps1) + 0*(x>eps1);
        xc <- 1*(log(x)<log(eps0)) +
            (1 - (log(x) - log(eps0))/(log(eps1) - log(eps0)))*(log(x)>=log(eps0) & log(x)<log(eps1)) +
            0*(log(x)>log(eps1));
        xc[x==0] <- 1;
        return(xc)
    }

    N <- nrow(K1Q.X)

    diag.K1Q.X <- diag(K1Q.X);
    isDiagonal <- all(K1Q.X==diag(diag.K1Q.X));

    ## For diagonal matrix, sort the diagonal and check to see if we have near repeated roots.
    if (isDiagonal) {
        diag.K1Q.X <- -sort(-diag.K1Q.X);
        K1Q.X <- diag(diag.K1Q.X);

        hasNearUnitRoot <- !all(abs(diag.K1Q.X)>eps1); ## Only applicable for diagonal
        hasNearRepeatedRoot <- !all(abs(diff(diag.K1Q.X))>eps1); ## Only applicable for diagonal
        isTypicalDiagonal <- isDiagonal && !hasNearRepeatedRoot && !hasNearUnitRoot;
    } else {
        isTypicalDiagonal <- FALSE
    }

    ## If we have a near repeated root, add a constnat above the diagonal. This
    ## representative of the equivalence class gives easier inversion for latent
    ## states vs. yields.  By varying the constant

    if (isDiagonal && !isTypicalDiagonal) {
        diff.diag <- abs(diff(diag.K1Q.X))
        super.diag <- cutoff.fun(diff.diag)
        K1Q.X[1:(N-1),2:N] <- K1Q.X[1:(N-1),2:N] +
            if (length(super.diag) == 1) {super.diag} else {diag(super.diag)}
    }

    super.diag = diag(K1Q.X[-N,-1]);
    m1 <- max(which(cumprod(c(1,super.diag))>0))

    list(K1Q.X=K1Q.X, isTypicalDiagonal=isTypicalDiagonal, m1=m1)
}

jsz.loadings.rho0cP <- function(W, K1Q.X, rho0.cP, Sigma.cP, mats, dt) {
    ## like jsz.loadings but parameterized in terms of rho0.cP instead of kinfQ
    ## gives not only loadings and rotate (cP) parameters,
    ## but also how AX and AcP load on kinfQ (alpha0.X, alpha1.X, alph0.cP, alpha1.cP)
    J <- length(mats)
    N <- nrow(K1Q.X)
    rho0d <- 0
    rho1d <- rep(1,N)
    mats.periods <- round(mats/dt)
    M <- max(mats.periods)

    adjK1QX <- jszAdjustK1QX(K1Q.X)
    K1Q.X <- adjK1QX$K1Q.X
    m1 <- adjK1QX$m1
    ## m1 <- 1
    PhiQ <- K1Q.X + diag(N)

    K0Q.X <- matrix(0,N,1)
    K0Q.X[m1] <- 1     # 1 instead of kinfQ

    ## loads.X.prelim <- gaussian.loadings(mats.periods, K0Q.X, K1Q.X, matrix(0, N, N), rho0d*dt, rho1d*dt, dt)
    loads.X.prelim <- affineLoadingsCpp(K0Q.X, PhiQ,  matrix(0, N, N), rho0d*dt, rho1d*dt, mats.periods, dt)
    BX <- loads.X.prelim$B
    alpha0.X <- loads.X.prelim$A    ### added
    WBXp <- W %*% t(BX)  # N*N
    WBXpinv <- solve(WBXp) # N*N
    rho1.cP <- t(WBXpinv) %*% rep(1,N)
    WBXpinvp <- t(WBXpinv)
    BcP <- WBXpinvp %*% BX # N*J
    Sigma.X <- WBXpinv %*% Sigma.cP %*% WBXpinvp # (W*BX')\Sigma.cP/(BX*W');

    ## Now with Sigma_X in hand, compute loadings for AX
    ## loads.X <- gaussian.loadings(mats.periods, K0Q.X, K1Q.X, Sigma.X, rho0d*dt, rho1d*dt, dt)
    loads.X <- affineLoadingsCpp(K0Q.X, PhiQ, Sigma.X, rho0d*dt, rho1d*dt, mats.periods, dt)
    AX1 <- loads.X$A  # 1*J         ### changed -- AX1 instead of AX
    ## AX1 gives the intercept with K0Q_X all zeros except 1 in the m1-th entry.
    ## So AX = alpha0_X*kinf + alpha1_X which alpha1_X = AX1 - alpha0_X
    alpha1.X <- AX1 - alpha0.X;     ### added

    ## Now find kinfQ that corresponds to desired rho0_cP:
    ## rt = 1'*Xt
    ## cPt = (W*alpha0_X')*kinf + (W*alpha1_X') + (W*BX')*Xt
    ## rt = 1'*(W*BX')^(-1)*[cPt -(W*alpha0_X')*kinf - (W*alpha1_X')]
    ## --> rho0_cP = -1'*(W*BX')^(-1)*(W*alpha0_X')*kinf - 1'*(W*BX')^(-1)*(W*alpha1_X')
    ## rho0_cP = beta0*kinfQ + beta1 => kinfQ = (rho0cP - beta1)/beta0
    beta0 <- -matrix(1,1,N)%*%WBXpinv%*%W%*%t(alpha0.X)
    beta1 <- -matrix(1,1,N)%*%WBXpinv%*%W%*%t(alpha1.X)
    kinfQ <- drop((rho0.cP - beta1)/beta0)
    ## kinfQ <- as.numeric(-(rho0.cP + a1)/a0);  ## -a0*kinf - a1 = rho0_cP
    ## kinfQ.alt <- drop((rho0.cP + t(rho1.cP) %*% W %*% t(alpha1.X))/(-t(rho1.cP) %*% W %*% t(alpha0.X)))
    ## stopifnot(all.equal(kinfQ.alt, kinfQ))
    K0Q.X[m1] <- kinfQ;

    ## now we can get correct A-loadings
    AX <- alpha0.X*kinfQ + alpha1.X
    C <- diag(J) - t(BX) %*% WBXpinv %*%W
    alpha0.cP <- t(C %*% t(alpha0.X))
    alpha1.cP <- t(C %*% t(alpha1.X))
    AcP <- alpha0.cP*kinfQ + alpha1.cP
    ## stopifnot(all.equal(AcP, AX %*% t(C)))

    ## compute rotated model parameters - as in jsz.loadings
    K1Q.cP <- WBXp %*% K1Q.X %*% WBXpinv
    WAXp <- W %*% t(AX)  # N*1
    K0Q.cP <- WBXp %*% K0Q.X - K1Q.cP %*% WAXp
    ## stopifnot(all.equal(drop(-t(WAXp) %*% rho1.cP), rho0.cP))
    list(AX=AX, BX=BX, AcP=AcP, BcP=BcP, K0Q.cP=K0Q.cP, K1Q.cP=K1Q.cP, rho0.cP=rho0.cP, K0Q.X=K0Q.X, K1Q.X=K1Q.X, rho1.cP=rho1.cP, alpha0.X=alpha0.X, alpha1.X=alpha1.X, alpha0.cP=alpha0.cP, alpha1.cP=alpha1.cP, m1=m1, beta0=beta0, beta1=beta1)
}

##################################################
## other helper functions

persistence <- function(mod) {
    ## VAR dynamics
    cat("P-eigenvalues:\n")
    print(abs(eigen(mod$Phi)$values))
    cat("Q-eigenvalues:\n")
    print(abs(eigen(mod$PhiQ)$values))
}

getYrn <- function(tau, Ptilde, Phi, sigtau2, Omega.Ptilde, rho0.cP, rho1.cP, mats, dt) {
    ## calculate risk-neutral yields
    ## risk factors: tau, Ptilde
    N <- nrow(Phi)
    mu.tilde <- rep(0, 1+N)
    Phi.tilde <- rbind(c(1, rep(0, N)),
                        cbind(0, Phi))
    Omega <- rbind(c(sigtau2, rep(0, N)),
                   cbind(0, Omega.Ptilde))
    delta.0 <- 0
    delta.1 <- c(1, rho1.cP)
    loads <- gaussian.loadings(round(mats/dt), mu.tilde, Phi.tilde - diag(N+1), Omega, delta.0*dt, delta.1*dt, dt)
    ## # without convexity
    ## loads <- gaussian.loadings(round(mats/dt), mu.tilde, Phi.tilde - diag(N+1), matrix(0, N+1, N+1), delta.0*dt, delta.1*dt, dt)
    Ztilde <- cbind(tau, Ptilde)
    rep(1, nrow(data)) %*% loads$A + Ztilde %*% loads$B
}

sim_unsp_model <- function(T, mod) {
    N <- nrow(mod$BcP)
    J <- ncol(mod$BcP)
    Ptilde <- matrix(0, T, N)
    tau <- numeric(T)
    utilde <- matrix(rnorm(T*N), T, N) %*% chol(mod$Sigma.tilde)
    eta <- rnorm(T, mean=0, sd=mod$sigma.tau)
    for (t in 2:T) {
        Ptilde[t,] <- mod$Phi %*% drop(Ptilde[t-1,]) + utilde[t,]
        tau[t] <- tau[t-1] + eta[t]
    }
    cP <- rep(1, T) %*% rbind(mod$Pbar) + tau %*% rbind(mod$gamma) + Ptilde
    Yhat <- rep(1, T) %*% mod$AcP + cP %*% mod$BcP
    Ysim <- Yhat + matrix(rnorm(T*J, mean=0, sd=mod$sigma.e), T, J)
    list(Y = Ysim, Ptilde = Ptilde, cP = cP, istar = tau)
}

sim_jsz_model <- function(T, mod) {
    N <- mod$N
    cP <- matrix(0, T, N)
    J <- ncol(mod$BcP)
    errors <- matrix(rnorm(T*N), T, N) %*% chol(mod$Omega.cP)
    for (t in 2:T)
        cP[t,] <- mod$mu + mod$Phi %*% drop(cP[t-1,]) + errors[t,]
    list(Y = rep(1, T) %*% mod$AcP + cP %*% mod$BcP + matrix(rnorm(T*J, 0, mod$sigma.e), T, J), cP = cP)
}

makePbar <- function(p, rho0, rho1) {
    N <- length(rho1)
    c(p, (-rho0 - crossprod(p, rho1[-N]))/rho1[N])
}

makegamma <- function(a, rho1) {
    ## gamma a.k.a A2
    N <- length(rho1)
    c(a, (1 - crossprod(a, rho1[-N]))/rho1[N])
}

check_restrictions <- function(P0, gamma, rho0, rho1) {
    a <- drop(rho0 + crossprod(rho1, P0))
    b <- drop(crossprod(rho1, gamma))
    if (!(isTRUE(all.equal(a, 0))))
        cat("not exactly zero:", a, "\n")
    if (!(isTRUE(all.equal(b, 1))))
        cat("not exactly one:", b, "\n")
    isTRUE(all.equal(c(a,b), c(0,1)))
}

##################################################
## OSE model

scale.ose <- c(10,10,1,1,1, 1, .1,.1,10,100,100,10,100,100,.1)
theta2pars_ose <- function(theta, WN, mats, dt) {
    ## convert theta vector to list of individual parameters
    N <- nrow(WN)
    len <- 1 + N + N-1 + N-1 + N*(N+1)/2 + 1
    ##  kinfQ + dlamQ + p + a + L + sigma.tau  -- Phi, sige2 concentrated out
    if (length(theta) != len)
        stop("length of theta is ", length(theta), " and not ", len, " as it should be")
    theta <- theta/scale.ose
    pars <- list(kinfQ = theta[1], dlamQ = theta[2:(N+1)]); theta <- tail(theta, -(N+1))
    pars$lamQ <- 1 + cumsum(pars$dlamQ)
    pars$p <- head(theta, N-1); theta <- tail(theta, -(N-1))
    pars$a <- head(theta, N-1); theta <- tail(theta, -(N-1))
    L <- matrix(0, N, N); L[lower.tri(L, diag=TRUE)] <- head(theta, 6); theta <- tail(theta, -6)
    pars$Sigma  <- L %*% t(L)
    stopifnot(length(theta)==1)
    pars$sigma.tau <- exp(theta[1])
    pars
}
pars2theta_ose <- function(pars) {
    ## convert list of individual parameters to theta vector
    dlamQ <- c(pars$lamQ[1]-1,diff(pars$lamQ))
    L <- t(chol(pars$Sigma))
    rval <- c(pars$kinfQ, dlamQ, pars$p, pars$a, L[lower.tri(L, diag=TRUE)], log(pars$sigma.tau))
    rval * scale.ose
}

startval_ose <- function(Y, istar, mod.jsz, WN, mats, dt) {
    cP <- Y %*% t(WN)
    ## VAR with intercept for Ptilde
    nobs <- nrow(cP); stopifnot(nobs==length(istar))
    sigtau2 <- var(diff(istar))
    obj <- function(theta) {
        p <- head(theta, N-1)
        a <- tail(theta, N-1)
        if (any(abs(p)>1)) return(C_PENAL)
        -sum(llk_ose(Y, istar, WN, mod.jsz$kinfQ, mod.jsz$lamQ, p, a, mod.jsz$Omega.cP, sigma.tau=sqrt(sigtau2), mats, dt)$llk)
    }
    ## optimal gamma if it wasn't constrained
    gamma <- rowSums(WN)
    Pbar <- colMeans(cP - istar %*% rbind(gamma))
    ## search for constrained optimal gamma
    theta.start <- c(Pbar[-N], gamma[-N])
    while (obj(theta.start) == C_PENAL)
        theta.start <- c(Pbar[-N], rnorm(2))
    res <- optim(theta.start, obj, method="L-BFGS-B")
    p <- head(res$par, N-1)
    a <- tail(res$par, N-1)
    Pbar <- makePbar(p, mod.jsz$rho0.cP, mod.jsz$rho1.cP)
    gamma <- makegamma(a, mod.jsz$rho1.cP)
    Ptilde <- cP - rep(1, nobs) %*% rbind(Pbar) - istar %*% rbind(gamma)
    res <- ar.ols(Ptilde, aic=FALSE, order.max=1, demean=FALSE, intercept=FALSE)
    Phi <- res$ar[,,]
    Sigma.tilde <- mod.jsz$Omega.cP - outer(gamma, gamma)*sigtau2
    if (invalidCovMat(Sigma.tilde)) {
        print(gamma)
        print(eigen(Sigme.tilde)$values)
        stop()
    }
    c(mod.jsz[c("kinfQ", "lamQ", "dlamQ", "sigma.e")], list(p=p, Pbar=Pbar, a=a, gamma=gamma, Sigma = mod.jsz$Omega.cP, Sigma.tilde = Sigma.tilde, sigtau2=sigtau2, Phi=Phi, sigma.tau = sqrt(sigtau2), sige2 = mod.jsz$sigma.e^2, Ptilde=Ptilde))
}

llk_ose <- function (Y, istar, W, kinfQ, lamQ, p, a, Sigma, sigma.tau, mats, dt, sigma.e=NA) {
    ## Compute likelihood for DTSM with unspanned shifting endpoint that is observed
    ## Setup
    T <- nrow(Y)-1
    J <- ncol(Y)
    N <- nrow(W)
    stopifnot(length(p)==N-1)
    stopifnot(length(a)==N-1)
    stopifnot(ncol(W)==J)
    stopifnot(nrow(Y) == length(istar))
    stopifnot(length(lamQ)==N)
    stopifnot(all.equal(dim(Sigma), c(N,N)))
    cP <- Y %*% t(W) # (T+1)*N

    ## loadings
    loads <- jsz.loadings(W, diag(lamQ-1), kinfQ, Sigma, mats, dt)

    ## Q-LIKELIHOOD:
    yields.m <- rep(1,T+1) %*% loads$AcP + cP %*% loads$BcP # (T+1)*J, model-implied yields
    yield.errors <- Y[2:(T+1),] - yields.m[2:(T+1),]; # T*J
    square_orthogonal_yield.errors <- yield.errors^2; # T*J
    ## Compute optimal sigma.e if it is not supplied
    if (is.na(sigma.e))
        sigma.e <- sqrt( sum(square_orthogonal_yield.errors)/(T*(J-N)) )
    llkQ <- -.5*(J-N)*log(2*pi*sigma.e^2) -.5*rowSums(square_orthogonal_yield.errors)/sigma.e^2

    ## P-LIKELIHOOD:
    gamma <- makegamma(a, loads$rho1.cP)
    Pbar <- makePbar(p, loads$rho0.cP, loads$rho1.cP)
    Sigma.tilde <- Sigma - outer(gamma, gamma) * sigma.tau^2
    if (invalidCovMat(Sigma.tilde)) return(list(llk=-C_PENAL))
    Ptilde <- cP - rep(1, T+1) %*% rbind(Pbar) - istar %*% rbind(gamma)
    var1 <- ar.ols(Ptilde, order=1, aic=FALSE, demean=FALSE, intercept=FALSE)
    Phi <- var1$ar[,,]
    if (max(abs(eigen(Phi)$values))>0.99)
        Phi <- makeStationary(Phi)
    innovations = t(Ptilde[2:(T+1),]) - Phi %*% t(Ptilde[1:T,]) # N*T
    llkP1 <- -.5*N*log(2*pi) - .5*log(det(Sigma.tilde)) - .5*colSums(innovations*solve(Sigma.tilde, innovations))
    llkP2 <- -.5*log(2*pi*sigma.tau^2) -.5*diff(istar)^2/sigma.tau^2
    stopifnot(length(llkP1)==T)
    stopifnot(length(llkP2)==T)
    llkP  <- llkP1 + llkP2

    list(llk = t(llkQ + llkP), Phi=Phi, sigma.e=sigma.e, llkP=llkP, llkQ=llkQ, cP=cP, Ptilde=Ptilde, AcP = loads$AcP, BcP = loads$BcP, Pbar = Pbar, gamma = gamma, istar = istar, Sigma.tilde = Sigma.tilde, rho0.cP = loads$rho0.cP, rho1.cP = loads$rho1.cP, K0Q.cP = loads$K0Q.cP, K1Q.cP = loads$K1Q.cP)
}

checkPars_ose <- function(pars) {
    ## check restrictions on parameter space
    valid <- TRUE
    if (pars$kinfQ < 0) valid <- FALSE
    ## eigenvalues of Phi.Q not explosive and decreasing
    if (any(pars$lamQ > 1)) valid <- FALSE
    if (any(pars$dlamQ > 0)) valid <- FALSE
    if (any(pars$lamQ < 0)) valid <- FALSE
    ## diagonal elements of Sigma positive and bounded away from zero
    if (any(diag(pars$L)<.Machine$double.eps^0.5)) valid <- FALSE
    if (any(diag(pars$Sigma)<.Machine$double.eps)) valid <- FALSE
    if (invalidCovMat(pars$Sigma)) valid <- FALSE
    if (pars$sigma.tau > 0.1) valid <- FALSE
    valid
}

obj_ose <- function(theta, istar, Y, WN, mats, dt) {
    ## objective function -- negative likelihood
    rval <- C_PENAL
    ## print(round(theta, 2))
    pars <- theta2pars_ose(theta, WN, mats, dt)
    ## check
    if (checkPars_ose(pars)) {
        rval <- -sum(llk_ose(Y, istar, WN, pars$kinfQ, pars$lamQ, pars$p, pars$a, pars$Sigma, pars$sigma.tau, mats, dt)$llk)
    }
    rval
}

