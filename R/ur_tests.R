## Functions for different unit root tests

persistence_stats <- function(z, nvar=1) {
    rho <- acf(z, lag.max=1, plot=FALSE)$acf[2]
    rval <- sprintf("%4.2f", c(sd(z), rho))
    rval <- c(rval, sprintf("%4.1f", log(0.5)/log(abs(rho))))
    rval <- c(rval, myADF(z, nvar), myPP(z, nvar)) #
    rval <- c(rval, sprintf("%4.2f", lfst(z, q=8)$pval))
    return(rval)
}

## Phillips-Perron Z_alpha Test
## based on urca::ur.pp
## changes:
## - only Z-alpha
## - no trend
## - allow for "none" -- AR(1) without intercept
my.ur.pp <- function(x, type=c("Z-alpha"), model=c("none", "constant"), lags=c("short", "long"), use.lag=NULL){
    x <- na.omit(as.vector(x))
    n <- length(x)
    y <- x[-1]
    y.l1 <- x[-n]
    n <- n-1
    lags <- match.arg(lags)
    model <- match.arg(model)
    type <- match.arg(type)
    if(!(is.null(use.lag))){
        lmax <- as.integer(use.lag)
        if(lmax < 0){
            warning("\nuse.lag has to be positive and integer; lags='short' used.")
            lmax <- trunc(4*(n/100)^0.25)}
    }else if(lags == "short"){
        lmax <- trunc(4*(n/100)^0.25)
    }else if(lags == "long"){
        lmax <- trunc(12*(n/100)^0.25)}
    if (model=="none") {
        test.reg <- summary(lm(y ~ y.l1 - 1))
        res <- residuals(test.reg)
        s <- 1/n*(sum(res^2))
        myy <- (1/n^2)*sum(y^2)
        idx <- 1:lmax
        coprods <- sapply(idx, function(l) t(res[-c(1:l)])%*%(res[-c((n-l+1):n)]))
        weights <- 1 - idx/(lmax+1)
        sig <- s + (2/n)*(t(weights)%*%coprods)
        alpha <- coef(test.reg)[1]
        teststat <- n*(alpha-1)-0.5*(sig-s)/myy
    } else if(model=="constant"){
        test.reg <- summary(lm(y ~ y.l1))
        res <- residuals(test.reg)
        s <- 1/n*(sum(res^2))
        myybar <- (1/n^2)*sum((y-mean(y))^2)
        idx <- 1:lmax
        coprods <- sapply(idx, function(l) t(res[-c(1:l)])%*%(res[-c((n-l+1):n)]))
        weights <- 1 - idx/(lmax+1)
        sig <- s + (2/n)*(t(weights)%*%coprods)
        alpha <- coef(test.reg)[2, 1]
        teststat <- n*(alpha-1)-0.5*(sig-s)/myybar
    }
    return(as.numeric(teststat))
}

cvsMacKinnon <- function(T, nvar) {
    ## MacKinnon critical values for 1%, 5%, 10%
    ## nvar=1 means no cointegrating regression, only one variable
    ## nvar=2 means residual from regression on one variable
    if (nvar==1) {
        b0 <- c(-3.43035, -2.86154, -2.56677)
        b1 <- c(-6.5393, -2.8903, -1.5384)
        b2 <- c(-16.786, -4.234, -2.809)
        b3 <- c(-79.433, -40.040, 0)
    } else if (nvar==2) {
        b0 <- c(-3.89644, -3.33613, -3.04445)
        b1 <- c(-10.9519, -6.1101, -4.2412)
        b2 <- c(-22.527, -6.823, -2.720)
        b3 <- rep(0, 3)
    } else if (nvar==3) {
        b0 <- c(-4.29374, -3.74066, -3.45218)
        b1 <- c(-14.4354, -8.5631, -6.2143)
        b2 <- c(-33.195, -10.852, -3.718)
        b3 <- c(47.433, 27.982, 0)
    } else {
        stop("not implemented for more than 3 variables in cointegration regression")
    }
    cvs <- sapply(1:3, function(i) b0[i] + b1[i]/T + b2[i]/T^2 + b3[i]/T^3)
    return(cvs)
}

myADF <- function(z, nvar) {
    ## Augmented Dickey-Fuller
    ## - general to specific procedure
    ## - use correct critical values
    kmax <- 4
    k <- kmax + 1; pval <- Inf
    type <- ifelse(nvar==1, "drift", "none") # no constant if series is cointegration residual
    while (pval > 0.1 & k > 0) {
        k <- k-1
        rval <- urca::ur.df(z, type, lags=k)
        pval <- tail(rval@testreg$coefficients[,4], 1)
    }
    result <- sprintf("%4.2f", rval@teststat[1])
    T <- length(rval@testreg$residuals)
    cvs <- cvsMacKinnon(T, nvar)
    ## print(cvs)
    sig <- (rval@teststat[1] < cvs)
    if (any(sig))
        result <- paste(c(result, rep("*", sum(sig))), collapse="")
    ## cat("ADF: lags = ", rval@lags, ", obs = ", T, ", N = ", nvar, "\n", sep="")
    return(result)
}

myPP <- function(z, nvar) {
    ## Phillips-Perron test
    if (nvar==1) {
        ## only one variable, no spurious regression
        rval <- my.ur.pp(z, type="Z-alpha", model="constant", lags="short")
        cvs <- c(-20.3, -14, -11.2) # Hamilton p. 762, T=250
    } else if (nvar>1) {
        rval <- my.ur.pp(z, type="Z-alpha", model="none", lags="short")
        if (nvar==2) {
            cvs <- c(-28.3, -20.5, -17) # Hamilton p. 765, T=500
        } else if (nvar==3) {
            cvs <- c(-34.2, -26.1, -22.2)
        }
    }
    result <- sprintf("%4.2f", rval)
    sig <- (rval < cvs)
    if (any(sig))
        result <- paste(c(result, rep("*", sum(sig))), collapse="")
    return(result)
}

##################################################
## Mueller-Watson low-frequency stationarity test

computePsi <- function(T, q) {
    Psi_j <- function(s, j)
        sqrt(2)*cos(j*s*pi)
    Psi <- function(s)
        Psi_j(s, j=1:q)
    t(vapply(((1:T)-1/2)/T, Psi, numeric(q)))
}

computeXs <- function(x, q) {
    ## standardized DCT coefficients
    T <- length(x)
    Psi_T <- computePsi(T, q)
    X_T <- t(Psi_T) %*% x / T # cosine-weighted averages
    X_T/sqrt(sum(X_T^2))
}

lfst <- function(x, q, Nsim=50000) {
    ## LFST test
    X <- computeXs(x, q)
    ## numerator (H0): I(0) model
    num <- sum(X^2)
    ## denominator (Ha): local-level model
    b <- 1/10
    chol_Ha_inv <- diag((1+1/(b*(1:q)*pi)^2)^(-1/2))
    Xa <- chol_Ha_inv %*% X
    denom <- sum(Xa^2)
    lfst <- num/denom
    Xsim <- matrix(rnorm(Nsim*q), q, Nsim)
    Xsimo <- Xsim
    Xsima <- chol_Ha_inv %*% Xsim
    lfst_sim <- colSums(Xsimo^2)/colSums(Xsima^2)
    pval <- mean(lfst_sim > lfst)
    list(stat=lfst, pval=pval)
}

