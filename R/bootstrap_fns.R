## functions for Bauer and Hamilton (2018) bootstrap for testing the spanning hypothesis

coef2array <- function(B, N, lag) {
    stopifnot(all.equal(dim(B), c(N, N * lag), check.attributes=FALSE))
    aperm(structure(B, dim = c(N, N, lag)), c(3, 1, 2))
}

coef2companion <- function(B) {
    N <- nrow(B)
    rbind(B, cbind(diag(ncol(B)-N), matrix(0, ncol(B)-N, N)))
}

getEVs <- function(B)
    eigen(coef2companion(B))$values

## following functions adapted from VAR.etp
## source code for VAR.etp package
## https://github.com/cran/VAR.etp/tree/master/R

VAR.ys2 <- function (x, b, p, e, type) {
    ## simulate y
    ## faster version of VAR.ys (not expanding matrix y but pre-allocating)
    n <- nrow(x)
    k <- nrow(b)
    b0 <- b[, 1:(k * p), drop = FALSE]
    if (type == "const")
        b1 <- b[, (k * p) + 1]
    if (type == "const+trend") {
        b1 <- b[, (k * p) + 1]
        b2 <- b[, (k * p) + 2]
    }
    y <- matrix(0, n, k)
    colnames(y) <- colnames(x)
    y[1:p,] <- x[1:p, , drop = FALSE] ## initial values are first few observations of actual series
    for (i in (p + 1):n) {
        index <- 1:k
        d1 <- 0
        for (j in 1:p) {
            d1 <- d1 + b0[, index] %*% y[i - j,]
            index <- index + k
        }
        d1 <- d1 + e[i - p, ]
        if (type == "const")
            d1 <- d1 + b1
        if (type == "const+trend")
            d1 <- d1 + b1 + b2 * i
        y[i, ] <- d1
    }
    return(y)
}

VAR.adjust2 <- function (b, bias, p, type) {
    require(VAR.etp)
    ## modified version of VAR.adjust -- throw error if NAs/infinite
    k <- nrow(b)
    bs1 <- b - bias
    delta <- VAR.etp:::VAR.modul(bs1, p)[1]
    if (delta < 1)
        bs2 <- bs1[, 1:(k * p), drop = FALSE]
    if (delta >= 1) {
        delta1 <- 1
        while (delta >= 1) {
            delta1 <- delta1 - 0.01
            bias <- delta1 * bias
            bs2 <- b[, 1:(k * p), drop = FALSE] - bias[, 1:(k *
                p), drop = FALSE]
            if (is.nan(sum(bs2)) | is.infinite(sum(bs2))) {
                print(bs2)
                stop("no bueno")
                ## bs2 <- b[, 1:(p * k), drop = FALSE]
                ## break
            }
            delta <- VAR.etp:::VAR.modul(bs2, p)[1]
        }
    }
    if (type == "const" | type == "const+trend")
        bs2 <- cbind(bs2, bs1[, (p * k + 1):ncol(b), drop = FALSE])
    return(bs2)
}

VAR.est2 <- function (x, p, type = "const", coef.only=FALSE) {
    ## much faster version than VAR.est
    n <- nrow(x)
    k <- ncol(x)
    y <- t(x[(p + 1):n, ])
    z <- matrix(0, k*p, n-p)
    xp <- t(x)
    for (i in (p+1):n)
        z[, i-p] <- xp[,(i-1):(i-p)]
    if (type == "const")
        z <- rbind(z, 1)
    if (type == "const+trend")
        z <- rbind(z, 1, (p + 1):n)
    b <- tcrossprod(y, z) %*% solve(tcrossprod(z))
    rownames(b) <- colnames(x)
    colnames(b) <- VAR.etp:::VAR.names(x, p, type)
    if (coef.only) {
        return(list(coef=b))
    } else {
        e <- y - b %*% z
        rownames(e) <- colnames(x)  # necessary for k=1
        sigu <- cov(t(e))  # tcrossprod(e)/((n - p) - ncol(b))
        zz <- tcrossprod(z)/(n - p)
        tem1 = (n - p)^(-1) * solve(zz) %x% sigu
        tem2 = sqrt(diag(tem1))
        tem3 = matrix(tem2, nrow = k, ncol = ncol(b))
        tmat = b/tem3
        return(list(coef = b, resid = t(e), sigu = sigu, zzmat = zz,
                    zmat = z, tratio = tmat, p=p))
    }
}

VAR.Boot2 <- function(x, p, nb = 500, type = "const", seed=1) {
    ## reverse engineer VAR.Boot
    set.seed(seed)
    require(VAR.etp)
    n <- nrow(x)
    k <- ncol(x)
    var1 <- VAR.est(x, p, type)
    b <- var1$coef  # OLS estimates
    mat <- matrix(0, nrow = k, ncol = ncol(b))
    for (i in 1:nb) {
        ## es <- VAR.etp:::resamp(e)
        es <- var1$resid[sample(n-p, n, replace=TRUE),,drop=FALSE]
        ## xs <- VAR.etp:::VAR.ys(x, b, p, es, type)
        xs <- VAR.ys2(x, b, p, es, type)
        bs <- VAR.est2(xs, p, type, coef.only=TRUE)$coef
        mat <- mat + bs/nb
    }
    bias <- mat - b
    bs <- VAR.adjust2(b, bias, p, type)
    colnames(bs) <- VAR.etp:::VAR.names(x, p, type)
    es <- VAR.etp:::VAR.resid(x, bs, var1$zmat, p)
    colnames(es) <- rownames(b)
    sigu <- cov(es)  # t(es) %*% es/((n - p) - ncol(b))
    return(list(coef = bs, resid = es, sigu = sigu, Bias = bias, p=p))
}

getYields <- function(df) {
    stopifnot("yield.cols" %in% names(attributes(df)))
    stopifnot(!any(is.na(match(attr(df, "yield.cols"), names(df)))))
    as.matrix(df[attr(df, "yield.cols")])
}

getBootDGP <- function(X1.names, X2.names, data, BC=FALSE, nb=1000) {
    ## prepare factors/variables and yield loadings
    X1 <- as.matrix(data[X1.names])
    X2 <- as.matrix(data[X2.names])
    if (!("mats" %in% names(attributes(data))))
        stop("data frame does not contain attribute 'mats'")
    if (!("W" %in% names(attributes(data))))
        stop("data frame does not contain attribute 'W'")
    W <- attr(data, "W")
    W <- W[,1:ncol(X1)]
    N1 <- ncol(X1); N2 <- ncol(X2)

    ## determine size of measurement error
    Yhat <- X1 %*% t(W)
    Y <- getYields(data)
    errors <- Y - Yhat
    sigma <- sqrt(mean(errors^2))

    dgp <- list(sigma=sigma, W=W, mats=attr(data, "mats"), X1.names=X1.names, X2.names=X2.names)
    dgp$X1.init <- X1[1,]
    dgp$X2.init <- X2[1,]

    cat("# DGP for Bootstrap \n")
    cat("VAR(1)\n")
    cat("orthogonal residuals\n")
    cat("Measurement error:", sigma, " percent\n")
    cat(ifelse(BC, "with", "without"), "bias-correction\n")

    ## VAR estimation - for X1 and X2 separately
    if (BC) {
        ## bootstrap bias correction
        rval1 <- VAR.Boot2(X1, 1, nb)
    } else {
        ## OLS
        rval1 <- VAR.est2(X1, 1)
    }
    dgp$B1 <- rval1$coef[, -ncol(rval1$coef), drop=FALSE]
    dgp$Phi1 <- coef2array(dgp$B1, N1, 1)
    dgp$resid1 <- rbind(matrix(NA, 1, N1), rval1$resid)
    dgp$mu1 <- rval1$coef[, ncol(rval1$coef), drop=FALSE]
    cat("max. abs. eigenvalue VAR for X1:", max(abs(getEVs(dgp$B1))), "\n")

    ## X2
    if (BC) {
        ## bootstrap bias correction
        rval2 <- VAR.Boot2(X2, 1, nb)
    } else {
        ## OLS
        rval2 <- VAR.est2(X2, 1)
    }
    dgp$B2 <- rval2$coef[, -ncol(rval2$coef), drop=FALSE]
    dgp$Phi2 <- coef2array(dgp$B2, N2, 1)
    dgp$resid2 <- rbind(matrix(NA, 1, N2), rval2$resid)
    dgp$mu2 <- rval2$coef[, ncol(rval2$coef), drop=FALSE]
    cat("max. abs. eigenvalue VAR for X2:", max(abs(getEVs(dgp$B2))), "\n")
    dgp
}

getReturn <- function(Y, n, h, mats) {
    ## n and mats in years
    ## h in quarters
    T <- nrow(Y)
    nmh <- ifelse(h==1, n, n-1) # approximate n-h - yield if quarterly change
    -(n-h/4)*Y[(1+h):T, mats==nmh] + n*Y[1:(T-h), mats==n] - h/4*Y[1:(T-h), mats==h/4]
}

simulateData <- function(dgp, T, h) {
    simVAR <- function(T, mu, Phi, e, X.init) {
        if (is.vector(e))
            dim(e) <- c(length(e), 1)
        N <- max(length(mu), length(mu))
        lags <- max(dim(Phi)[1], dim(Phi)[1])
        Tsim <- lags + T
        Xsim <- matrix(0, Tsim, N)
        Xsim[1,] <- X.init
        for (t in (lags+1):Tsim) {
            tmpsum <- mu + e[t,]
            for (j in 1:lags)
                tmpsum <- tmpsum + Phi[j,,] %*% Xsim[t-j,]
            Xsim[t,] <- tmpsum
        }
        Xsim[(lags+1-1):(lags+T-1),]
    }
    ## bootstrap residuals
    Tsim <- 1 + h + T
    resids1 <- na.omit(dgp$resid1)
    resids1 <- resids1[sample(1:nrow(resids1), Tsim, replace=TRUE),,drop=FALSE]
    resids2 <- na.omit(dgp$resid2)
    resids2 <- resids2[sample(1:nrow(resids2), Tsim, replace=TRUE),,drop=FALSE]
    ## (a) X1 and yields
    X1sim <- simVAR(T+h, dgp$mu1, dgp$Phi1, resids1, dgp$X1.init)
    Yhat <- X1sim %*% t(dgp$W)
    J <- length(dgp$mats)
    errors <- matrix(rnorm(J*(T+h), mean=0, sd=dgp$sigma), T+h, J)
    Ysim <- Yhat + errors
    colnames(Ysim) <- paste0("y", dgp$mats)
    ## predictors
    if (is.vector(X1sim)) {
        xdat1 <- as.matrix(X1sim[1:T])
    } else {
        xdat1 <- X1sim[1:T,]
    }
    colnames(xdat1) <- paste0("PC", 1:ncol(xdat1))
    ## (b) returns
    n <- dgp$mats[dgp$mats>1] ## for all maturities longer than one year
    returns <- sapply(n, function(i) getReturn(Ysim, i, h, dgp$mats))
    colnames(returns) <- paste0("xr", n)
    ## put together data set
    df <- data.frame(returns, xdat1, Ysim[1:T,])
    xrname <- "xr"
    df[[xrname]] <- rowMeans(returns)
    ## (c) X2
    X2sim <- simVAR(T, dgp$mu2, dgp$Phi2, resids2, dgp$X2.init)
    if (is.vector(X2sim)) {
        xdat2 <- as.matrix(X2sim[1:T])
    } else {
        xdat2 <- X2sim[1:T,]
    }
    colnames(xdat2) <- dgp$X2.names
    df <- data.frame(df, xdat2)
    df
}

bootstrapTest <- function(fmla1, fmla2, data, dgp, h, M=1000, vcovfn, adjR2=FALSE) {
    require(sandwich)
    require(lmtest)
    regnames <- attr(terms(fmla2), "term.labels")
    K <- length(regnames)
    K1 <- length(attr(terms(fmla1), "term.labels"))
    K2 <- K - K1
    depvar <- deparse(fmla1[[2]])
    R2fn <- if (adjR2) {
        function(mod) summary(mod)$adj.r.squared
    } else {
        function(mod) summary(mod)$r.squared
    }

    ## tables with results
    tblCoef <- matrix(NA, 5, K + 1)  ## 3 rows for data, 2 rows for bootstrap
    colnames(tblCoef) <- c(regnames, "Wald")
    rownames(tblCoef) <- c("Coefficient", "HAC statistic", "HAC $p$-value",
                           "Bootstrap 5\\% c.v.", "Bootstrap $p$-value")
    tblR2 <- matrix(NA, 3, 6)
    rownames(tblR2) <- c("Data", "Bootstrap mean", "95% bootstrap interval")
    colnames(tblR2) <- c("R^2_1", "", "R^2_2", "", "R^2_2 - R^2_1", "")

    ## data
    lm1 <- lm(fmla1, data=data)
    lm2 <- lm(fmla2, data=data)
    T <- length(lm2$residuals) + ifelse(any(grep("PC", depvar)), 1, 0)
    tblCoef[1, 1:K] <- lm2$coef[-1]
    ## HAC inference
    SEs <- sqrt(diag(vcovfn(lm2)))
    tstats <- (lm2$coef/SEs)[-1]
    tblCoef[2, 1:K] <- tstats
    tblCoef[3, 1:K] <- pnorm(abs(tstats), lower.tail=FALSE)*2
    rval <- waldtest(lm1, lm2, vcov=vcovfn, test="Chisq")
    stopifnot(K2 == rval$Df[2])  # check degrees of freedom
    tblCoef[2, K+1] <- rval$Chi[2]
    tblCoef[3, K+1] <- rval$Pr[2]
    ## R^2
    tblR2[1, c(1,3,5)] <- c(R2fn(lm1), R2fn(lm2), R2fn(lm2) - R2fn(lm1))

    ## bootstrapping
    cat("# Simulating bootstrap samples: T =", T, ", M =", M, "...\n")
    statsHAC <- matrix(NA, M, K+1)
    pvalsHAC <- matrix(NA, M, K+1)
    R2.r <- numeric(M)
    R2.ur <- numeric(M)
    for (b in 1:M) {
        simData <- simulateData(dgp, T, h)
        mod1 <- lm(fmla1, data=simData)
        mod2 <- lm(fmla2, data=simData)
        R2.r[b] <- R2fn(mod1)
        R2.ur[b] <- R2fn(mod2)
        ## HAC
        SEs <- sqrt(diag(vcovfn(mod2)))
        statsHAC[b, 1:K] <- abs(mod2$coef/SEs)[-1]
        pvalsHAC[b, 1:K] <- pnorm(statsHAC[b, 1:K], lower.tail=FALSE)*2
        rval <- waldtest(mod1, mod2, vcov=vcovfn, test="Chisq")
        statsHAC[b, K+1] <- rval$Chi[2]
        pvalsHAC[b, K+1] <- rval$Pr[2]
    }
    if (lm2$df != mod2$df)
        stop("degrees of freedom not the same in actual and simulated data")
    tblCoef[4, ] <- apply(statsHAC, 2, quantile, .95)  ## bootstrap critical values
    tblCoef[5, ] <- colMeans(statsHAC > rep(abs(tblCoef[2, ]), each=M))
    tblCoef[4:5, 1:K1] <- NA
    tblR2[2,c(1,3,5)] <- colMeans(cbind(R2.r, R2.ur, R2.ur-R2.r))      ## means
    tblR2[3, 1:2] <- quantile(R2.r, c(.025, .975))
    tblR2[3, 3:4] <- quantile(R2.ur, c(.025, .975))
    tblR2[3, 5:6] <- quantile(R2.ur-R2.r, c(.025, .975))
    list(tblCoef=tblCoef, tblR2=tblR2)
}
