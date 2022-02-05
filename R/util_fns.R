## utilities functions

myDMtest <- function (e1, e2, alternative = c("two.sided", "less", "greater"),
                      h = 1, power = 2, lrvar = c("HH", "NW")) {
    ## Diebold-Mariano test
    ## based on dm.test from forecast package
    alternative <- match.arg(alternative)
    d <- c(abs(e1))^power - c(abs(e2))^power
    lrvar <- match.arg(lrvar)
    if (lrvar == "HH") {
        ## Hansen-Hodrick
        d.cov <- acf(d, na.action = na.omit, lag.max = h - 1, type = "covariance",
                     plot = FALSE)$acf[, , 1]
        dv <- sum(c(d.cov[1], 2 * d.cov[-1]))/length(d)
    } else {
        ## Newey-West
        L <- floor(1.5*(h-1)) # need to use 1.5 x more lags, instead of L = h-1
        d.cov <- acf(d, na.action = na.omit, lag.max = L, type = "covariance",
                     plot = FALSE)$acf[, , 1]
        dv <- sum(c(d.cov[1], 2*(1 - (1:L)/(L+1)) * d.cov[-1]))/length(d)
        ## this is identical to NeweyWest(lm(d ~ 1), lag=L, prewhite=FALSE)
    }
    if (dv > 0)
        STATISTIC <- mean(d, na.rm = TRUE)/sqrt(dv)
    else if (h == 1)
        stop("Variance of DM statistic is zero")
    else {
        cat("Variance for h =", h, "is negative, using horizon h=", h-1, "\n")
        return(myDMtest(e1, e2, alternative, h = h-1, power))
        ## warning("Variance is negative, using NW")
        ## return(myDMtest(e1, e2, alternative, h, power, lrvar="NW"))
    }
    n <- length(d)
    k <- ((n + 1 - 2 * h + (h/n) * (h - 1))/n)^(1/2)
    STATISTIC <- STATISTIC * k
    names(STATISTIC) <- "DM"
    if (alternative == "two.sided")
        PVAL <- 2 * pt(-abs(STATISTIC), df = n - 1)
    else if (alternative == "less")
        PVAL <- pt(STATISTIC, df = n - 1)
    else if (alternative == "greater")
        PVAL <- pt(STATISTIC, df = n - 1, lower.tail = FALSE)
    PARAMETER <- c(h, power)
    names(PARAMETER) <- c("Forecast horizon", "Loss function power")
    structure(list(statistic = STATISTIC, parameter = PARAMETER,
        alternative = alternative, p.value = PVAL, method = "Diebold-Mariano Test",
        data.name = c(deparse(substitute(e1)), deparse(substitute(e2)))),
        class = "htest")
}

drawNormal <- function(mu, Omega)
    ## draw from multivariate normal
    mu + t(chol(Omega)) %*% rnorm(length(mu))

get_optim <- function(theta, obj, ..., trace=0) {
    ## more reliable numerical optimization
    obj <- match.fun(obj)
    val <- obj(theta, ...)
    if (trace>0)
        cat('Value at starting point:', val, '\n')
    i <- 1; improvement <- -Inf
    while (improvement < -.1) {
        res <- optim(theta, obj, gr=NULL, ..., control=list(maxit=5000))
        improvement <- res$value - val
        val <- res$value
        theta <- res$par
        if (trace>0)
            cat('iteration ', i,', value = ', val,'\n')
        i <- i + 1
    }
    if (trace>0)
        cat('improvement = ', improvement, ' -- proceed to final step\n')
    res <- optim(theta, obj, gr=NULL, ..., control=list(maxit = 50000))
    theta <- res$par
    val <- res$value
    if (trace>0) {
        cat('final Nelder-Mead step, value = ', res$value, "\n")
        cat("Convergence:", res$convergence, "\n")
        cat("Message:", res$message, "\n")
    }
    ## now iterate gradient-based until no more improvement
    improvement <- -Inf; i <- 1
    while (improvement < -.1) {
        improvement <- tryCatch({
            res <- optim(theta, obj, gr=NULL, ..., method="L-BFGS-B")
            improvement <- res$value - val
            val <- res$value
            theta <- res$par
            if (trace>0) {
                cat('Gradient-based, iteration ', i,', value = ', val,'\n')
                cat("Convergence:", res$convergence, "\n")
                cat("Message:", res$message, "\n")
            }
            improvement
        }, error = function(err) {
            if (trace>0)
                cat("Error in get_optim, L-BFGS-B:", err$message, "\n")
            return(0)
        })
        i <- i+1
    }
    res$par
}

logdinvgamma <- function(x, alpha, beta)
    alpha * log(beta) - lgamma(alpha) - (alpha + 1) * log(x) - (beta/x)

matrix.power <- function(X, a) {
    rval <- eigen(X)
    rval$vectors %*% diag(rval$values^a) %*% solve(rval$vectors)
}

makePD <- function(A) {
    ## make matrix positive definite
    if (min(eigen(A)$values)<0) {
        cat('matrix not PD, making adjustment...\n')
        D <- diag(eigen(A)$values)
        V <- eigen(A)$vectors
        D[D<0] <- 0+1e-6
        return(V %*% D %*% t(V))
    } else {
        return(A)
    }
}

plot_recessions <- function(date_range, yrange) {
    recessions <- rbind(c(196912, 197011),
                        c(197311, 197503),
                        c(198001, 198007),
                        c(198107, 198211),
                        c(199007, 199103),
                        c(200103, 200111),
                        c(200712, 200906))
    makeDate <- function(yyyymm) {
        s <- as.character(yyyymm*100+1)
        as.Date(s, "%Y%m%d")
    }
    for (i in which(recessions[,2] > min(as.numeric(format(date_range, "%Y%m"))))) {
        fromDate <- makeDate(recessions[i, 1])
        toDate <- makeDate(recessions[i, 2])
        polygon(x = c(fromDate, fromDate, toDate, toDate),
                y = c(yrange, rev(yrange)),
                density=NA, col=adjustcolor("gray",alpha.f=0.5), border=NA)
    }
}

invalidCovMat <- function(Omega) {
    ## is symmetric matrix not positive definite? then it's not a valid covariance matrix
    ## vals <- .Internal(La_rg(Omega, TRUE))$values # this is faster, knowing that it's real and symmetric
    vals <- eigen(Omega, symmetric=TRUE, only.values=TRUE)$values
    min(vals) < .Machine$double.eps^0.5
}

makeCovMat <- function(omega, N=3) {
    ## make symmetric covariance matrix from lower triangular elements
    ltri <- lower.tri(matrix(0, N, N), diag=TRUE)
    stopifnot(length(omega)==N*(N+1)/2)
    Om <- matrix(0, N, N)
    Om[ltri] <- omega
    Om[upper.tri(Om)] <- t(Om)[upper.tri(Om)]
    Om
}

makeStationary <- function(Phi, step=0.001, max.eigen=0.99) {
    ## make VAR mean-reversion matrix stationary -- ensure largest eigenvalue does not exceed 0.99
    ## like Kilian (1998, REStat)
    delta <- 1
    Phi.new <- Phi
    while (max(abs(eigen(Phi.new)$values))>max.eigen) {
        delta <- delta - step
        Phi.new <- delta * Phi
    }
    return(Phi.new)
}

makePCs <- function(Y) {
    W <- eigen(cov(Y))$vectors[,1:3]
    W <- W %*% diag(sign(W[nrow(W),])) # make sure PC1 and PC2 correspond to level and slope
    PCs <- Y %*% W
    attr(PCs, "W") <- W
    PCs
}

excess_returns <- function(Y, mats, h=4) {
    ## calculate excess bond returns in quarterly yield data
    ## Y - annualized yields (if in percent, returns will be in percent)
    ## mats - maturities in years
    stopifnot(h==4 | h==1)
    stopifnot(ncol(Y)==length(mats))
    nobs <- nrow(Y)
    xrn <- matrix(NA, nobs, sum(mats>1))
    for (n in mats[mats>1]) {
        nmh <- ifelse(h==4, n-1, n)  ## if h=1, approximate n-1/4 year bond yield with n-year bond yield
        xrn[, which(mats[mats>1] == n)] <- c(-(n-h/4)*Y[(1+h):nobs, mats==nmh] + n*Y[1:(nobs-h), mats==n] - h/4*Y[1:(nobs-h), mats==h/4], rep(NA, h))
    }
    xrn
}

predictReturns <- function(Y, istar, h) {
    ## take yields and i-star and return restricted and unrestricted R^2
    stopifnot(h==4 | h==1)
    xrn <- excess_returns(Y, mats, h)
    xr <- rowMeans(xrn)
    PC <- makePCs(Y)
    mod1 <- lm(xr ~ PC)
    mod2 <- lm(xr ~ PC + istar)
    c(summary(mod1)$r.squared, summary(mod2)$r.squared)
}

logrationorm <- function(x1, x2, mu, Sigma)
    ## log ratio of normal densities evaluated at two different vectors
    -0.5*crossprod(x1-mu, solve(Sigma, x1-mu)) + 0.5*crossprod(x2-mu, solve(Sigma, x2-mu))

