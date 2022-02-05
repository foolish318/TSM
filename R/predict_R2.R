## R^2 of predictive regressions for excess bond returns
## - alternative regression specifications
## - with small-sample bootstrap intervals for R^2 under the spanning null hypothesis
## Online Appendix Table 5

library(xtable) # xtable()
library(VAR.etp) # for bootstrap
source("R/data_fns.R")
source("R/bootstrap_fns.R")
source("R/util_fns.R")

data <- loadData()
data$istar <- data$pistar.ptr + data$rstar.realtime
yield.cols <- attr(data, "yield.cols")
mats <- as.numeric(substr(yield.cols, 2, 5))
Y <- data.matrix(data[,yield.cols])

M <- 5000

getR2 <- function(data) {
    sapply(list(lm(xr ~ PC1 + PC2 + PC3, data),
                lm(xr ~ PC1 + PC2 + PC3 + pistar.ptr, data),
                lm(xr ~ PC1 + PC2 + PC3 + pistar.ptr + rstar.realtime, data),
                lm(xr ~ PC1 + PC2 + PC3 + istar, data),
                lm(xr ~ PC1d1 + PC2d1 + PC3d1, data),
                lm(xr ~ PC1d2 + PC2d2 + PC3d2, data),
                lm(xr ~ PC1d3 + PC2d3 + PC3d3, data)),
           function(mod) summary(mod)$r.squared)
}

nms <- c("Yields only",
         "Yields and $\\pi^\\ast_t$",
         "Yields, $\\pi^\\ast_t$ and $r^\\ast_t$",
         "Yields and $i^\\ast_t$",
         "Yields detrended by $\\pi^\\ast_t$",
         "Yields detrended by $\\pi^\\ast_t$ and $r^\\ast_t$",
         "Yields detrended by $i^\\ast_t$")
N <- length(nms)  # number of specs
tbl <- data.frame(matrix(NA, N*2, 5))
tbl[,1] <- sprintf("%-45s", c(t(cbind(nms, ""))))
colnames(tbl) <- c("Specification", "Original sample R^2", "Delta R^2", "1985-2015 R^2", "Delta R^2")

h <- 1  ## one-quarter holding period
## h <- 4   ## annual holding period
xrn <- excess_returns(Y, mats, h) # construct excess returns
data$xr <- rowMeans(xrn)

data[paste0("PC", 1:3, "d1")] <- NA
data[paste0("PC", 1:3, "d2")] <- NA
data[paste0("PC", 1:3, "d3")] <- NA

col <- 2
for (subsample in c(FALSE, TRUE)) {
    if (subsample) {
        cat("Post-1985 subsample:\n")
        ind <- data$yyyymm >= 198501
    } else {
        cat("Full sample:\n")
        ind <- 1:nrow(data)
    }
    cat("Date range:", range(data$yyyymm), "\n")
    W <- eigen(cov(Y[ind,]))$vectors[,1:3]
    attr(data, "W") <- W
    attr(data, "mats") <- mats
    attr(data, 'yield.cols') <- yield.cols

    ## PCs of yields
    data[paste0("PC", 1:3)] <- Y %*% W # will only use data[ind,]

    ## PCs of detrended yields
    data[ind, paste0("PC", 1:3, "d1")] <- makePCs(lsfit(data$pistar.ptr[ind], Y[ind,])$residuals)
    data[ind, paste0("PC", 1:3, "d2")] <- makePCs(lsfit(cbind(data$pistar.ptr[ind], data$rstar.realtime[ind]), Y[ind,])$residuals)
    data[ind, paste0("PC", 1:3, "d3")] <- makePCs(lsfit(data$istar[ind], Y[ind,])$residuals)

    ## data
    R2 <- getR2(data[ind,])
    tbl[(1:N)*2-1, col] <- sprintf("%4.2f", R2)
    tbl[(2:N)*2-1, col+1] <- sprintf("%4.2f", R2[-1]-R2[1])
    tbl[1:2, col+1] <- ""

    ## bootstrapping
    dgp <- getBootDGP(c("PC1", "PC2", "PC3"), c("pistar.ptr", "rstar.realtime"), data, BC=TRUE)
    cat("# Simulating bootstrap samples: T =", length(ind), ", M =", M, "...\n")
    R2boot <- matrix(NA, M, N)
    for (b in 1:M) {
        simData <- simulateData(dgp, length(ind), h)
        simData$istar <- simData$pistar.ptr + simData$rstar.realtime
        Ysim <- as.matrix(simData[,yield.cols])
        simData[paste0("PC", 1:3)] <- makePCs(Ysim)
        simData[paste0("PC", 1:3, "d1")] <- makePCs(lsfit(simData$pistar.ptr, Ysim)$residuals)
        simData[paste0("PC", 1:3, "d2")] <- makePCs(lsfit(cbind(simData$pistar.ptr, simData$rstar.realtime), Ysim)$residuals)
        simData[paste0("PC", 1:3, "d3")] <- makePCs(lsfit(simData$istar, Ysim)$residuals)
        R2boot[b,] <- getR2(simData)
    }
    tmp <- apply(R2boot, 2, quantile, c(.025, .975))
    tbl[(1:N) * 2, col] <- apply(tmp, 2, function(x) sprintf("[%4.2f, %4.2f]", x[1], x[2]))
    dR2boot <- R2boot[,2:ncol(R2boot)] - R2boot[,1]
    tmp <- apply(dR2boot, 2, quantile, c(.025, .975))
    tbl[(2:N) * 2, col+1] <- apply(tmp, 2, function(x) sprintf("[%4.2f, %4.2f]", x[1], x[2]))

    col <- col+2
}

## construct and print table for paper
print(tbl)

filename <- "tables/returns_R2.tex"  ## quarterly returns
print(xtable(tbl),
      include.rownames=FALSE, include.colnames=FALSE, only.contents=TRUE,
      sanitize.text.function=function(x){x}, hline.after=8,
      file = filename)

