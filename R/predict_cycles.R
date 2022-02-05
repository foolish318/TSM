## Predict excess returns with PCs of detrended yields

library(xtable) # xtable()
library(VAR.etp) # for bootstrap
source("R/data_fns.R")
source("R/util_fns.R")

data <- loadData()
data$istar <- data$pistar.ptr + data$rstar.realtime
yield.cols <- attr(data, "yield.cols")
mats <- as.numeric(substr(yield.cols, 2, 5))
Y <- data.matrix(data[,yield.cols])

## h <- 4
## vcovfn <- function(mod) sandwich::NeweyWest(mod, lag=6, prewhite=FALSE)

h <- 1
vcovfn <- function(mod) sandwich::vcovHC(mod, type="HC0")

xrn <- excess_returns(Y, mats, h) # construct annual excess returns
data$xr <- rowMeans(xrn)
data[paste0("PC", 1:3)] <- NA

tbl <- data.frame(matrix(NA, 7, 9))
row.nms <- c("PC1", "", "PC2", "", "PC3", "", "$R^2$")
tbl[,1] <- sprintf("%-1s", row.nms)
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

    for (j in 1:4) {
        if (j==1) {
            PCs <- makePCs(Y[ind,])
        } else if (j==2) {   ## pi-star only
            PCs <- makePCs(lsfit(data$pistar.ptr[ind], Y[ind,])$residuals)
        } else if (j==3) {   ## pi-star and r-star
            PCs <- makePCs(lsfit(cbind(data$pistar.ptr[ind], data$rstar.realtime[ind]), Y[ind,])$residuals)
        } else if (j==4) {   ## i-star
            PCs <- makePCs(lsfit(data$istar[ind], Y[ind,])$residuals)
        }
        W <- attr(PCs, "W")
        sc <- c(sum(W[,1]), W[nrow(W),2]-W[1,2], W[nrow(W),3]-2*W[mats==2,3]+W[1,3])
        data[ind, paste0("PC", 1:3)] <- PCs %*% diag(1/sc)
        mod <- lm(xr ~ PC1 + PC2 + PC3, data[ind,])
        tbl[c(1,3,5), col] <- sprintf("%4.2f", mod$coef[-1])
        tbl[c(2,4,6), col] <- sprintf("(%4.2f)", sqrt(diag(vcovfn(mod)))[-1])
        tbl[7, col] <- sprintf("%4.2f", summary(mod)$r.squared)
        col <- col+1
    }
}

print(tbl)
filename <- "tables/returns_cycles.tex"
print(xtable(tbl, digi=2),
      include.rownames=FALSE, include.colnames=FALSE, only.contents=TRUE,
      sanitize.text.function=function(x){x}, hline.after=NULL,
      file = filename)
