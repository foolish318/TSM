## Predict excess returns with yields and trend proxies
## Table 2: predictive regressions with pi-star, r-star, i-star

library(xtable) # xtable()
library(VAR.etp) # for bootstrap
source("R/data_fns.R")
source("R/util_fns.R")
source("R/bootstrap_fns.R")

bootstrap_flag <- TRUE
M <- 5000

data <- loadData()
data$istar <- data$pistar.ptr + data$rstar.realtime
yield.cols <- attr(data, "yield.cols")
mats <- as.numeric(substr(yield.cols, 2, 5))
Y <- data.matrix(data[,yield.cols])

rstar.nms <- c("rstar.filt", "rstar.realtime", "rr.ewma")
rstar.desc <- c("filtered", "real-time", "mov.~avg.")

row.nms <- c("PC1", "", "PC2", "", "PC3", "",
             "$\\pi_t^\\ast$", "", "",
             "$r_t^\\ast$", "", "",
             "$i_t^\\ast$", "", "",
             "$R^2$",
             "Memo: $r^\\ast$")

## annual returns
## h <- 4
## vcovfn <- function(mod) sandwich::NeweyWest(mod, lag=6, prewhite=FALSE)

## quarterly returns
h <- 1
vcovfn <- function(mod) sandwich::vcovHC(mod, type="HC0")

xrn <- excess_returns(Y, mats, h)
data$xr <- rowMeans(xrn)

for (subsample in c(FALSE, TRUE)) {
    if (subsample) {
        cat("Post-1985 subsample:\n")
        ind <- data$yyyymm >= 198501
    } else {
        cat("Full sample:\n")
        ind <- 1:nrow(data)
    }

    ## create PCs
    PCs <- makePCs(Y)
    data[paste0("PC", 1:3)] <- PCs
    W <- attr(PCs, "W")
    ## scaled PCs (can't use in bootstrap where orthonormal vector are needed)
    sc <- c(sum(W[,1]), W[nrow(W),2]-W[1,2], W[nrow(W),3]-2*W[mats==2,3]+W[1,3])
    data[paste0("PC", 1:3, "sc")] <- PCs %*% diag(1/sc)
    ## for bootstrap
    attr(data, 'W') <- W
    attr(data, 'mats') <- mats
    attr(data, 'yield.cols') <- yield.cols
    fmla.0 <- xr ~ PC1 + PC2 + PC3

    ## table for results
    tbl <- data.frame(matrix(NA, length(row.nms), 7))
    tbl[,1] <- sprintf("%-15s", row.nms)
    col <- 2
    for (i in 1:(3+length(rstar.nms))) {
        if (i==1) {
            fmla <- xr ~ PC1sc + PC2sc + PC3sc
            regrows <- c(1,3,5)
            tbl[c(7:15, 17), col] <- ""
        } else if (i==2) {
            ## pi-star only
            fmla <- xr ~ PC1sc + PC2sc + PC3sc + pistar.ptr
            regrows <- c(1,3,5,7)
            if (bootstrap_flag) {
                fmla.a <- xr ~ PC1 + PC2 + PC3 + pistar.ptr
                dgp <- getBootDGP(c("PC1", "PC2", "PC3"), "pistar.ptr", data[ind,], BC=TRUE)
                rval <- bootstrapTest(fmla.0, fmla.a, data[ind,], dgp, h, M, vcovfn)
                tbl[9, col] <- sprintf("[%4.2f]", rval$tblCoef[5, 4])
            }
            tbl[c(10:15, 17), col] <- ""
        } else if ((i>2) & (i <= 2+length(rstar.nms))) {
            ## pi-star and r-star
            fmla <- formula(paste0("xr ~ PC1sc + PC2sc + PC3sc + pistar.ptr + ", rstar.nms[i-2]))
            regrows <- c(1,3,5,7,10)
            tbl[nrow(tbl), col] <- rstar.desc[i-2]
            if (bootstrap_flag) {
                fmla.a <- formula(paste0("xr ~ PC1 + PC2 + PC3 + pistar.ptr + ", rstar.nms[i-2]))
                dgp <- getBootDGP(c("PC1", "PC2", "PC3"), c("pistar.ptr", rstar.nms[i-2]), data[ind,], BC=TRUE)
                rval <- bootstrapTest(fmla.0, fmla.a, data[ind,], dgp, h, M, vcovfn)
                tbl[c(9,12), col] <- sprintf("[%4.2f]", rval$tblCoef[5, 4:5])
            }
            tbl[13:15, col] <- ""
        } else {
            ## i-star
            fmla <- xr ~ PC1sc + PC2sc + PC3sc + istar
            regrows <- c(1,3,5,13)
            tbl[nrow(tbl), col] <- "real-time"
            if (bootstrap_flag) {
                fmla.a <- xr ~ PC1 + PC2 + PC3 + istar
                dgp <- getBootDGP(c("PC1", "PC2", "PC3"), "istar", data[ind,], BC=TRUE)
                rval <- bootstrapTest(fmla.0, fmla.a, data[ind,], dgp, h, M, vcovfn)
                tbl[15, col] <- sprintf("[%4.2f]", rval$tblCoef[5, 4])
            }
            tbl[7:12, col] <- ""
        }
        mod <- lm(fmla, data[ind,])
        tbl[regrows, col] <- sprintf("%4.2f", mod$coef[-1])
        tbl[regrows+1, col] <- sprintf("(%4.2f)", sqrt(diag(vcovfn(mod)))[-1])
        tbl[nrow(tbl)-1, col] <- sprintf("%4.2f", summary(mod)$r.squared)
        col <- col+1
    }
    print(tbl)

    filename <- ifelse(subsample, "tables/returns_subsample.tex", "tables/returns.tex")
    print(xtable(tbl, digi=2),
          include.rownames=FALSE, include.colnames=FALSE, only.contents=TRUE,
          sanitize.text.function=function(x){x}, hline.after=NULL,
          file = filename)
}

