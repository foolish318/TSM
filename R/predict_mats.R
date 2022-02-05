## Predict excess returns for individual bond maturities
## Online Appendix Table 4

library(xtable) # xtable()
library(VAR.etp) # for bootstrap
source("R/data_fns.R")
source("R/util_fns.R")
source("R/bootstrap_fns.R")
M <- 5000 # number of bootstrap replications

data <- loadData()
data$istar <- data$pistar.ptr + data$rstar.realtime
yield.cols <- attr(data, "yield.cols")
mats <- as.numeric(substr(yield.cols, 2, 5))
Y <- data.matrix(data[,yield.cols])

##################################################
filename <- "tables/returns_mats.tex"

### quarterly/annual returns
## h <- 4
## vcovfn <- function(mod) sandwich::NeweyWest(mod, lag=6, prewhite=FALSE)
h <- 1
vcovfn <- function(mod) sandwich::vcovHC(mod, type="HC0")

### sample period
## cat("Post-1985 subsample:\n")
## ind <- data$yyyymm >= 198501
cat("Full sample:\n")
ind <- 1:nrow(data)
##################################################

pred.mats <- c(2, 5, 7, 10, 15)

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

xrn <- excess_returns(Y, mats, h)
data$xr <- rowMeans(xrn)
colnames(xrn) <- paste0("xr", 2:15)
data[paste0("xr", pred.mats)] <- xrn[,paste0("xr", pred.mats)]

## table for results
tbl <- data.frame(matrix(NA, 2*length(pred.mats), 10))
tbl[,1] <- c(t(cbind(paste0(pred.mats, "y"), "")))

dgp <- getBootDGP(c("PC1", "PC2", "PC3"), "istar", data[ind,], BC=TRUE)

for (mat in pred.mats) {
    xrnm <- paste0("xr", mat)
    row <- which(tbl[,1] == paste0(mat, "y"))
    fmla.0 <- formula(paste(xrnm, "~ PC1 + PC2 + PC3"))

    ## yields only
    fmla <- formula(paste(xrnm, "~ PC1sc + PC2sc + PC3sc"))
    mod <- lm(fmla, data[ind,])
    tbl[row, 2:4] <- sprintf("%4.2f", mod$coef[-1])
    tbl[row+1, 2:4] <- sprintf("(%4.2f)", sqrt(diag(vcovfn(mod)))[-1])
    tbl[row, 5] <- sprintf("%4.2f", summary(mod)$r.squared)
    tbl[row+1, 5] <- ""

    ## i-star
    fmla <- formula(paste(xrnm, "~ PC1sc + PC2sc + PC3sc + istar"))
    mod <- lm(fmla, data[ind,])
    tbl[row, 6:9] <- sprintf("%4.2f", mod$coef[-1])
    tbl[row+1, 6:9] <- sprintf("(%4.2f)", sqrt(diag(vcovfn(mod)))[-1])
    tbl[row, 10] <- sprintf("%4.2f", summary(mod)$r.squared)
    ## bootstrap
    fmla.a <- formula(paste(xrnm, "~ PC1 + PC2 + PC3 + istar"))
    rval <- bootstrapTest(fmla.0, fmla.a, data[ind,], dgp, h, M, vcovfn)
    tbl[row+1, 10] <- sprintf("[%4.2f]", rval$tblCoef[5, 4])
}
print(tbl)

print(xtable(tbl, digi=2),
      include.rownames=FALSE, include.colnames=FALSE, only.contents=TRUE,
      sanitize.text.function=function(x){x}, hline.after=NULL,
      file = filename)


