## Plot and analyze detrended ten-year yield
## Figure 3

graphics.off()
library(dynlm)
source("R/data_fns.R")
source("R/util_fns.R")

data <- loadData()
data$rstar <- data$rstar.realtime

demean <- function(x) x - mean(x)
p <- 4 # lags for DOLS
data$d.pistar.ptr <- c(NA, diff(data$pistar.ptr))
data$d.rstar <- c(NA, diff(data$rstar))
data$istar <- data$pistar.ptr + data$rstar
data$d.istar <- c(NA, diff(data$istar))

pdf("figures/cycles.pdf", width=7, height=5.5, pointsize=10)
cols <- c("black", "red", "steelblue")
par(mar=c(3,4,1,1), mgp=c(2,.6,0))
yrange <- range(demean(data$y10))
plot(data$date, demean(data$y10), type="l", col=cols[1], lwd=2, ylim=yrange, ylab="Percent", xlab="", xaxs = "i", yaxs = "i")
plot_recessions(range(data$date), yrange)
lines(data$date, demean(data$y10), col=cols[1], lwd=2)

## (1) simple difference
data$resid1 <- data$y10 - data$pistar.ptr - data$rstar
lines(data$date, data$resid1, col=cols[2], lwd=2)

## (2) cointegration residual using i*
fmla <- y10 ~ istar + L(d.istar, -p:p)
mod <- dynlm(fmla, data=ts(data))
data$resid2 <- data$y10 - drop(cbind(1, data$istar) %*% mod$coef[1:2])
lines(data$date, data$resid2, col=cols[3], lwd=2)
abline(h=0, lty=2)
legend("topright", c("Ten-year yield, demeaned",
                     expression(paste("Difference between yield and i*")),
                     expression(paste("Cointegration residual between yield and i*"))),
       col=cols, lwd=2, bg="white")
dev.off()

## numerical change in these series
tmp <- data[c("yyyymm", "y10", "resid1", "resid2")]
t1 <- which(data$yyyymm==199903)
t2 <- which(data$yyyymm==201803)
print(round(rbind(tmp[t1,], tmp[t2,], tmp[t2,]-tmp[t1,]), 2))
ind1 <- which(data$yyyymm%/%100 == 1999)
ind2 <- (nrow(data)-3):nrow(data)
tmp2 <- rbind(colMeans(tmp[ind1,]),
              colMeans(tmp[ind2,]))
tmp2 <- rbind(tmp2, tmp2[2,]-tmp2[1,])
print(round(tmp2, 2))
