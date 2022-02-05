## Plot various estimates of r-star
## Online Appendix Figures 1 and 2

graphics.off()
source("R/data_fns.R")
source("R/util_fns.R")

data <- loadData()
nms.sm <- attr(data, "nms.sm")
nms.filt <- attr(data, "nms.filt")
nms.realtime <- attr(data, "nms.realtime")
cat("Date range:", range(data$yyyymm), "\n")
cols <- palette()

## Online Appendix Figure 1
pdf("figures/app_rstar_external.pdf", width=7, height=4, pointsize=12)
par(mfrow=c(1,3), mar=c(2,2,2,1), mgp=c(2,.6,0))
nms1 <- c("rstar.dn.sm", "rstar.jm.sm", "rstar.lw.sm", "rstar.kiley.sm")
nms1d <- c("Del Negro et al. (2017)", "Johannsen and Mertens (2016)", "Laubach and Williams (2016)", "Kiley (2015)")
nms2 <- c("rstar.lw", "rstar.hlw", "rstar.kiley")
nms2d <- c("Laubach and Williams (2016)", "Holston et al. (2017)", "Kiley (2015)")
nms3 <- c("rstar.dn", "rstar.jm")
nms3d <- c("Del Negro et al. (2017)", "Johannsen and Mertens (2016)")
yrange <- range(as.matrix(data[c(nms1,nms2,nms3)]))
## left panel: smoothed
plot(data$date, data[[nms1[1]]], type="l", col=cols[1], lwd=2, ylim=yrange, xlab="", ylab="Percent", xaxs = "i", yaxs = "i")
title("Smoothed estimates", font.main=1, cex.main=1.1)
for (i in 2:length(nms1))
    lines(data$date, data[[nms1[i]]], col=cols[i], lwd=2)
legend("topright", nms1d, lwd=2, lty=1, col=cols, cex=0.9, bg="white")
## middle panel: filtered
plot(data$date, data[[nms2[1]]], type="l", col=cols[1], lwd=2, ylim=yrange, xlab="", ylab="Percent", xaxs = "i", yaxs = "i")
title("Filtered estimates", font.main=1, cex.main=1.1)
for (i in 2:length(nms2))
    lines(data$date, data[[nms2[i]]], col=cols[i], lwd=2)
legend("topright", nms2d, lwd=2, lty=1, col=cols, cex=0.9, bg="white")
## right panel: real-time
plot(data$date, data[[nms3[1]]], type="l", col=cols[1], lwd=2, ylim=yrange, xlab="", ylab="Percent", xaxs = "i", yaxs = "i")
title("Real-time estimates", font.main=1, cex.main=1.1)
for (i in 2:length(nms3))
    lines(data$date, data[[nms3[i]]], col=cols[i], lwd=2)
legend("topright", nms3d, lwd=2, lty=1, col=cols, cex=0.9)
dev.off()

## Online Appendix Figure 2
pdf("figures/app_rstar_internal.pdf", width=5.5, height=3.5, pointsize=9)
par(mfrow=c(1,2), mar=c(2,2,2,1), mgp=c(2,.6,0))
nms1 <- c("rstar.uc.sm", "rstar.proxies.sm", "rstar.ssm.sm")
nms.d <- c("UC model", "Proxies model", "SSM", "Moving average")
nms2 <- c("rstar.uc", "rstar.proxies", "rstar.ssm", "rr.ewma")
yrange <- range(as.matrix(data[c(nms1,nms2)]))
## left panel: smoothed
plot(data$date, data[[nms1[1]]], type="l", col=cols[1], lwd=2, ylim=yrange, xlab="", ylab="Percent", xaxs = "i", yaxs = "i")
title("Smoothed estimates", font.main=1, cex.main=1.1)
for (i in 2:length(nms1))
    lines(data$date, data[[nms1[i]]], col=cols[i], lwd=2)
## right panel: real-time
plot(data$date, data[[nms2[1]]], type="l", col=cols[1], lwd=2, ylim=yrange, xlab="", ylab="Percent", xaxs = "i", yaxs = "i")
for (i in 2:length(nms2))
    lines(data$date, data[[nms2[i]]], col=cols[i], lwd=2)
title("Real-time estimates", font.main=1, cex.main=1.1)
legend("bottomright", nms.d, lwd=2, lty=1, col=cols, cex=0.9, bg="white")
dev.off()

