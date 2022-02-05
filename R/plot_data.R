## Plot the data - ten-year yield, pi*, r*
## Figures 1 and 2

graphics.off()
source("R/data_fns.R")
source("R/util_fns.R")

data <- loadData()
cat("Sample period: ", range(data$yyyymm), "\n")

data$istar <- data$pistar.ptr + data$rstar.mean

## Figure 1: 10-year yield, pi*, and r*
pdf("figures/data1.pdf", width=7, height=5.5, pointsize=10)
cols <- c("black", "red", "green3", "steelblue", "darkgreen")
par(mar=c(3,4,1,1), mgp=c(2,.6,0))
yrange <- range(c(0, data$y10))
plot(data$date, data$y10, type="l", col=cols[1], lwd=2, ylim=yrange, xlab="", ylab="Percent", yaxp=c(2,14,6), xaxs = "i", yaxs = "i")
plot_recessions(range(data$date), yrange)
lines(data$date, data$y10, lwd=2, col=cols[1])
lines(data$date, data$pistar.ptr, col=cols[2], lwd=2)
lines(data$date, data$rstar.mean, col=cols[3], lwd=2)
lines(data$date, data$istar, col=cols[4], lwd=2)
legend("topright", c("Ten-year yield",
                     expression(paste('Trend Inflation, ',pi, '*')),
                     "Equilibrium real short rate, r*",
                     "Equilibrium short rate, i*"), col=cols, lwd=2, bg="white")
box()
dev.off()

## Figure 2: alternative measures of r*
rstars.sm <- as.matrix(data[attr(data, "nms.sm")])
nms.all <- attr(data, "nms.all")
rstars.all <- as.matrix(data[nms.all])
cat("Smoothed r-star estimates:", attr(data, "nms.sm"), "\n")
cat("Filtered r-star estimates:", attr(data, "nms.filt"), "\n")
cat("Real-time r-star estimates:", attr(data, "nms.realtime"), "\n")
cols <- palette()
pdf("figures/data2.pdf", width=7, height=5, pointsize=10)
par(mfrow=c(1,2),
    mar = c(2,3,2,0.5),
    mgp=c(2,.6,0))
yrange <- range(rstars.sm, rstars.all, na.rm=TRUE)
## left panel: smoothed estimates
plot(data$date, data$rstar.sm, type="l", col=cols[1], ylim=yrange, xlab="", ylab="Percent", main="", xaxs = "i", yaxs = "i")
range.min <- apply(rstars.sm, 1, min)
range.max <- apply(rstars.sm, 1, max)
polygon(x = c(data$date, rev(data$date)), y = c(range.min, rev(range.max)), density=NA, col="lightblue")
lines(data$date, data$rstar.sm, col=cols[1], lwd=2)
legend("topright", c("Avg. smoothed r*"), col=cols, lwd=2)
title("Smoothed estimates", font.main=1, cex.main=1.1)
abline(h=0, lty=2)
box()
## right panel: filtered/real-time
par(mar = c(2,2,2,1.5))
plot(data$date, data$rstar.filt, type="l", col=cols[1], ylim=yrange, xlab="", ylab="", main="", xaxs = "i", yaxs = "i")
range.min <- apply(rstars.all, 1, min)
range.max <- apply(rstars.all, 1, max)
polygon(x = c(data$date, rev(data$date)), y = c(range.min, rev(range.max)), density=NA, col="lightblue")
lines(data$date, data$rstar.filt, col=cols[1], lwd=2)
lines(data$date, data$rstar.realtime, col=cols[2], lwd=2)
lines(data$date, data$rr.ewma, col=cols[3], lwd=2)
legend("topright", c("Avg. filtered r*", "Avg. real-time r*", "Moving-average r*"), col=cols, lwd=2)
title("Filtered and real-time estimates", font.main=1, cex.main=1.1)
abline(h=0, lty=2)
box()
dev.off()

## report decrease in r-star
cat("Change since 2000:\n")
i1 <- which(data$yyyymm==200003)
i2 <- nrow(data)
X <- data[c("yyyymm", attr(data, "nms.all"), "rstar.mean", "y10")]
print(round(rbind(X[i1,], X[i2,], X[i2,]-X[i1,]), 2))

