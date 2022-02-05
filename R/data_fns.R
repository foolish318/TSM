## Functions to load and process the data

make_quarterly <- function(data) {
    month <- as.numeric(format(data$date, "%m"))
    stopifnot(!any(is.na(month)))
    data[month %% 3 == 0,]
}

loadLW <- function() {
    ## Laubach-Williams
    lw <- read.csv("data/rstar/lw.csv" , na.strings = "#N/A")
    lw$Date <- as.Date(lw$Date, format="%m/%d/%Y")
    lw$yyyymm <- as.numeric(format(lw$Date, "%Y%m"))
    names(lw)[2:3] <- c("rstar.lw", "rstar.lw.sm")
    lw[c("yyyymm", "rstar.lw", "rstar.lw.sm")]
}

loadHLW <- function() {
    hlw <- read.csv("data/rstar/hlw.csv")
    hlw$Date <- as.Date(hlw$Date, "%m/%d/%Y")
    hlw$yyyymm <- as.numeric(format(hlw$Date, "%Y%m"))+2
    names(hlw)[2] <- "rstar.hlw"
    hlw[c("yyyymm", "rstar.hlw")]
}

Matlab2yyyymm <- function(date) {
    ## transform Matlab's dates (Kiley/JM), 2012.00 = 2012:Q1, 2012.75 = 2012:Q4
    year <- floor(date)
    x <- date - year
    quarter <- match(x, x[1:4])
    year*100 + quarter*3
}

loadKiley <- function() {
    ## Kiley
    kiley <- read.csv("data/rstar/kiley.csv", header=FALSE)
    names(kiley) <- c("date", "rstar.kiley.sm", "rstar.kiley")
    kiley$yyyymm <- Matlab2yyyymm(kiley$date)
    kiley <- kiley[c("yyyymm", "rstar.kiley", "rstar.kiley.sm")]
    kiley
}

loadJM <- function() {
    jm <- read.csv("data/rstar/jm_realtime.csv")
    jm$yyyymm <- Matlab2yyyymm(jm$Date)
    names(jm)[1] <- "rstar.jm"
    jm <- jm[c("yyyymm", "rstar.jm")]
    jm2 <- read.csv("data/rstar/jm_smoothed.csv")
    jm2$yyyymm <- Matlab2yyyymm(jm2$Date)
    names(jm2)[1] <- "rstar.jm.sm"
    jm <- merge(jm, jm2[c("yyyymm", "rstar.jm.sm")])
}

loadPTR <- function() {
    ptr <- read.csv("data/pistar_PTR.csv")
    ptr$yyyymm <- ptr$Year*100+(ptr$Quarter)*3
    ptr <- ptr[c("yyyymm", "pistar_PTR")]
    names(ptr)[2] <- "pistar.ptr"
    ptr
}

loadPCEPI <- function(core=TRUE, yoy=TRUE) {
    ## PCE price index - seasonally adjusted
    if (core) {
        pce <- read.csv("data/PCEPILFE.csv")
    } else {
        pce <- read.csv("data/PCEPI.csv")
    }
    names(pce) <- c("date", "pce")
    pce$date <- as.Date(pce$date)
    pce <- make_quarterly(pce)
    pce$yyyymm <- as.numeric(format(pce$date, "%Y%m"))
    if (yoy) {
        pce$pi <- 100*c(rep(NA, 4), diff(log(pce$pce), 4)) # year-over-year inflation
    } else {
        pce$pi <- 400*c(NA, diff(log(pce$pce))) # annualized quarterly inflation
    }
    pce <- na.omit(pce)
    pce[c("date", "yyyymm", "pi")]
}

loadTB <- function() {
    ## three-month T-bill rate - starts in 1934-01
    tb <- read.csv("data/TB3MS.csv") # quarterly
    names(tb) <- c("date", "tb3m")
    tb$date <- as.Date(tb$date)
    tb <- make_quarterly(tb)
    tb$yyyymm <- as.numeric(format(tb$date, "%Y%m"))
    tb <- tb[tb$yyyymm >= 195201,] # the year after the Treasury accord
    tb
}

loadEPRR <- function() {
    ## quarterly ex-post real rate = nominal rate - realized year-over-year inflation
    ##  nominal: 3m T-bill rate
    ##  inflation: Core PCE inflation
    inf <- loadPCEPI(core=TRUE, yoy=TRUE)
    tb <- loadTB()
    rr <- merge(inf, tb[c("yyyymm", "tb3m")], all.x = TRUE)
    rr$rr <- rr$tb3m - rr$pi
    rr <- na.omit(rr)
    rr
}

ewma <- function(x, v = 0.95, y0) {
## ewma <- function(x, v = 0.98, y0) {
    ## calculate exponentially-weighted moving average as in Cieslak-Povala
    ## v = 0.95 is the default smoothing parameter, appropriate for quarterly data
    tau <- numeric(length(x))
    if (missing(y0)) {
        tau[1] <- x[1]
    } else {
        tau[1] <- y0
    }
    for (t in 2:length(x))
        tau[t] <- tau[t-1] + (1-v)*(x[t] - tau[t-1])
    tau
}

loadData <- function(start = 197112, end = 201803) {
    ## (quarterly) data set for analysis in the paper

    ## yields
    data <- read.csv("data/yields.csv") # end-of-quarter yields
    yield.cols <- names(data)[-1]
    data$yyyymm <- floor(data$Date/100)
    data$year <- floor(data$yyyymm/100)
    data$month <- data$yyyymm-data$year*100
    data$date <- as.Date(as.character(data$Date), format="%Y%m%d")
    data$Date <- NULL
    data <- subset(data, yyyymm <= end & yyyymm >= start) # right away select subsample
    Y <- data.matrix(data[,yield.cols])
    w <- eigen(cov(Y))$vectors[,1]
    data$level <- Y %*% w / sum(w)

    ## pi-star estimate: PTR
    data <- merge(data, loadPTR(), all.x = TRUE)

    ## r-star estimates
    data <- merge(data, loadLW(), all.x = TRUE)
    data <- merge(data, loadHLW(), all.x = TRUE)
    data <- merge(data, loadKiley(), all.x = TRUE)
    data <- merge(data, loadJM(), all.x = TRUE)

    ## ex-post real rates and moving average
    rr <- loadEPRR()
    rr <- subset(rr, yyyymm >= 196111)
    rr$rr.ewma <- ewma(rr$rr, v=0.98)
    data <- merge(data, rr[c("yyyymm", "rr", "rr.ewma")], all.x=TRUE)

    ## Del Negro et al
    dn <- read.csv("data/rstar/delnegro_realtime.csv")
    names(dn) <- c("yyyymm", "rstar.dn", "X", "Y")
    data <- merge(data, dn[c("yyyymm", "rstar.dn")], all.x = TRUE)
    dn <- read.csv("data/rstar/delnegro.csv")
    names(dn) <- c("yyyymm", "rstar.dn.sm", "rstar.dn.lb", "rstar.dn.ub")
    data <- merge(data, dn[c("yyyymm", "rstar.dn.sm")], all.x = TRUE)

    ## Proxies model
    proxies <- read.csv("data/rstar/proxies_realtime.csv")
    proxies$rstar.proxies <- proxies$rstar
    data <- merge(data, proxies[c("yyyymm", "rstar.proxies")], all.x = TRUE)
    proxies <- read.csv("data/rstar/proxies.csv")
    proxies$rstar.proxies.sm <- proxies$rstar
    data <- merge(data, proxies[c("yyyymm", "rstar.proxies.sm")], all.x = TRUE)

    ## UC model
    uc <- read.csv("data/rstar/uc_realtime.csv")
    uc$rstar.uc <- uc$rstar
    data <- merge(data, uc[c("yyyymm", "rstar.uc")], all.x = TRUE)
    uc <- read.csv("data/rstar/uc.csv")
    uc$rstar.uc.sm <- uc$rstar.mcmc
    data <- merge(data, uc[c("yyyymm", "rstar.uc.sm")], all.x = TRUE)

    ## State-space model
    ssm <- read.csv("data/rstar/ssm.csv")
    ssm$rstar.ssm.sm <- ssm$rstar
    ssm$pistar.ssm.sm <- ssm$pistar
    data <- merge(data, ssm[c("yyyymm", "rstar.ssm.sm", "pistar.ssm.sm")], all.x = TRUE)
    ssm <- read.csv("data/rstar/ssm_realtime.csv")
    ssm$rstar.ssm <- ssm$rstar
    ssm$pistar.ssm <- ssm$pistar
    data <- merge(data, ssm[c("yyyymm", "rstar.ssm", "pistar.ssm")], all.x = TRUE)

    ## make last observation of MCMC real-time estimates same as full-sample estimates
    ## (differences only due to random MCMC sampling)
    nobs <- nrow(data)
    nms <- c("rstar.dn", "rstar.uc", "rstar.proxies", "rstar.ssm", "pistar.ssm")
    data[nobs, nms] <- data[nobs, paste(nms, "sm", sep=".")]

    ## average of real-time estimates
    nms.realtime <- c("rstar.jm", "rstar.dn", "rstar.uc", "rstar.proxies", "rstar.ssm", "rr.ewma")
    data$rstar.realtime <- rowMeans(data.matrix(data[nms.realtime]), na.rm=TRUE)

    ## average of filtered estimates
    nms.filt <- c("rstar.lw", "rstar.kiley", "rstar.hlw")
    data$rstar.filt <- rowMeans(data.matrix(data[nms.filt]), na.rm=TRUE)

    ## average of all filtered and real-time estimates
    nms.all <- unique(c(nms.filt, nms.realtime))
    data$rstar.mean <- rowMeans(data.matrix(data[nms.all]), na.rm=TRUE)

    ## average of smoothed estimates
    nms.sm <- c("rstar.lw.sm", "rstar.jm.sm", "rstar.dn.sm", "rstar.kiley.sm", # existing estimates - HLW n/a
                "rstar.uc.sm", "rstar.proxies.sm", "rstar.ssm.sm")             # our estimates
    data$rstar.sm <- rowMeans(as.matrix(data[nms.sm]), na.rm=TRUE)

    attr(data, "yield.cols") <- yield.cols
    attr(data, "nms.realtime") <- nms.realtime
    attr(data, "nms.filt") <- nms.filt
    attr(data, "nms.sm") <- nms.sm
    attr(data, "nms.all") <- nms.all
    stopifnot(all(data$month %% 3 == 0)) ## make sure it's quarterly
    rownames(data) <- NULL
    data
}

loadBlueChip <- function() {
    ## Blue Chip long-range forecasts of ten-year Treasury yield
    ## Expected file format:
    ##  Column 1: Forecast date (format %m/%d/%Y)
    ##  Columns 2-6: Years of annual forecasts (integer)
    ##  Columns 7-11: Annual forecasts (numeric)
    ##  Column 12: Long-range forecast (not used, average over horizons of about 6-10 years)
    bc <- read.csv("data/bluechip_10y.csv", header=TRUE)
    bcyearcols <- paste0("year", 1:5)
    bcfcols <- paste0("f", 1:5)
    names(bc) <- c("Date", bcyearcols, bcfcols, "lr")
    bc$Date <- as.Date(bc$Date, "%m/%d/%Y")
    bc$year <- as.numeric(format(bc$Date, "%Y"))
    bc$month <- as.numeric(format(bc$Date, "%m"))
    bc$yyyymm <- bc$year*100+bc$month
    bc
}
