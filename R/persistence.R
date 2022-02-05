## Persistence of interest rates, detrended interest rates, trends
## Online Appendix Table 2

library(xtable)
library(urca)
source("R/data_fns.R")
source("R/ur_tests.R")

data <- loadData()
data$istar <- data$pistar.ptr + data$rstar.realtime

rstar.nms <- c("rstar.filt", "rstar.realtime", "rr.ewma")
rstar.tex <- c("r_t^{\\ast, F}", "r_t^{\\ast, RT}", "r_t^{\\ast, MA}")
yield.nms <- c("y2", "y10", "level")
yield.tex <- c("y^{(8)}_t", "y^{(40)}_t", "L_t")

nms <- c("Series", "SD", "ACF(1)", "Half-life", "ADF", "PP", "LFST")
tbl <- data.frame(matrix(NA, (length(yield.nms)+1)*(2+length(rstar.nms)), length(nms)))
colnames(tbl) <- nms
row <- 1
for (i in seq_along(yield.nms)) {
    ## yield
    tbl[row, 1] <- paste0("$", yield.tex[i], "$")
    tbl[row, -1] <- persistence_stats(data[[yield.nms[i]]])
    row <- row+1
    ## yield minus pistar
    tbl[row, 1] <- paste0("$", yield.tex[i], " - \\pi_t^\\ast$")
    tbl[row, -1] <- persistence_stats(data[[yield.nms[i]]] - data$pistar.ptr)
    row <- row+1
    ## yield minus pistar minus rstar
    for (j in seq_along(rstar.nms)) {
        tbl[row, 1] <- paste0("$", yield.tex[i], "- \\pi_t^\\ast - ", rstar.tex[j], "$")
        tbl[row, -1] <- persistence_stats(data[[yield.nms[i]]] - data$pistar.ptr - data[[rstar.nms[j]]])
        row <- row+1
    }
}
## macro trends
trend.nms <- c("pistar.ptr", rstar.nms, "istar")
trend.tex <- c("\\pi_t^\\ast", rstar.tex, "i^\\ast_t = \\pi_t^\\ast + r_t^{\\ast, RT}")
for (i in seq_along(trend.nms)) {
    tbl[row, 1] <- paste0("$", trend.tex[i], "$")
    tbl[row, -1] <- persistence_stats(data[[trend.nms[i]]])
    row <- row+1
}

tbl[,1] <- sprintf(paste0("%-", max(nchar(tbl[,1])),"s"), tbl[,1])
print(tbl)

sink("tables/persistence.tex")
print(xtable(tbl, digi=2),
      include.rownames=FALSE, include.colnames=FALSE, only.contents=TRUE,
      sanitize.text.function=function(x){x}, hline.after=c(5,10,15))
sink()

