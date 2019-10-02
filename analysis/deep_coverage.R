require("Hmisc")
require("argweaver")
require("rphast")
source("../scripts/functions.R")

dir <- "mig250_div1000_2afr/regions"
chrs <- c("autosome", "X")

migCol <- list(
    hToD="green",
    hToN="red",
    sToA="black",
    sToN="orange",
    sToD="blue")

x <- summarizeRealResults(dir, chrs=chrs)
x <- x[x$ind != "any",]
doPDF <- FALSE

width <- 6
height <- 3.5
cex <- 1.0

## Figure 6
pdfout <- "ancientBarplot.pdf"

if (doPDF) pdf(pdfout, width=width, height=height)

if (names(dev.cur()) != "pdf") {
#    dev.off()
    x11(width=width, height=height)
}
par(mar=c(4,4,1,1),mgp=c(2.5,1,0))
bp <- makeBarPlotFromSummary(x, cex=cex)
if (names(dev.cur()) == "pdf") dev.off()
