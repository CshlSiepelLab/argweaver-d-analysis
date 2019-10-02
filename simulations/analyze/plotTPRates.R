require("argweaver")
require("rphast")

setwd("~/arghmm/SGDP/argweaver-d-analysis/simulations/generate/deepSims")
source("../../../scripts/functions.R")

doPDF <- FALSE


## want to make main roc plot for paper
## want to include new results with 50kya Sup->Afr event

cutoffs <- c(seq(from=0, to=1, length.out=11), 1.01)
cutoffs <- 0.5

runs1 <- c("mig50_div1500", "mig50_div1000", "mig150_div1500", "mig150_div1000", "mig250_div1500", "mig250_div1000", "mig350_div1500", "mig350_div1000")
runs2 <- runs1

usedirs <- list()
for (run1 in runs1) {
    for (run2 in runs2) {
        usedirs[[sprintf("%s.%s", run1, run2)]] <- list(run1=run1, subdir=sprintf("%s.2afr.5e-9_recomb", run2), fpdir="nomig")
    }
}

rocStats <- list()
migs <- c("hToN", "sToD", "sToA")
fpmigs <- c("hToD", "hToN", "sToD", "sToA", "sToN")
for (currname in names(usedirs)) {
    if (is.element(currname, names(rocStats))) next
    cat("usedir ", currname, " run1=", usedirs[[currname]]$run1, " subdir=", usedirs[[currname]]$subdir, "\n")
    rocStats[[currname]]$tprate <- matrix(nrow=3, ncol=length(cutoffs), dimnames=list(migs, 1:length(cutoffs)))
    rocStats[[currname]]$fprate <- matrix(nrow=5, ncol=length(cutoffs), dimnames=list(fpmigs, 1:length(cutoffs)))
    for (j in 1:length(cutoffs)) {
        cutoff <- cutoffs[j]
        cat("cutoff = ", cutoff, "\n")
        rv1 <- readSimResults(usedirs[[currname]]$run1,
                              subdir=usedirs[[currname]]$subdir, cutoff=cutoff)
        if (is.null(rv1)) {
            rocStats[[currname]] <- NULL
            break
        }
        tmp1 <- tapply(rv1$all$tprate, rv1$all$mig, mean)
        tmp1 <- tmp1[is.element(names(tmp1), migs)]
        rocStats[[currname]]$tprate[names(tmp1),j] <- tmp1
        rocStats[[currname]]$confusionMat <- rv1$confusionMat
        if (usedirs[[currname]]$fpdir != usedirs[[currname]]$run1)
            rv1 <- readSimResults(usedirs[[currname]]$fpdir,
                                  subdir=usedirs[[currname]]$subdir, cutoff=cutoff)
        if (!is.null(rv1)) {
            tmp1 <- tapply(rv1$all$fprate, rv1$all$mig, mean)
            rocStats[[currname]]$fprate[names(tmp1),j] <- tmp1
        } ## otherwise they are still NA
    }
}


migtimecol <- list("50"="orange","150"="red", "250"="blue","350"="purple","0"="black")

midpt <- which(cutoffs==0.5)
tps <- matrix(nrow=3, ncol=length(rocStats))
migs <-  rownames(rocStats[[1]]$tprate)
rownames(tps) <- migs
fps <- tps
migtimes1 <- c()
divtimes1 <- c()
migtimes2 <- c()
divtimes2 <- c()
r1 <- c()
r2 <- c()
cols <- c()
for (i in 1:length(rocStats)) {
    r1[i] <- strsplit(names(rocStats)[i], ".", fixed=TRUE)[[1]][1]
    r2[i] <- strsplit(names(rocStats)[i], ".", fixed=TRUE)[[1]][2]
    migtimes1[i] <- as.numeric(sub("mig", "", strsplit(r1[i], "_")[[1]][1]))
    divtimes1[i] <- as.numeric(sub("div", "", strsplit(r1[i], "_")[[1]][2]))
    migtimes2[i] <- as.numeric(sub("mig", "", strsplit(r2[i], "_")[[1]][1]))
    divtimes2[i] <- as.numeric(sub("div", "", strsplit(r2[i], "_")[[1]][2]))
    cols[i] <- migtimecol[[as.character(migtimes1[i])]]
    if (divtimes1[i] == 1e6) cols[i] <- makeTransparent(cols[i], 0.5)
    for (mig in migs) {
        tps[mig,i] <- rocStats[[i]]$tprate[mig,midpt]
        fps[mig,i] <- rocStats[[i]]$fprate[mig,midpt]
        if ((migtimes1[i] == 50 || migtimes2[i] == 50) && mig != "sToA") {
            tps[mig,i] <- NA
            fps[mig,i] <- NA
        }
    }
}
o <- order(migtimes1, divtimes1, migtimes2, divtimes2)
o <- order(migtimes1, divtimes1, migtimes2, divtimes2)
tps <- tps[,o]
fps <- fps[,o]
migtimes1 <- migtimes1[o]
migtimes2 <- migtimes2[o]
divtimes1 <- divtimes1[o]
divtimes2 <- divtimes2[o]
cols <- cols[o]

## What about 4x2 matrix
doDetailed <- TRUE
doFPrate <- FALSE

## When doFPrate is FALSE this produces figure 5
if (doPDF) {
if (doDetailed) {
    pdf(sprintf("sim%sratesFull.pdf", if (doFPrate) {"FP"} else {"TP"}),
        width=4, height=4)
} else {
    pdf(sprintf("sim%sratesSummary.pdf", if (doFPrate) {"FP"} else {"TP"}), width=4, height=4)
}
}


if (names(dev.cur()) != "pdf") {
    if (names(dev.cur()) != "null device") dev.off()
    x11(width=4, height=4)
}
migtimecol <- list("50"="orange","150"="red", "250"="purple","350"="blue")
currMigNames <- list(hToN="Hum to\nNea",
                     sToD="Sup to\nDen",
                     sToA="Sup to\nAfr")
migtimes <- c(50,150,250,350)
divtimes <- c(1000, 1500)
par(mar=c(4,4,2,1), mgp=c(2.5,1,0))
par(mfrow=c(1,1))
f <- migtimes1==150 & divtimes1==1500
tmp <- as.numeric(t(tps[,f]))
f2 <- !is.na(tmp)
if (!doDetailed) f2 <- f2 & divtimes2[f]==1000
layout(t(matrix(1:8, nrow=2, ncol=4)), widths=c(0.56,0.46), heights=c(0.28,0.23,0.23,0.31))
for (migtime1 in migtimes) {
    if (migtime1 == migtimes[1]) topmargin <- 1.5 else topmargin <- 0
    if (migtime1==tail(migtimes,1)) {
        labels <- as.character(currMigNames[rownames(tps)])
        bottommargin <- 3.5
    } else {
        labels <- rep("", nrow(tps))
        bottommargin <- 1
    }
    for (divtime1 in divtimes) {
        if (divtime1 == divtimes[1]) {
            leftmargin <- 5
            ylab <- ""
            rightmargin <- 0
        } else {
            leftmargin <- 0
            ylab <- ""
            rightmargin <- 1
        }
        par(mar=c(bottommargin, leftmargin,topmargin,rightmargin))
        f <- migtimes1==migtime1 & divtimes1==divtime1
        if (doFPrate) {
            tmp <- as.numeric(t(fps[,f]))
            ylim <- c(0, 0.01)
        } else {
            tmp <- as.numeric(t(tps[,f]))
            ylim <- c(0,1)
        }
        migcat <- as.character(sapply(rownames(tps), rep, ncol(tps[,f])))[f2]
        col <- as.character(migtimecol[rep(as.character(migtimes2[f]), nrow(tps))])[f2]
        spacing <- rep(0, length(migcat))
        spacing[1] <- 1.5
        for (i in 2:length(migcat)) if (migcat[i] != migcat[i-1]) spacing[i] <- 1
        bp <- barplot(tmp[f2], beside=TRUE, ylim=ylim, axes=FALSE, ylab=ylab,
                names.arg=rep("", length(tmp[f2])), col=col, space=spacing, las=2)
        if (divtime1 == divtimes[1]) {
            if (doFPrate) {
                axis(side=2, at=c(0, 0.005, 0.01), las=2)
            } else {
                axis(side=2, at=c(0, 0.5, 1.0), las=2)
            }
        }
        if (doDetailed) {
            barplot(tmp[f2], beside=TRUE, add=TRUE, axes=FALSE, ylab="", names.arg=rep("", length(tmp[f2])), col="black", space=spacing, density=ifelse(divtimes2[f2]==1500000,30,0))
        }
        if (migtime1 == migtimes[1]) {
            mtext(bquote(t[div] ~ "=" ~ .(divtime1/1e3) ~ "Mya"), cex=0.9)
        }
        if (doFPrate) {
            abline(h=c(0,0.005,0.01), lty=3)
        } else {
            abline(h=c(0,0.25,0.5,0.75, 1.0),lty=3)
        }
        abline(h=0)
        par(xpd=TRUE)
        par(xpd=FALSE)
        if (migtime1==migtimes[1] & divtime1==divtimes[1]) {
            par(xpd=TRUE)
            mtext("True Positive Rate", side=2, line=3, at=-1.25)
            par(xpd=FALSE)
        }
        if (divtime1 == divtimes[1]) {
            par(xpd=TRUE)
            text(x=par("usr")[2], y=0.85, bquote(t[mig] ~ "="), cex=1.4, adj=c(1,0.5))
            par(xpd=FALSE)
        }
        if (divtime1 == divtimes[2]) {
            par(xpd=TRUE)
            text(x=par("usr")[1], y=0.86, sprintf("%s kya", migtime1), cex=1.4, adj=c(-0.1,0.5))
            par(xpd=FALSE)
        }
        col <- migtimecol[[as.character(migtime1)]]
        if (divtime1 == 1.5e3) col <- makeTransparent(col, 0.15) else col <- makeTransparent(col,0.08)
        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], border=NA, col=col)
        if (migtime1 == tail(migtimes,1)) {
            tmp <- tapply(bp[,1], migcat, mean)
            mtext(currMigNames[names(tmp)], side=1, line=1.5, at=tmp, cex=0.8)
        }
    }
}
if (names(dev.cur())=="pdf") dev.off()
