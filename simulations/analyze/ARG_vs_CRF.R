require("argweaver")
require("rphast")
source("../../scripts/functions.R")

indir <- "../generate/recentSims"

## first read CRF output
crfPred <- data.frame()
for (i in 1:100) {
    snps <- read.table(sprintf("%s/%d/crf/snp", indir, i))
    midpts <- as.numeric(sprintf("%.0f", (snps[-1,4]+snps[-nrow(snps),4])/2))
    tmp <- read.table(sprintf("%s/%d/crf/out/anc", indir, i))
    inds <- as.character(read.table(sprintf("%s/%d/crf/eur.ind", indir, i))[,1])
    if (nrow(snps) != nrow(tmp)) stop("error")
    crfPred <- rbind(crfPred,
                     data.frame(chrom=snps[,2], chromStart=c(0,midpts), chromEnd=c(midpts, 2e6),pos=snps[,4], 
                        eur1=tmp[,which(inds=="eur1")],
                        eur2=tmp[,which(inds=="eur2")],
                        eur3=tmp[,which(inds=="eur3")],
                        eur4=tmp[,which(inds=="eur4")]))
}

subdir <- "eurAfr_simple_model/5e-9_recomb_constMut"
argPred <- data.frame()
for (i in 1:100) {
    argPred <- rbind(argPred, readArgSummary(sprintf("%s/%d/%s/out.migStatsHap.bed.gz", indir, i, subdir)))
}
names(argPred)[4:7] <- sprintf("eur%i", 1:4)

## read true regions
trueRegions <- list()
for (ind in sprintf("eur%i", 1:4)) {
    trueRegions[[ind]] <- data.frame()
    for (i in 1:100) {
        tmp <- NULL
        try(tmp <- read.table(sprintf("%s/%d/%s.nToX.txt", indir, i, ind)), silent=TRUE)
        if (!is.null(tmp)) trueRegions[[ind]] <- rbind(trueRegions[[ind]], tmp)
    }
    names(trueRegions[[ind]]) <- c("chrom", "chromStart", "chromEnd")
}
for (ind in names(trueRegions)) trueRegions[[ind]] <- flatten.feat(bedToFeat(trueRegions[[ind]]))

allRegions <- flatten.feat(bedToFeat(argPred))
cutoffs <- seq(from=0, to=1.01, length.out=100)
argPerf <- data.frame()
for (maxlen in c(NA, 50000, 25000, 10000)) {
    useAll <- allRegions
    if (is.na(maxlen)) {
        useTrue <- trueRegions
    } else {
        useTrue <- list()
        for (ind in names(trueRegions)) {
            f <- (trueRegions[[ind]]$end - trueRegions[[ind]]$start <= maxlen)
            useTrue[[ind]] <- trueRegions[[ind]][f,]
            useAll <- coverage.feat(allRegions, useAll, trueRegions[[ind]][!f,], get.feats=TRUE, not=c(FALSE, FALSE, TRUE))
        }
    }
for (cutoff in cutoffs) {
    tp <- 0
    fp <- 0
    tn <- 0
    fn <- 0
    for (ind in names(trueRegions)) {
        predRegions <- bedToFeat(argPred[argPred[,ind] >= cutoff,])
        tp <- tp + coverage.feat(useAll, predRegions, useTrue[[ind]])
        fp <- fp + coverage.feat(allRegions, useAll, predRegions, useTrue[[ind]],
                                 not=c(FALSE, FALSE, FALSE, TRUE))
        tn <- tn + coverage.feat(allRegions, useAll, useTrue[[ind]], predRegions,
                                 not=c(FALSE, FALSE, TRUE, TRUE))
        fn <- fn + coverage.feat(allRegions, useAll, useTrue[[ind]], predRegions,
                                 not=c(FALSE, FALSE, FALSE, TRUE))
    }
    argPerf <- rbind(argPerf, data.frame(maxlen=maxlen, maxcutoff=cutoff, tp=tp, fp=fp, tn=tn, fn=fn, tprate=tp/(tp+fn), fprate=fp/(fp+tn)))
}
}



crfFeat <- feat(crfPred$chrom, start=crfPred$pos, end=crfPred$pos)
crfSnps <- sprintf("%s:%s", crfFeat$seqname, crfFeat$end)
crfPerf <- data.frame()
for (maxlen in c(NA, 50000, 25000, 10000)) {
trueSnps <- list()
negSnps <- list()
tooLongSnps <- c()
for (ind in names(trueRegions)) {
    cat(maxlen, ind,"\n")
    if (is.na(maxlen)) {
        trueSnps[[ind]] <- coverage.feat(trueRegions[[ind]], crfFeat, get.feats=TRUE)
    } else {
        f <- trueRegions[[ind]]$end - trueRegions[[ind]]$start + 1 <= maxlen
        trueSnps[[ind]] <- coverage.feat(trueRegions[[ind]][f,], crfFeat, get.feats=TRUE)
        tooLong <- coverage.feat(trueRegions[[ind]][!f,], crfFeat, get.feats=TRUE)
        tooLongSnps <- unique(c(tooLongSnps,sprintf("%s:%s", tooLong$seqname, tooLong$end)))
    }
    trueSnps[[ind]] <- sprintf("%s:%s", trueSnps[[ind]]$seqname, trueSnps[[ind]]$end)
    negSnps[[ind]] <- crfSnps[!is.element(crfSnps, trueSnps[[ind]])]
}
if (!is.na(maxlen)) {
    for (ind in names(negSnps)) {
        cat("Removing long elements from negSnps: ", ind, "\n")
        negSnps[[ind]] <- negSnps[[ind]][!is.element(negSnps[[ind]], tooLongSnps)]
    }
}
for (cutoff in cutoffs) {
    cat(maxlen, cutoff,"\n")
    tp <- 0
    fp <- 0
    tn <- 0
    fn <- 0
    for (ind in names(trueRegions)) {
        f <- crfPred[,ind] >= cutoff
        predSnps <- sprintf("%s:%s", crfPred[f,"chrom"], crfPred[f,"pos"])
        notPred <- crfSnps[!is.element(crfSnps, predSnps)]
        tp <- tp + sum(is.element(predSnps, trueSnps[[ind]]))
        tn <- tn + sum(is.element(notPred, negSnps[[ind]]))
        fp <- fp + sum(is.element(predSnps, negSnps[[ind]]))
        fn <- fn + sum(is.element(notPred, trueSnps[[ind]]))
    }
    crfPerf <- rbind(crfPerf, data.frame(maxlen=maxlen, cutoff=cutoff, tp=tp, fp=fp, tn=tn, fn=fn, tprate=tp/(tp+fn), fprate=fp/(fp+tn)))
    print(crfPerf[nrow(crfPerf),])
}
}
    

## figure 2
pdf("ARG_vs_CRF_ROC.pdf", width=8, height=4)

if (names(dev.cur()) != "pdf") {
    if (names(dev.cur()) != "null device") dev.off()
    x11(width=8, height=4)
}
par(mfrow=c(1,2), mar=c(4,4,2,1), mgp=c(2.25,1,0))
plot(argPerf$fprate, argPerf$tprate, xlim=c(0, 0.1), col="red", type="n", xlab="FP Rate", ylab="TP Rate")
lty <- 1
for (maxlen in c(NA, 50e3, 25e3, 10e3)) {
    if (is.na(maxlen)) f <- is.na(argPerf$maxlen) else f <- !is.na(argPerf$maxlen) & argPerf$maxlen == maxlen
    if (sum(f) > 0) 
        lines(argPerf[f,"fprate"], argPerf[f,"tprate"], col="red", lty=lty, lwd=2)
    if (is.na(maxlen)) {
        f <- is.na(crfPerf$maxlen)
    } else f <- !is.na(crfPerf$maxlen) & crfPerf$maxlen == maxlen
    if (sum(f) > 0) 
        lines(crfPerf[f,"fprate"], crfPerf[f,"tprate"], lty=lty, lwd=2, col="blue")
    lty <- lty+1
}
#tmp <- legend(x="bottomright", legend=sapply(c("", ""), rep, 4), col=sapply(c("red", "blue"), rep, 4), lty=c(1,2,3,4), lwd=2, ncol=2, bty="n")
#text(c("all", "<= 50kb", "<= 25kb", "<=10kb"), x=tmp$rect$left, y=unique(tmp$text$y), pos=2)
                                        #text(c("ARG", "CRF"), x=seq(from=tmp$rect$left, to=tmp$rect$left + tmp$rect$w, length.out=5)[c(2,4)], y=tmp$rect$top)
#rect(xleft=tmp$rect$left-tmp$rect$w*0.575, xright=par("usr")[2], ybottom=par("usr")[3], ytop=0.4)
legend(x="bottomright", c("all", "<= 50kb", "<= 25kb", "<= 10kb"), lwd=2, lty=c(1,2,3,4))
mtext("A", side=2, line=2, at=1.1, las=2, cex=1.3)

pred0 <- list(arg=argPred, crf=crfPred)
pred0$arg$score <- apply(pred0$arg[,c("eur1","eur2","eur3","eur4")],1,max)
pred0$crf$score <- apply(pred0$crf[,c("eur1","eur2","eur3","eur4")],1,max)

## Look at lengths of predictions.
cutoff <- 0.5
trueLens <- unlist(sapply(trueRegions, function(x) {lengths.feat(x)}))
myPlotCdf(trueLens, col="black", lwd=2, xlab="Length", ylab="CDF", xlim=c(0, 6e5))
for (i in 1:length(pred0)) {
    f <- pred0[[i]]$score >= 0.5
    if (names(pred0)[i] == "crf") col <- "blue" else col <- "red"
    myPlotCdf(lengths.feat(flatten.feat(bedToFeat(pred0[[i]][f,]))), col=col, add=TRUE, lwd=2)
}
legend(x="bottomright", c("True", "ARG", "CRF"), lwd=2, col=c("black", "red", "blue"))
mtext("B", side=2, line=2, at=1.1, las=2, cex=1.3)

if (names(dev.cur()) == "pdf") {
dev.off()
}
