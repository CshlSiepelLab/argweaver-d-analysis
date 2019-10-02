require("spatstat")

bedToFeat <- function(f) {
    feat(seqname=f[,1], start=f[,2]+1, end=f[,3])
}

lengths.feat <- function(x) {
    x$end - x$start + 1
}

readStarch <- function(f, ...) {
    if (!file.exists(f)) stop(f, " does not exist")
    rtry <- try(x <- read.table(pipe(sprintf("unstarch %s", f)), header=FALSE, stringsAsFactors=FALSE, ...))
    if (inherits(rtry, "try-error")) {
        stop("Error reading ", f)
    }
    names(x)[1:3] <- c("chrom", "chromStart", "chromEnd")
    x
}



## if x is given those points will be used on x axis; otherwise will space evenly across range of x
myPlotCdf <- function(dist, x=NULL, weights=rep(1/length(dist), length(dist)), xlim=NULL, npoints=100, add=FALSE, type="l", normalize=TRUE, ...) {
    if (is.null(x)) {
        if (is.null(xlim)) {
            if (add) {
                xlim <- par("usr")[1:2]
            } else {
                xlim <- range(dist, na.rm=TRUE, finite=TRUE)
            }
        }
        x <- seq(from=xlim[1], to=xlim[2], length.out=npoints)
    }
    y <- ewcdf(dist, weights/sum(weights))(x)
    if (!normalize) {
        y <- y * sum(dist)
    } 
    if (!add) {
        plot(x, y, type=type, xlim=xlim, ...)
    } else lines(x, y, type=type, ...)
    invisible(list(x=x, y=y))
}

readSimResults <- function(dir="mig250_div1000", subdir="mig250_div1000.2afr.5e-9_recomb",
                           cutoff=0.5, afrs=c("afr1", "afr2")) {
    ## first read true regions
    if (dir != "nomig" && dir != "nomig_X") {
    trueFiles <- data.frame(ind=c("alt", "vin", "den", afrs),
                            mig=c("hToN", "hToN", "sToD", rep("sToA", length(afrs))))
    trueFiles$file <- sprintf("%s/true_regions/%s.%s.bed", dir, trueFiles$mig, trueFiles$ind)
    trueRegions <- list()
    for (i in 1:nrow(trueFiles)) {
        file <- trueFiles[i,"file"]
        mig <- trueFiles[i,"mig"]
        ind <- trueFiles[i,"ind"]
        if (!file.exists(file)) {
            if (mig == "hToN") {
                file <- sprintf("%s/true_regions/%s.%s.bed", dir, "aToN", ind)
            }
            if (mig == "sToD") {
                file <- sprintf("%s/true_regions/%s.merged.bed", dir, mig)
            }
            if (!file.exists(file)) {
                next
                stop("Error finding ", file)
            }
        }
        tmpfeat <- NULL
        try(tmpfeat <- read.feat(file), silent=TRUE)
        if (!is.null(tmpfeat)) {
            trueRegions[[sprintf("%s.%s", mig, ind)]] <- flatten.feat(tmpfeat)
        }
    }
    } else trueRegions <- NULL
    ## now read pred regions; need hom/het for each ind
    allRegions <- NULL
    pred <- list()
    predfiles <- list.files(sprintf("%s/combined_pred/%s/", dir, subdir),
                            pattern="out\\..*.bed.starch", full.names=TRUE)
    if (length(predfiles) == 0) return(NULL)
    for (predfile in predfiles) {
        bn <- basename(predfile)
        mig <- strsplit(bn, ".", fixed=TRUE)[[1]][2]
        if (mig == "sToH") mig <- "sToA"
        if (mig == "aToN") mig <- "hToN"
        if (mig == "aToD") mig <- "hToD"
        type <- strsplit(bn, ".", fixed=TRUE)[[1]][4]
        if (type == "bed") next
        if (type != "hom" && type != "het") stop("Error got type ", type)
        ind <- strsplit(bn, ".", fixed=TRUE)[[1]][3]
        str <- sprintf("%s.%s.%s", mig, ind, type )
        tmppred <- readStarch(predfile)
        if (is.null(allRegions)) allRegions <- flatten.feat(bedToFeat(tmppred))
        pred[[str]] <- flatten.feat(bedToFeat(tmppred[tmppred[,4] >= cutoff,]))
    }

    ## reduce trueRegions to regions where we made predictions, in case all ARGweaver wasn't run on all replicates
    if (!is.null(trueRegions)) {
        for (trueIdx in names(trueRegions)) {
            trueRegions[[trueIdx]] <- coverage.feat(trueRegions[[trueIdx]], allRegions, get.feats=TRUE)
        }
    }
    ## get regions that are not introgressed anywhere
    if (!is.null(trueRegions)) {
        nullRegions <- allRegions
        for (trueIdx in names(trueRegions)) {
            nullRegions <- coverage.feat(allRegions, nullRegions, trueRegions[[trueIdx]], not=c(FALSE, FALSE, TRUE), get.feats=TRUE)
        }
        allTrue <- coverage.feat(allRegions, nullRegions, not=c(FALSE, TRUE), get.feats=TRUE)
    } else {
        nullRegions <- allRegions
        allTrue <- NULL
    }
    ## want just a TP rate and a FP rate for each pred.ind
    ## then average over inds with same pred
    rv <- data.frame()
    for (predIdx in names(pred)) {
        if (!grepl(".hom$", predIdx)) next
        mig <- strsplit(predIdx, ".", fixed=TRUE)[[1]][1]
        ind <- strsplit(predIdx, ".", fixed=TRUE)[[1]][2]
        hetIdx <- sprintf("%s.%s.het", mig, ind)
        if (!is.element(hetIdx, names(pred))) stop("Error")
        trueEl <- if (is.null(trueRegions)) NULL else trueRegions[[sprintf("%s.%s", mig, ind)]]

        ## deal with situation when same region is called both hom and het - this mostly
        ## happens when cutoff is zero - but can potentially happen when cutoff <= 0.5
        ## default to homozygous
        pred[[hetIdx]] <- coverage.feat(allRegions, pred[[hetIdx]], pred[[predIdx]], not=c(FALSE, FALSE, TRUE), get.feats=TRUE)
        bothpred <- coverage.feat(pred[[predIdx]], pred[[hetIdx]], get.feats=TRUE, or=TRUE)
        if (is.null(trueEl)) {
            tp <- 0
            fn <- 0
        } else {
            tp <- coverage.feat(trueEl, bothpred)
            fn <- coverage.feat(allRegions, trueEl, bothpred, not=c(FALSE, FALSE, TRUE))
        }

        fp <- 0
        tmphet <- if (is.null(allTrue)) pred[[hetIdx]] else overlap.feat(pred[[hetIdx]], allTrue, overlapping=FALSE)
        if (!is.null(tmphet)) fp <- fp + coverage.feat(tmphet, nullRegions)*0.5
        tmphom <- if (is.null(allTrue)) pred[[predIdx]] else overlap.feat(pred[[predIdx]], allTrue, overlapping=FALSE)
        if (!is.null(tmphom)) fp <- fp + coverage.feat(tmphom, nullRegions)
        tn <- coverage.feat(allRegions, nullRegions, bothpred, not=c(FALSE, FALSE, TRUE)) + coverage.feat(pred[[hetIdx]], nullRegions)*0.5
        tmprv <- data.frame(mig=mig, ind=ind, tp=tp, fn=fn, fp=fp, tn=tn)
        rv <- rbind(rv, tmprv)
    }
    rv$tprate <- rv$tp/(rv$tp+rv$fn)
    rv$fprate <- rv$fp/(rv$fp + rv$tn)
    ## want a confusion matrix showing how often we classify one type of migration for another
    combinedTrue <- list()
    for (trueMig in names(trueRegions)) {
        mig <- strsplit(trueMig, ".", fixed=TRUE)[[1]][1]
        if (is.null(combinedTrue[[mig]])) {
            combinedTrue[[mig]] <- trueRegions[[trueMig]]
        } else combinedTrue[[mig]] <- coverage.feat(combinedTrue[[mig]], trueRegions[[trueMig]], get.feats=TRUE, or=TRUE)
    }
    oldCombinedTrue <- combinedTrue
    ## now remove overlapping regions
    for (mig1 in names(combinedTrue)) {
        for (mig2 in names(combinedTrue)) {
            if (mig1 == mig2) next
            combinedTrue[[mig1]] <- coverage.feat(allRegions, combinedTrue[[mig1]], oldCombinedTrue[[mig2]], not=c(FALSE, FALSE, TRUE), get.feats=TRUE)
        }
    }
    combinedTrue[["none"]] <- nullRegions
    ## (OK there don't seem to be overlapping regions anyway, I removed them earlier when I
    ## created the trueRegions..)
    combinedPred <- list()
    predMigs <- c()
    notCalled <- allRegions
    for (str in names(pred)) {
        if (str == "none") next
        mig <- strsplit(str, ".", fixed=TRUE)[[1]][1]
        ind <- strsplit(str, ".", fixed=TRUE)[[1]][2]
        het <- strsplit(str, ".", fixed=TRUE)[[1]][3]
        str2 <- sprintf("%s.%s", mig, het)
        if (is.null(combinedPred[[str2]])) {
            combinedPred[[str2]] <- pred[[str]]
        } else {
            combinedPred[[str2]] <- coverage.feat(combinedPred[[str2]], pred[[str]], or=TRUE, get.feats=TRUE)
        }
        notCalled <- coverage.feat(allRegions, notCalled, combinedPred[[str2]], not=c(FALSE, FALSE, TRUE), get.feats=TRUE)
        predMigs <- unique(c(predMigs, mig))
    }
    ## now make sure homo/het don't overlap
    for (mig in predMigs) {
        hetstr <- sprintf("%s.het", mig)
        homstr <- sprintf("%s.hom", mig)
        combinedPred[[hetstr]] <- coverage.feat(allRegions, combinedPred[[hetstr]], combinedPred[[homstr]], not=c(FALSE, FALSE, TRUE), get.feats=TRUE)
    }
    confusionMat <- matrix(nrow=length(combinedTrue), ncol=length(predMigs)+1)
    rownames(confusionMat) <- names(combinedTrue)
    colnames(confusionMat) <- c(predMigs,"none")
    for (truemig in names(combinedTrue)) {
        confusionMat[truemig,"none"] <- 0
        for (predmig in predMigs) {
            confusionMat[truemig, predmig] <- (
                coverage.feat(combinedTrue[[truemig]], combinedPred[[sprintf("%s.hom", predmig)]]) + coverage.feat(combinedTrue[[truemig]], combinedPred[[sprintf("%s.het", predmig)]]))
        }
        confusionMat[truemig,"none"] <- coverage.feat(combinedTrue[[truemig]], notCalled)
    }
    rv2 <- data.frame()
    for (mig in unique(rv$mig)) {
        w <- which(rv$mig == mig)
        tprate <- mean(rv[w,"tprate"])
        fprate <- mean(rv[w,"fprate"])
        if (mig == "hToD" || mig == "aToD") migname <- "humToDen"
        if (mig == "hToN" || mig == "aToN") migname <- "humToNea"
        if (mig == "sToD") migname <- "supToDen"
        if (mig == "sToA" || mig=="sToH") migname <- "supToAfr"
        if (mig == "sToN") migname <- "supToNea"
        rv2 <- rbind(rv2, data.frame(mig=mig, tprate=tprate, fprate=fprate))
    }
    rv2

    list(summary=rv2, all=rv, confusionMat=confusionMat)
}



## want to return coverage for each mig, each individual, split into autosome/X,
## het/hom. Report averages per haplotype (so divide het coverage by 2)
summarizeRealResults <- function(dir, cutoff=0.5, chrs=c("autosome", "X"), dir2=NULL) {
    rv <- data.frame()
    covfile <- sprintf("%s/covered.bed", dir)
    if (!file.exists(covfile)) stop("Error cov not found")
    tmpcov <- read.feat(covfile)
    if (!is.null(dir2)) {
        covfile2 <- sprintf("%s/covered.bed", dir2)
        if (!file.exists(covfile2)) stop("Error cov2 not found")
        tmpcov2 <- read.feat(covfile2)
        tmpcov <- coverage.feat(tmpcov2, tmpcov, get.feats=TRUE)
    }
    cov <- list()
    for (chr in chrs) {
        if (chr == "autosome") {
            cov[[chr]] <- tmpcov[tmpcov$seqname != "X",]
        } else if (chr == "all") {
            cov[[chr]] <- tmpcov
        } else cov[[chr]] <- tmpcov[tmpcov$seqname == chr,]
    }
    files <- c(list.files(dir, pattern=sprintf(".*.het.%.1f.bed", cutoff), full.names=TRUE),
               list.files(dir, pattern=sprintf(".*.hom.%.1f.bed", cutoff), full.names=TRUE),
               list.files(dir, pattern=sprintf(".*.either.%.1f.bed", cutoff), full.names=TRUE))
    print(files)
    for (file in files) {
        tmp <- Sys.readlink(file)
        if (tmp != "") next
        bn <- basename(file)
        mig <- strsplit(bn, ".", fixed=TRUE)[[1]][1]
        ind <- strsplit(bn, ".", fixed=TRUE)[[1]][2]
        hetstr <- strsplit(bn, ".", fixed=TRUE)[[1]][3]
        tmp <- flatten.feat(read.feat(file))
        if (!is.null(dir2)) {
            file2 <- sprintf("%s/%s", dir2, bn)
            if (!file.exists(file2)) {
                cat("Cannot find ", file2, "\n")
            } else {
                tmp2 <- read.feat(file2)
                tmp <- flatten.feat(coverage.feat(tmp, tmp2, get.feats=TRUE))
            }
        }
        region <- list()
        for (chr in chrs) {
            if (chr == "autosome") {
                region[[chr]] <- tmp[tmp$seqname != "X",]
            } else if (chr == "all") {
                region[[chr]] <- tmp
            } else region[[chr]] <- tmp[tmp$seqname == chr,]
        }
        for (chr in chrs) {
            if (is.null(region[[chr]])) {
                meancov <- 0
            } else {
                meancov <- coverage.feat(region[[chr]])/coverage.feat(cov[[chr]])
            }
            rv <- rbind(rv, data.frame(mig=mig, ind=ind, het=hetstr, chr=chr, meancov=meancov,
                                       numel=nrow(region[[chr]]),
                                       meanlen=coverage.feat(region[[chr]])/nrow(region[[chr]]),
                                       medlen=median(lengths.feat(region[[chr]])), stringsAsFactors=FALSE))
        }
    }
    rv
}


makeBarPlotFromSummary <- function(x, stat="meancov", cex=0.9, add=FALSE, groupChrs=FALSE, doLegend=FALSE, ...) {
    ## I've been inconsistent about using "a" for African or "h" for human in these models where the
    ## only humans are Africans.. so change a's to h's if necessary
    x$mig <- ifelse(x$mig=="aToN", "hToN",
             ifelse(x$mig=="aToD", "hToD", as.character(x$mig)))
    x$migInd <- sprintf("%s.%s", x$mig, x$ind)
    migInds <- unique(x$migInd)
    chrs <- unique(x$chr)
    cov <- matrix(ncol=length(migInds), nrow=length(chrs))
    rownames(cov) <- unique(x$chr)
    colnames(cov) <- migInds
    covhom <- cov
    for (migind in migInds) {
        for (chr in chrs) {
            whet <- which(x$migInd==migind & x$chr==chr & x$het=="het")
            if (length(whet) != 1) stop("error1", whet, migind, chr, sep=" ")
            whom <- which(x$migInd==migind & x$chr==chr & x$het=="hom")
            if (length(whom) != 1) stop("Error2")
            if (stat == "meanlen") {
                x[x$numel==0,"meanlen"] <- 0
                cov[chr,migind] <- (x[whet,stat]*x[whet,"numel"]+x[whom,stat]*x[whom,"numel"])/(x[whet,"numel"]+x[whom,"numel"])
            } else if (stat == "meancov" || stat == "numel") {
                ## this should work for meancov or numel
                cov[chr,migind] <- x[whet,stat]*0.5+x[whom,stat]
                covhom[chr,migind] <- x[whom,stat]
            } else cat("Do not know stat ", stat)
        }
    }
    migs <- sapply(strsplit(migInds,".", fixed=TRUE), function(x) {x[1]})
    inds <- sapply(strsplit(migInds,".", fixed=TRUE), function(x) {x[2]})
    cols <- as.character(migCol[migs])
    if (stat == "meancov") cols <- makeTransparent(cols)
    cols <- sapply(cols, rep, length(chrs))
    ylabs <- list(meancov="Fraction introgressed",
                  numel="Number of elements",
                  meanlen="Mean length")
    ylim <- c(0, max(cov, na.rm=TRUE))
    if (ylim[2] > 0.038 && ylim[2] < 0.040) ylim[2] <- 0.040
    indName <- ifelse(inds=="den", "Den",
               ifelse(inds=="alt", "Alt",
               ifelse(inds=="vin", "Vin",
               ifelse(inds=="afr1", "San",
               ifelse(inds=="afr2", "Mandenka",
               ifelse(inds=="Basque_2F", "Basque",
               ifelse(inds=="Khomani_San_1F", "San",
               ifelse(inds=="Mandenka_2F", "Mandenka",
               ifelse(inds=="Mende_2F", "Mende",
               ifelse(inds=="Luhya_1F", "Luhya",
               ifelse(inds=="BantuKenya_2F", "BantuKenya",
               ifelse(inds=="Gambian_2F", "Gambian",
               ifelse(inds=="Mbuti_2F", "Mbuti",
               ifelse(inds=="Yoruba_1F", "Yoruba",
               ifelse(inds=="Papuan_1F", "Papuan",
               ifelse(inds=="Denisova", "Den",
               ifelse(inds=="Altai", "Alt",
               ifelse(inds=="Vindija", "Vin", "???"))))))))))))))))))
    migStart <- substr(migs, 1, 1)
    migSource <- ifelse(migStart=="h", "Hum",
                     ifelse(migStart=="s", "Sup",
                     ifelse(migStart=="d", "Den",
                     ifelse(migStart=="n", "Nea", "???"))))
    cat(inds, "\n")
    spacing <- numeric()
    if (!groupChrs) {
        for (i in 1:length(migs)) {
            if (i==1 || migs[i] == migs[i-1]) {
                spacing <- c(spacing, 0.5,rep(0, length(chrs)-1))
            } else spacing <- c(spacing, 1.0, rep(0,0, length(chrs)-1))
        }
        if (length(chrs) == 2) density <- c(NA, 30) else density <- rep(c(NA, 20, 40), length.out=length(chrs))
        print(cov)
        bp <- barplot(cov, beside=TRUE,
                      col=cols, space=spacing, density=density,
                      ylab=ylabs[[stat]], yaxs="i", names.arg=rep("", ncol(cov)), ylim=ylim, add=add, yaxt=if(add) {"n"} else NULL)
        if (stat=="meancov") {
            barplot(covhom, beside=TRUE, add=TRUE,
                    col=cols, space=spacing, density=density,
                    ylab=ylabs[[stat]], yaxs="i", names.arg=rep("", ncol(cov)), yaxt="n")
        }
        if (stat == "meancov" && !add)
            abline(h=seq(from=0, to=0.1, by=0.01), lty=3)
    } else {
        density <- rep(NA, length(migInds))
        nextDens <- 30
        for (i in 2:length(migInds)) {
            if (migs[i] == migs[i-1]) {
                density[i] <- nextDens
                nextDens <- nextDens+30
            } else nextDens <- 30
        }
        bp <- barplot(t(cov), beside=TRUE, col=t(cols), ylab=ylabs[[stat]], ylim=ylim, density=density, ...)
        barplot(t(covhom), beside=TRUE, add=TRUE, col=t(cols), yaxt="n", ylab="", density=density,...)
        tmp <- cov
        if (stat == "meancov" && !add)
            abline(h=seq(from=0, to=0.1, by=0.01), lty=3)
        if (stat == "meanlen" && !add)
            abline(h=seq(from=50000, to=1e6, by=50000), lty=3)
        if (doLegend) {
            legend(x="topleft", sprintf("%sTo%s", migSource, indName),
                   fill=as.character(migCol[migs]), density=density, ncol=2, bg="white")
        }
        return(list(bp=bp, mig=migSource, ind=indName))
    }
    if (!add) {
        mtext(text=sprintf("%s->\n%s", migSource, indName),
              at=apply(bp, 2, mean), line=-0.5, side=1, padj=1, cex=cex)
    }
    invisible(bp)
}
