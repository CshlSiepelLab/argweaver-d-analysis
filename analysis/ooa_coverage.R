## Code to generate figure 4 (coverage of CRF and ARGweaver calls in modern humans)

source("../scripts/functions.R")
require("Hmisc")
require("argweaver")
require("rphast")


crf <- list()
crf$Basque <- list()
crf$Basque$nea.homo <- read.feat("CRF_results/S_Basque-2/neandertal.set1.homo.bed")
crf$Basque$nea.het <- read.feat("CRF_results/S_Basque-2/neandertal.set1.het.bed")
crf$Basque$den.homo <- read.feat("CRF_results/S_Basque-2/denisova.set2.homo.bed")
crf$Basque$den.het <- read.feat("CRF_results/S_Basque-2/denisova.set2.het.bed")
crf$Basque$regions <- read.feat("hg19_regions.bed")
crf$Papuan <- list()
crf$Papuan$nea.homo <- read.feat("CRF_results/S_Papuan-1/neandertal.set2.homo.bed")
crf$Papuan$nea.het <- read.feat("CRF_results/S_Papuan-1/neandertal.set2.het.bed")
crf$Papuan$den.homo <- read.feat("CRF_results/S_Papuan-1/denisova.set2.homo.bed")
crf$Papuan$den.het <- read.feat("CRF_results/S_Papuan-1/denisova.set2.het.bed")
crf$Papuan$regions <- read.feat("hg19_regions.bed")


x <- list()
x$Basque <- list()
x$Basque$nea.homo <- read.feat("ooaModel_nosuper_1fra2afr/regions2/nToH.Basque_2F.hom.0.5.bed")
x$Basque$nea.het <- read.feat("ooaModel_nosuper_1fra2afr/regions2/nToH.Basque_2F.het.0.5.bed")
x$Basque$den.homo <- read.feat("ooaModel_nosuper_1fra2afr/regions2/dToH.Basque_2F.hom.0.5.bed")
x$Basque$den.het <- read.feat("ooaModel_nosuper_1fra2afr/regions2/dToH.Basque_2F.het.0.5.bed")
x$Basque$regions <- read.feat("ooaModel_nosuper_1fra2afr/regions2/covered.bed")
x$Papuan <- list()
x$Papuan$nea.homo <- read.feat("ooaModel_nosuper_1pap2afr/regions2/nToH.Papuan_1F.hom.0.5.bed")
x$Papuan$nea.het <- read.feat("ooaModel_nosuper_1pap2afr/regions2/nToH.Papuan_1F.het.0.5.bed")
x$Papuan$den.homo <- read.feat("ooaModel_nosuper_1pap2afr/regions2/dToH.Papuan_1F.hom.0.5.bed")
x$Papuan$den.het <- read.feat("ooaModel_nosuper_1pap2afr/regions2/dToH.Papuan_1F.het.0.5.bed")
x$Papuan$regions <- read.feat("ooaModel_nosuper_1pap2afr/regions2/covered.bed")
x$San <- list()
x$San$nea.homo <- read.feat("ooaModel_nosuper_1pap2afr/regions2/nToH.Khomani_San_1F.hom.0.5.bed")
x$San$nea.het <- read.feat("ooaModel_nosuper_1pap2afr/regions2/nToH.Khomani_San_1F.het.0.5.bed")
x$San$den.homo <- read.feat("ooaModel_nosuper_1pap2afr/regions2/dToH.Khomani_San_1F.hom.0.5.bed")
x$San$den.het <- read.feat("ooaModel_nosuper_1pap2afr/regions2/dToH.Khomani_San_1F.het.0.5.bed")
x$San$regions <- read.feat("ooaModel_nosuper_1pap2afr/regions2/covered.bed")
x$Mandenka <- list()
x$Mandenka$nea.homo <- read.feat("ooaModel_nosuper_1fra2afr/regions2/nToH.Mandenka_2F.hom.0.5.bed")
x$Mandenka$nea.het <- read.feat("ooaModel_nosuper_1fra2afr/regions2/nToH.Mandenka_2F.het.0.5.bed")
x$Mandenka$den.homo <- read.feat("ooaModel_nosuper_1fra2afr/regions2/dToH.Mandenka_2F.hom.0.5.bed")
x$Mandenka$den.het <- read.feat("ooaModel_nosuper_1fra2afr/regions2/dToH.Mandenka_2F.het.0.5.bed")
x$Mandenka$regions <- read.feat("ooaModel_nosuper_1fra2afr/regions2/covered.bed")

usedata <- list(ARG=x,CRF=crf)
##cols <- list(nea=rgb(1,0,0,0.25), den=rgb(0,0,1,0.25))
cols <- list(nea="red", den="blue")
mat <- matrix(nrow=2, ncol=24)
rownames(mat) <- c("both", "only")
colnames(mat) <- sprintf("c%i", 1:24)
colmat <- matrix(nrow=2, ncol=24)
rownames(colmat) <- c("homo", "het")
idx <- 1
for (ind in c("Basque", "Papuan", "Mandenka", "San")) {
    for (chr in c("autosome", "X")) {
        for (arc in c("nea", "den")) {
             if (! is.null(usedata$CRF[[ind]])) {
                 tmp <- list()
                 tmp$het1 <- usedata$ARG[[ind]][[sprintf("%s.het", arc)]]
                 tmp$het2 <- usedata$CRF[[ind]][[sprintf("%s.het", arc)]]
                 tmp$hom1 <- usedata$ARG[[ind]][[sprintf("%s.homo", arc)]]
                 tmp$hom2 <- usedata$CRF[[ind]][[sprintf("%s.homo", arc)]]
                 tmp$cov <- usedata$ARG[[ind]]$regions
                 if (chr == "autosome") {
                     for (i in names(tmp)) {tmp[[i]] <- tmp[[i]][tmp[[i]]$seqname!="X",]}
                 } else {
                     for (i in names(tmp)) {tmp[[i]] <- tmp[[i]][tmp[[i]]$seqname=="X",]}
                 }
                 for (i in names(tmp)) {
                     if (i != "cov") tmp[[i]] <- coverage.feat(tmp[[i]], tmp$cov, get.feats=TRUE)
                 }
                 mat["both",idx] <- (coverage.feat(tmp$hom1, tmp$hom2, tmp$cov)+
                                     coverage.feat(tmp$het1, tmp$het2, tmp$cov)/2+
                                     coverage.feat(tmp$hom1, tmp$het2, tmp$cov)/2+
                                     coverage.feat(tmp$het1, tmp$hom2, tmp$cov)/2)/coverage.feat(tmp$cov)
                 mat["both",idx+1] <- mat["both",idx]
                 mat["only",idx] <- (coverage.feat(tmp$cov, tmp$hom1, tmp$het2, tmp$hom2, not=c(FALSE, FALSE, TRUE, TRUE)) +
                                     coverage.feat(tmp$cov, tmp$het1, tmp$het2, tmp$hom2, not=c(FALSE, FALSE, TRUE, TRUE))/2)/coverage.feat(tmp$cov)
                 mat["only",idx+1] <- (coverage.feat(tmp$cov, tmp$hom2, tmp$het1, tmp$hom1, not=c(FALSE, FALSE, TRUE, TRUE)) +
                                     coverage.feat(tmp$cov, tmp$het2, tmp$het1, tmp$hom1, not=c(FALSE, FALSE, TRUE, TRUE))/2)/coverage.feat(tmp$cov)
                 colnames(mat)[idx] <- sprintf("%s.%s.%s.ARG", ind, arc, chr)
                 colnames(mat)[idx+1] <- sprintf("%s.%s.%s.CRF", ind, arc, chr)
                 idx <- idx+2
             }  else {
                 for (method in c("ARG", "CRF")) {
                     if (is.null(usedata[[method]][[ind]])) next
                     ## use ARG here instead of method because ARG coverage is strictly less than crf
                     cov <- usedata[["ARG"]][[ind]]$regions
                     het <- usedata[[method]][[ind]][[sprintf("%s.het", arc)]]
                     hom <- usedata[[method]][[ind]][[sprintf("%s.homo", arc)]]
                     if (chr == "autosome") {
                         cov <- cov[cov$seqname != "X",]
                     } else {
                         cov <- cov[cov$seqname == "X",]
                     }
                     mat["only", idx] <- (coverage.feat(hom,cov)+0.5*coverage.feat(het,cov))/coverage.feat(cov)
                     mat["both",idx] <- 0
                     colnames(mat)[idx] <- sprintf("%s.%s.%s.%s", ind, arc, chr, method)
                     idx <- idx+1
                 }
             }
         }
    }
}


## figure 4
pdf("ooaBarPlot.pdf", width=5, height=3.5)

if (names(dev.cur()) != "pdf") {
    dev.off()
    x11(width=5, height=3.5)
}
par(mar=c(4,4,1,1))
add <- FALSE
spacing <- c(1,0,0,0,  # basque/autosome
             0.5,0,0,0,  # basque/X
             1,0,0,0,   #pap/autosome
             0.5,0,0,0,       #pap/X
             1,0,0.5,0,1,0,0.5,0)
argvals <- numeric(ncol(mat))
crfvals <- numeric(ncol(mat))
for (i in 1:ncol(mat)) {
    if (grepl("ARG", colnames(mat)[i])) {
        argvals[i] <- mat["both",i] + mat["only",i]
        crfvals[i] <- mat["both",i]
    } else {
        argvals[i] <- mat["both",i]
        crfvals[i] <- mat["both",i]+mat["only",i]
    }
}
usecols <- ifelse(grepl(".nea.", colnames(mat), fixed=TRUE), cols[["nea"]],cols[["den"]])
bp <- barplot(argvals, col=makeTransparent(usecols), space=spacing, xaxt="n", yaxt="n")
bp <- barplot(crfvals, col=usecols, angle=30, density=25, add=TRUE, space=spacing)
axis(side=2, cex=0.7)
mtext("Coverage", side=2, line=2.5, cex=1.1)
mtext(c("Basque", "Papuan", "Mandenka", "San"), side=1, line=1, cex=0.9, 
      at=c(mean(bp[1:8,]),mean(bp[9:16,]),mean(bp[17:20,]),mean(bp[21:24,])))
mtext(c("autosome","X"), side=1, line=0.0, cex=0.7,
      at=c(mean(bp[1:4,]),mean(bp[5:8,]),mean(bp[9:12,]),mean(bp[13:16,]),mean(bp[17:18,]),mean(bp[19:20,]),mean(bp[21:22,]),mean(bp[23:24,])))
abline(h=seq(from=0, to=0.04, by=0.005), lty=3)
legend(x="topright", legend=c("ARG-Nea", "CRF-Nea", "ARG-Den", "CRF-Den"),
       fill=c(makeTransparent(cols$nea), cols$nea, makeTransparent(cols$den), cols$den),
       angle=30, density=c(NA,25,NA,25), bg="white")
if (names(dev.cur())=="pdf") dev.off()

