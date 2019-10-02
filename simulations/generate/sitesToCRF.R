require("argweaver")
outgroup <- "chimp"
targetStr <- "eur"
meanRho <- 1.358879e-8
args <- commandArgs(trailingOnly=TRUE)
indir <- args[1]


## x is subset of sites$sites array with individuals we want to output
## will group genotypes by adjacent individuals
writeGenoFileOld <- function(x, outgroup, file) {
    str <- rep("", nrow(x))
    for (i in seq(from=1, to=ncol(x), by=2)) {
        geno <- ifelse(x[,i]=="N" | x[,i+1]=="N", 9,
                       (x[,i]!=outgroup) + (x[,i+1]!=outgroup))
        str <- sprintf("%s%d", str, geno)
    }
    write(str, file=file, sep="\n", append=FALSE)
    invisible(NULL)
}

writeGenoFile <- function(x, outgroup, file) {
    str <- rep("", nrow(x))
    for (i in 1:ncol(x)) {
        geno <- ifelse(x[,i]=="N", 9,
                       x[,i]!=outgroup)
        str <- sprintf("%s%d", str, geno)
    }
    write(str, file=file, sep="\n", append=FALSE)
    invisible(NULL)
}

## x is subset of sites$sites array with individuals we want to output
## will group genotypes by adjacent individuals
writeFreqFile <- function(x, outgroup, file) {
    count <- rep(0, nrow(x))
    n <- rep(0, nrow(x))
    for (i in 1:ncol(x)) {
        count <- count + (x[,i] != "N" & x[,i] != outgroup)
        n <- n + ( x[,i] != "N")
    }
    write.table(data.frame(snp=sprintf("snp%i", 1:nrow(x)),
                           freq=sprintf("%.3f", count/n)),
                row.names=FALSE, col.names=FALSE, file=file, append=FALSE, sep="\t",
                quote=FALSE)
    invisible(NULL)
}




outdir <- sprintf("%s/crf", indir)
system(sprintf("mkdir -p %s", outdir))

x <- readSites(sprintf("%s/sim.sites.gz", indir))

## first, remove sites that are only invariant in chimp. Also remove very rare
## sites with > 2 alleles

outcol <- which(names(x$sites)==outgroup)
numAllele <- apply(x$sites[,-outcol], 1, function(x) {length(table(as.character(x)))})
f <- numAllele == 2
x$sites <- x$sites[f,]
x$pos <- x$pos[f]
## now ref is chimp allele

## write ind file for admixed (target) individuals
indfile <- sprintf("%s/%s.ind", outdir, targetStr)
targetInds <- grep(targetStr, names(x$sites), value=TRUE)
df <- data.frame(ind=targetInds, gender="U", status=targetStr)
#f1 <- seq(from=1, to=nrow(df), by=2)
#df$ind <- sprintf("%s.%s", df[f1,"ind"], df[f1+1,"ind"])
#df <- df[f1,]
write.table(df, file=indfile, quote=FALSE, row.names=FALSE,
            sep="\t", col.names=FALSE)

## need to get physical coordinate
map <- read.table(sprintf("%s/recomb_map.bed.gz", indir))
names(map) <- c("chrom", "chromStart", "chromEnd", "rho")
tmp <- c(0,map$rho*(map$chromEnd - map$chromStart))
tmp <- tmp[-length(tmp)]
map$pos <- cumsum(tmp)/meanRho/1e8
map <- map[map$chromEnd - map$chromStart >= 1,]
tmp <- sapply(x$pos, function(x) {tail(which(x >= map$chromStart), n=1)})
genpos <- map[tmp,"pos"] + (x$pos - map[tmp,"chromStart"])*map[tmp,"rho"]/meanRho/1e8

write.table(data.frame(id=sprintf("snp%i", 1:length(x$pos)),
                       chrom=x$region$chrom,
                       genpos=sprintf("%.4g", genpos),
                       pos=x$pos),
            file=sprintf("%s/snp", outdir), quote=FALSE, row.names=FALSE,
            col.names=FALSE, sep="\t")

afrCols <- grepl("afr", names(x$sites))
writeGenoFile(x$sites[,afrCols],
              outgroup=x$sites[,outgroup],
              file=sprintf("%s/afr.geno", outdir))
writeFreqFile(x$sites[,afrCols], outgroup=x$sites[,outgroup],
              file=sprintf("%s/afr.freq", outdir))

writeGenoFile(x$sites[,grepl(targetStr,names(x$sites))],
              outgroup=x$sites[,outgroup],
              file=sprintf("%s/%s.geno", outdir, targetStr))

neaCols <- c("vin1", "vin2", "alt1", "alt2")
writeGenoFile(x$sites[,neaCols],
              outgroup=x$sites[,outgroup],
              file=sprintf("%s/nea.geno", outdir))
writeFreqFile(x$sites[,neaCols], outgroup=x$sites[,outgroup],
              file=sprintf("%s/nea.freq", outdir))


chr <- as.character(x$region$chr)
system(sprintf("sed 's/CHR/%s/g' par.caller.test > %s/par.caller.test", chr, outdir))
    
system(sprintf("cp -f mcle %s", outdir))


