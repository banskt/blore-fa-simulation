library(R2BGLiMS)
library(argparse)

parser <- ArgumentParser(description='Process some arguments')

parser$add_argument('--meta',
                    dest='metadir',
                    metavar='DIR',
                    type="character",
                    help='output directory of META')

parser$add_argument('--locusnames',
                    dest='locusnames',
                    metavar='FILE',
                    type="character",
                    help='name of file containing locus prefixes')

parser$add_argument('--locidir',
                    dest='locidir',
                    metavar='DIR',
                    type="character",
                    help='name of the directory containing reference genotype')

parser$add_argument('--lddir',
                    dest='lddir',
                    metavar='DIR',
                    type="character",
                    help='name of the directory containing reference LD matrices')

parser$add_argument('--outdir',
                    dest='outdir',
                    metavar='DIR',
                    type="character",
                    help='name of the output directory')

parser$print_help()
args = parser$parse_args(commandArgs(TRUE))

metadir = args$metadir
locusnames_file = args$locusnames
locidir = args$locidir
lddir = args$lddir
outdir = args$outdir

locusprefixes <- scan(locusnames_file, what="", sep="\n")

snps.regions = list()
Xref.regions = list()
marginal.betas = c()

for (locusprefix in locusprefixes) {
    metafile = sprintf("%s/%s.%s", metadir, locusprefix, "meta.out")
    metares = read.table(metafile, header=TRUE)
    this.marginal.betas = metares$beta

    this.snp.rsids = as.character(metares$rsid)
    names(this.marginal.betas) = this.snp.rsids

    matgenfile = sprintf("%s/%s.matgen", locidir, locusprefix)
    matgen = read.table(matgenfile, header=FALSE, sep = ",")
    nsample = dim(matgen)[2]
    this.Xref = t(matgen[, 4:dim(matgen)[2]])
    colnames(this.Xref) = as.character(matgen[,1])
    rownames(this.Xref) = c()

    full.rank = JAM_RankCheck(this.Xref)

    if (! full.rank) {
        ldfile = sprintf("%s/%s.LD", lddir, locusprefix)
        ldmetafile = sprintf("%s/%s.meta", lddir, locusprefix)
        ldmatrix = read.table(ldfile, header=FALSE)
        choose.snps = as.character(read.table(ldmetafile, header=TRUE)$RSID)
    }

    # Keep reducing LD threshold till we get a full rank genotype matrix
    eps = 0.05
    ldthres = 1.0
    while (! full.rank) {
        prune.snps = c()
        for (i in 1:(length(choose.snps) - 1) ) {
            for (j in (i+1):length(choose.snps)) {
                if (ldmatrix[i, j] ** 2 >= ldthres){
                   prune.snps = c(prune.snps, choose.snps[j])
                }
            }
        }
        prune.snps = unique(prune.snps)
        new.choose.snps = setdiff(choose.snps, prune.snps)
        this.Xref = this.Xref[, new.choose.snps]
        ldthres = ldthres - eps
        full.rank = JAM_RankCheck(this.Xref)
    }

    Xref.regions[[locusprefix]] = this.Xref
    snps.regions[[locusprefix]] = colnames(this.Xref)
    marginal.betas = c(marginal.betas, this.marginal.betas[colnames(this.Xref)])
}

names(snps.regions) = locusprefixes
names(Xref.regions) = locusprefixes

# All regions results
all.regions.results <- JAM(marginal.betas = marginal.betas,
                           X.ref = Xref.regions,
                           model.space.prior = list("a" = 1,
                                                    "b" = sum(lengths(snps.regions)),
                                                    "Variables" = unlist(snps.regions, recursive = FALSE, use.names = FALSE)
                                                    ),
                           n.mil = 1,
                           n = dim(Xref.regions[[1]])[1]
                           )

# Write results
for (locusprefix in locusprefixes) {
    outfilename = sprintf("%s/%s.out", outdir, locusprefix)
    outdata = data.frame(all.regions.results@posterior.summary.table[snps.regions[[locusprefix]], c("PostProb", "BF")])
    outdata = cbind("rsid" = rownames(outdata), outdata)
    write.table (outdata, file=outfilename, sep="\t", quote=FALSE, row.names = FALSE)
}

outfilename = sprintf("%s/benchmark.times", outdir)
benchmark.times = data.frame(all.regions.results@run.times)
colnames(benchmark.times) = c("write", "bglims", "processing")
write.table(benchmark.times, file=outfilename, sep="\t", quote=FALSE, row.names = FALSE)
#elapsed_time <- results@run.times$write.time + results@run.times$results.processing.time + results@run.times$bglims.time
