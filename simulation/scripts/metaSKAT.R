library(argparse)
library(MetaSKAT)
library(SKAT)

parser <- ArgumentParser(description='Process some arguments')

parser$add_argument('--studies',
                    dest='studynames',
                    metavar='COHORT',
                    type="character",
                    nargs='*',
                    help='list of study names')

parser$add_argument('--locusnames',
                    dest='locusnames',
                    metavar='FILE',
                    type="character",
                    help='name of file containing locus prefixes')

parser$add_argument('--outdir',
                    dest='outdir',
                    metavar='DIR',
                    type="character",
                    help='name of the output directory')

parser$add_argument('--outfile',
                    dest='outfile',
                    metavar='FILE',
                    type="character",
                    help='name of the output file')

parser$print_help()
args <- parser$parse_args(commandArgs(TRUE))

studynames <- args$studynames
nstudy <- length(studynames)
locusnames_file <- args$locusnames 
outdir <- args$outdir
outfile <- args$outfile
locusprefixes <- scan(locusnames_file, what="", sep="\n")

study_file_name <- function(outdir, study, locusprefix, ext) {
    filename <- sprintf("%s/%s/%s.%s", outdir, study, locusprefix, ext)
    return (filename)
}

file.create(outfile)

for (locusprefix in locusprefixes) {
    # Generate summary statistics for each study
    file.mat.vec <- rep("", nstudy)
    file.info.vec <- rep("", nstudy)
    
    for (study.indx in 1:nstudy) {
    
        study <- studynames[study.indx]
    
        bedfile   <- study_file_name(outdir, study, locusprefix, "bed")
        bimfile   <- study_file_name(outdir, study, locusprefix, "bim")
        famfile   <- study_file_name(outdir, study, locusprefix, "fam")
        setfile   <- study_file_name(outdir, study, locusprefix, "SetID")
        mssdfile  <- study_file_name(outdir, study, locusprefix, "MSSD")
        minfofile <- study_file_name(outdir, study, locusprefix, "MInfo")
    
        pheno <- read.table(famfile, header=FALSE)[,6]
        nsample <- length(pheno)
        null.obj <- SKAT_Null_Model(pheno ~ 1)
        Generate_Meta_Files(null.obj, bedfile, bimfile, setfile, mssdfile, minfofile, nsample)
        file.mat.vec[study.indx] <- mssdfile
        file.info.vec[study.indx] <- minfofile
    }
    
    # open files
    cohort.info <- Open_MSSD_File_2Read (file.mat.vec, file.info.vec)
    
    # get a p-value of the first set.
    #MetaSKAT_MSSD_OneSet (cohort.info, SetID=locusprefix)$p.value
    
    # get p-values of all sets
    result <- MetaSKAT_MSSD_ALL(cohort.info)
    
    resultline <- sprintf("%s %g", result$SetID, result$p.value)
    write(resultline, file=outfile, append=TRUE)
}
