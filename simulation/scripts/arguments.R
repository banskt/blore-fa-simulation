library(argparse)

parser <- ArgumentParser(description='Process some arguments')
parser$add_argument('--num',
                    dest='integers', 
                    metavar='N', 
                    type="integer", 
                    nargs='+',
                    help='integers for the accumulator')
parser$add_argument('--gen',
                    dest='genotype_files',
                    metavar='FILE',
                    type="character",
                    nargs='*',
                    help='list of genotype files')
parser$print_help()
# default args for ArgumentParser()$parse_args are commandArgs(TRUE)
# which is what you'd want for an Rscript but not for interactive use

#args <- parser$parse_args(c("--sum", "1", "2", "3"))
args <- parser$parse_args(commandArgs(TRUE))

integers <- args$integers
genfiles <- args$genotype_files
print(integers)
print(genfiles)
