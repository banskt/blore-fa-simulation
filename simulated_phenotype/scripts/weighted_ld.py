import numpy as np
import os
import argparse

def parse_args():

    parser = argparse.ArgumentParser(description='Calculate genomic inflation factor in the simulations')

    parser.add_argument('-l', '--locusprefix',
                        type=str,
                        dest='locusprefix',
                        metavar='STR',
                        help='prefix of locusname')

    parser.add_argument('-st', '--studies',
                        nargs='*',
                        type=str,
                        dest='studies',
                        metavar='STR',
                        help='list of study names')

    parser.add_argument('-stn', '--samples',
                        nargs='*',
                        type=int,
                        dest='samples',
                        metavar='INT',
                        help='list of numbers of samples in studies')

    parser.add_argument('-o', '--outdir',
                        type=str,
                        dest='outdir',
                        metavar='DIR',
                        help='path of the output directory')

    parser.add_argument('-fmeta', '--metafile',
                        type=str,
                        dest='metafile',
                        metavar='FILE',
                        help='file containing results of META for the given locus')

    parser.add_argument('-dld', '--lddir',
                        type=str,
                        dest='lddir',
                        metavar='DIR',
                        help='path of the directory containing LD of all studies')


    opts = parser.parse_args()
    return opts


opts = parse_args()
studies = opts.studies
samples = opts.samples
locusprefix = opts.locusprefix
lddir = opts.lddir
outdir = opts.outdir
metafile = opts.metafile

if not os.path.exists(outdir):
    os.makedirs(outdir)

# Order LD matrix according to results of META
allrsid = list()
with open(metafile, 'r') as mfile:
    next(mfile)
    for mline in mfile:
        allrsid.append(mline.split()[1].strip())

snptot = len(allrsid)

ld = [[] for study in studies]
for i, study in enumerate(studies):

    ldfile = os.path.join ( lddir, study, locusprefix + '.LD' )
    ldmatrix = np.loadtxt(ldfile)
    studyrsid = list()
    studymapfile = os.path.join( lddir, study, locusprefix + '.meta' )
    with open(studymapfile, 'r') as mfile:
        next(mfile)
        for line in mfile:
            studyrsid.append(line.split()[1].strip())

    rsq = np.zeros([snptot, snptot], dtype=np.float64)
    for snp in range(snptot):
        rsq[snp, snp] = 1.0

    for snp1 in range(snptot - 1):
        for snp2 in range(snp1 + 1, snptot):
            rsid1 = allrsid[snp1]
            rsid2 = allrsid[snp2]
            if all(x in studyrsid for x in [rsid1, rsid2]):
                indx1 = studyrsid.index(rsid1)
                indx2 = studyrsid.index(rsid2)
                rsq[snp1, snp2] = ldmatrix[indx1, indx2]
            else:
                rsq[snp1, snp2] = 0.0
            rsq[snp2, snp1] = rsq[snp1, snp2]
    ld[i] = rsq

ldmod = np.zeros_like(ld[0])
for i in range(len(studies)):
    ldmod += ld[i] * samples[i]

ldweighted = ldmod / sum(samples)
outfilename = os.path.join(outdir, locusprefix + '.LD') 
np.savetxt(outfilename, ldweighted)

rsidfilename = os.path.join(outdir, locusprefix + '.LD.rsid')
with open(rsidfilename, 'w') as mfile:
    for rsid in allrsid:
        mfile.write('{:s}\n'.format(rsid))
