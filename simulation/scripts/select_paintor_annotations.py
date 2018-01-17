from scipy import stats
import numpy as np
import collections
import argparse
import os

def parse_args():

    parser = argparse.ArgumentParser(description='Select top, roughly uncorrelated annotations for PAINTOR')

    parser.add_argument('--corr',
                        default='func_annot_corr.dat',
                        type=str,
                        dest='corr_filename',
                        metavar='FILE',
                        help='name of file for correlation of functional annotations')

    parser.add_argument('--names',
                        default='functional_annotation_names.dat',
                        type=str,
                        dest='fnames_filename',
                        metavar='FILE',
                        help='name of file for the names of functional annotations')

    parser.add_argument('--logBFdir',
                        type=str,
                        dest='annotbf_dir',
                        metavar='DIR',
                        help='directory of logBF from PAINTOR')

    parser.add_argument('--out',
                        default='selected_functional_annotation.txt',
                        type=str,
                        dest='outfilename',
                        metavar='FILE',
                        help='name of output file with names of selected functional annotations')

    parser.add_argument('--fmax',
                        default=5,
                        type=int,
                        dest='fmax',
                        metavar='INT',
                        help='number of functional annotations to be selected')

    parser.add_argument('--corrmax',
                        default=0.3,
                        type=float,
                        dest='corrmax',
                        metavar='real',
                        help='maximum allowed correlation')

    opts = parser.parse_args()
    return opts

opts = parse_args()
corr_filename = opts.corr_filename
fnames_filename = opts.fnames_filename
annotbf_dir = os.path.realpath(opts.annotbf_dir)
outfilename = opts.outfilename


# Dictionary of functional annotation correlation
fa_corr = np.loadtxt(corr_filename)
with open(fnames_filename, 'r') as mfile:
    next(mfile)
    fa_names = mfile.readlines()
fa_names = [x.strip() for x in fa_names]

fa_corr_dict = collections.defaultdict(lambda:0)
for i, a1 in enumerate(fa_names):
    for j, a2 in enumerate(fa_names):
        fa_corr_dict[a1, a2] = fa_corr[i,j]

# Dictionary of logBF for each annotation
logbf = collections.defaultdict(lambda:0)
filename = "%s/BF.%s" % (annotbf_dir, 'Base')
with open(filename, 'r') as mfile:
    logbf_h0 = float(mfile.readline().strip())

for a1 in fa_names:
    filename = "%s/BF.%s" % (annotbf_dir, a1)
    with open(filename, 'r') as mfile:
        logbf[a1] = float(mfile.readline().strip())

# Create a dictionary of p-values for each annotation
pval = collections.defaultdict(lambda:0)
for a1 in fa_names:
    LRT = -2 * (logbf_h0 - logbf[a1])
    pval[a1] = 1 - stats.chi2.cdf(LRT, 1)
    
# Ranked list of annotations
ranked_annot = sorted(pval, key=pval.get)

# Select top M annotations which are roughly uncorrelated
selected_annot = [ranked_annot[0]]
i = 1
while (len(selected_annot) < opts.fmax):
    next_annot = ranked_annot[i]
    # correlation with previous annotations
    corrlist = [fa_corr_dict[next_annot, x] for x in selected_annot]
    maxcorr = max(corrlist)
    if maxcorr < opts.corrmax:
        selected_annot.append(next_annot)
    i += 1

with open(outfilename, 'w') as mfile:
    for x in selected_annot:
        mfile.write("%s\n" % x)
