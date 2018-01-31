import argparse
import glob
import random
import collections
import os
import numpy as np
import scipy.stats

SNPINFO_FIELDS = ['rsid', 'bp_location', 'alt_allele', 'ref_allele', 'studies']
class SnpInfo(collections.namedtuple('_SnpInfo', SNPINFO_FIELDS)):
    __slots__ = ()

def parse_args():

    parser = argparse.ArgumentParser(description='Pick random SNPs for simulation')

    parser.add_argument('-dl', '--locidir',
                        type=str,
                        dest='locidir',
                        metavar='DIR',
                        help='path of the loci dosages directory')

    parser.add_argument('-df', '--featdir',
                        type=str,
                        dest='featdir',
                        metavar='DIR',
                        help='path of the features directory')

    parser.add_argument('-ds', '--simdir',
                        type=str,
                        dest='simdir',
                        metavar='DIR',
                        help='path of the simulation directory')

    parser.add_argument('-fl', '--locusnames',
                        type=str,
                        dest='locusnames',
                        metavar='FILE',
                        help='file with prefix of locusnames')

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

    parser.add_argument('-o', '--outfile',
                        default='causal.snplist',
                        type=str,
                        dest='outfile',
                        metavar='FILE',
                        help='name of file for causal SNPs list')

    parser.add_argument('-l', '--lmax',
                        default=100,
                        type=int,
                        dest='ncausal_loci',
                        metavar='INT',
                        help='number of causal loci')

    parser.add_argument('-p', '--prop',
                        default=5,
                        type=float,
                        dest='proportion',
                        metavar='FLOAT',
                        help='fraction of causal SNPs')

    parser.add_argument('-hg', '--heritability',
                        default=0.25,
                        type=float,
                        dest='sigma_herited_sq',
                        metavar='real',
                        help='variance (sigma^2) of narrow heritability ')

    parser.add_argument('-k', '--prevalence',
                        default=0.3,
                        type=float,
                        dest='prevalence',
                        metavar='real',
                        help='disease prevalence')

    parser.add_argument('-t', '--simtype',
                        default='fixed',
                        type=str,
                        dest='simtype',
                        metavar='STR',
                        help='type of beta distribution to be used for phenotype -- fixed, unimodal, bimodal, studentsT')


    opts = parser.parse_args()
    return opts


def read_locus_prefixes(filename):
    locusprefixes = list()
    with open(filename, 'r') as mfile:
        for line in mfile:
            mline = line.split()
            locusprefixes.append(mline[0].strip())
    return locusprefixes

def read_snps(locidir, studies, locusprefixes):
    nstudy = len(studies)
    snpinfo = list()
    for study in studies:
        snpstudy = list()
        for locusprefix in locusprefixes:
            mapfile = os.path.join(locidir, study, locusprefix + '.map')
            snplocus = list()
            with open(mapfile, 'r') as mfile:
                for line in mfile:
                    mline = line.split()
                    this_snp = SnpInfo(rsid = mline[1],
                                       bp_location = int(mline[3]),
                                       alt_allele = mline[4],
                                       ref_allele = mline[5],
                                       studies = [])
                    snplocus.append(this_snp)
            snpstudy.append(snplocus)
        snpinfo.append(snpstudy)
    return snpinfo

def read_features(locusprefixes, snppool, featuredir):
    featurelist = list()
    # Assert number of loci is equal to number of input files
    # Assert number of features is same in all loci
    for l, locusprefix in enumerate(locusprefixes):
        featurefile = os.path.join(featuredir, locusprefix + '.dhs')
        this_snpinfo = snppool[l]
        with open(featurefile, 'r') as mfile:
            featnames = mfile.readline().split()[3:]
            nfeat = len(featnames) + 1
            features = np.zeros((nfeat, len(this_snpinfo)))
            features[0, :] = 1
            for line in mfile:
                linesplit = line.split()
                rsid = linesplit[0]
                bp   = int(linesplit[2])
                vect = np.array(linesplit[3:])
                snp  = [x for x in this_snpinfo if x.rsid == rsid and x.bp_location == bp]
                if len(snp) == 0:
                    print ("%s not present in these cohorts" % rsid)
                elif len(snp) == 1:
                    indx = this_snpinfo.index(snp[0])
                    features[1:, indx] = vect
                    #features[:, indx] = vect
                else:
                    print ("Error: More than one SNP with same rsid")
        featurelist.append(features)
    return featurelist

def norm_binom(gt, f):
    gt = (gt - (2 * f)) / np.sqrt(2 * f * (1 - f))
    return gt

def read_genotype(gendir, snpinfo, nsample):
    ncol = nsample * 3 + 5
    nloci = len(snpinfo)
    genotypelist = list()
    for l, locusprefix in enumerate(locusprefixes):
        filename = os.path.join(gendir, locusprefix + '.gen')
        dosage = np.loadtxt(filename, usecols=range(5, ncol)).T
        nsnps = len(snpinfo[l])
        genotype = np.empty([nsample, nsnps], dtype=np.float64)
        for j in range(nsample):
            maj_ind = j * 3            # Allele AA, where A is the alternate allele
            het_ind = maj_ind + 1      # Allele AB, where B is the reference allele
            min_ind = het_ind + 1      # Allele BB
            genotype[j,:] = 2 * dosage[min_ind, :] + dosage[het_ind, :] # [AA, AB, BB] := [0, 1, 2]
        freq  = np.sum(genotype, axis=0) / (2 * genotype.shape[0])
        genotype = norm_binom(genotype, freq)
        genotypelist.append(genotype)
    return genotypelist
    
def pick_random_snps(studies, snpinfo, locus, cmax):
    nstudy = len(studies)
    snppool = list()
    for i in range(nstudy):
        snppool += [snp for snp in snpinfo[i][locus] if snp not in snppool]
    cval = random.randint(1, cmax)
    select = random.sample(snppool, cval)
    for snp in select:
        present_in = [studies[i] for i in range(nstudy) if snp in snpinfo[i][locus]]
        select[select.index(snp)] = snp._replace(studies = present_in)
    return select

def get_snppool(snpinfo):
    nstudy = len(snpinfo)
    snppool = list()
    nloci = len(snpinfo[0]) # Assumes same number of loci in each study
    for locus in range(nloci):
        snplocus = list()
        for i in range(nstudy):
            snplocus += [snp for snp in snpinfo[i][locus] if snp not in snplocus]
        snppool.append(snplocus)
    return snppool

def select_causal_snps(studies, snpinfo, snppool, locus, snp_score):
    select = list()
    prob = snp_score #/ np.sum(snp_score)
    none_selected = True
    while none_selected:
        for i, snp in enumerate(snppool):
            mrand = random.uniform(0, 1)
            if prob[i] > mrand:
                select.append(snp)
        if len(select) > 0:
            none_selected = False
        else:
            select = list()

    for snp in select:
        present_in = [study for i, study in enumerate(studies) if snp in snpinfo[i][locus]]
        select[select.index(snp)] = snp._replace(studies = present_in)
    return select


def select_causal_loci(nloci, nselect, featurelist, betadir, nfeat, nchoose, cfrac):
    choose = np.random.choice(nfeat, nchoose)
    betas = np.zeros(nfeat)
    beta0 = -np.log(-1 + (1 / cfrac))
    enrichment = np.random.uniform(low=2.0, high=8.0, size=nchoose)
    enrichbeta = - beta0 - np.log(((1 + np.exp(-beta0)) / enrichment) - 1)
    np.random.shuffle(enrichbeta)
    betas[choose] = enrichbeta
    betas[0] = beta0

    betafile = os.path.join(betadir, 'feature.betas')
    np.savetxt(betafile, betas)
    
    snp_prob = list()
    locus_prob = list()
    for locus in range(nloci):
        features = featurelist[locus]
        A = np.einsum('k, ki -> i', betas, features)
        pi = 1 / (1 + np.exp(-A))
        pi_locus = 1 - np.exp(np.sum(np.log(1 - pi)))
        locus_prob.append(pi_locus)
        snp_prob.append(pi)
    locus_prob = np.array(locus_prob)
    isort = np.argsort(locus_prob)[::-1]
    selected = isort[:nselect]
    return selected, snp_prob


opts = parse_args()
studies = opts.studies
samples = opts.samples

# proper files and paths
locidir = os.path.realpath(opts.locidir)
featdir = os.path.realpath(opts.featdir)
simdir  = os.path.realpath(opts.simdir)
sampledir = os.path.join(simdir, 'samples')
if not os.path.exists(sampledir):
    os.makedirs(sampledir)
betadir = os.path.join(sampledir, 'beta_values')
if not os.path.exists(betadir):
    os.makedirs(betadir)
effectfilename = 'snps.effectsize'
samplefilename = 'phenotypes.sample'

locusprefixes = read_locus_prefixes(opts.locusnames)
nloci = len(locusprefixes)

# Read all SNP info
snpinfo = read_snps(locidir, studies, locusprefixes)
snppool = get_snppool(snpinfo)
print ("Read all SNP info")

# Read all features
featurelist = read_features(locusprefixes, snppool, featdir)
print ("Read features")

# Select the causal loci
nfeat = featurelist[0].shape[0]
picks, snpscores = select_causal_loci(nloci, opts.ncausal_loci, featurelist, betadir, nfeat, 3, opts.proportion)
causality = np.zeros(nloci)
causality[picks] = 1
print ("Picked causal loci")

# Select causal SNPs per loci and write in causal.snplist
causal_snps = list()
ncausal = 0

fname = os.path.join(sampledir, opts.outfile)
with open(fname, 'w') as mfile:
    for i in range(nloci):
        mfile.write("#--------------------\n")
        mfile.write("Locus %i, Causality: %i \n" % (i+1, causality[i]))
        if causality[i] == 0:
            causal_snps.append(None)
        if causality[i] == 1:
            snps = select_causal_snps(studies, snpinfo, snppool[i], i, snpscores[i])
            causal_snps.append(snps)
            ncausal += len(snps)
            for snp in causal_snps[i]:
                mfile.write("%s %i \t %i %s\n" % (snp.rsid, snp.bp_location, len(snp.studies), snp.studies))

print ("Picked %d causal SNPs with maximum %d SNPs in single locus" % (ncausal, max([len(x) for x in causal_snps if x is not None])))


beta = np.zeros(ncausal)
if opts.simtype == 'fixed':
    beta = np.ones(ncausal)
elif opts.simtype == 'bimodal':
    mean1 = 0.5
    mean2 = -0.5
    mvar = 0.2
    for i in range(ncausal):
        mrandom = np.random.uniform()
        if mrandom <= 0.5:
            beta[i] = np.random.normal(mean1, mvar)
        else:
            beta[i] = np.random.normal(mean2, mvar)
elif opts.simtype == 'studentsT':
    beta = np.random.standard_t(1, size = ncausal)
elif opts.simtype == 'unimodal':
    beta = np.random.normal(0.5, 0.3)
else:
    beta = np.random.rand(ncausal)

beta *= np.sqrt( opts.sigma_herited_sq / np.sum(np.square(beta)) )

# Simulate phenotype for each study
for i, study in enumerate(studies):
    print ("Creating phenotype for cohort %s" % study)

    # Create folders
    outdir = os.path.join(sampledir, study)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Read genotype
    gendir = os.path.join(locidir, study)
    nsample = samples[i]
    target_cases = int(nsample / 2.1)
    gt = read_genotype(gendir, snpinfo[i], nsample)
    gt_tot = np.concatenate(gt, axis=1)
    allrsid = [snp.rsid for locus in snpinfo[i] for snp in locus]  # flatten out the snpinfo.rsid

    
    fname = os.path.join(outdir, effectfilename)
    mfile = open(fname, 'w')
    gtmask = np.zeros(gt_tot.shape[1], dtype=bool)
    betamask = np.zeros(beta.shape[0], dtype=bool)
    count = 0
    for locus in range(nloci):
        if causality[locus] == 1:
            for snp in causal_snps[locus]:
                if snp.rsid in allrsid:
                    k = allrsid.index(snp.rsid)
                    gtmask[k] = True
                    betamask[count] = True
                    mfile.write("%s \t %g\n" % (allrsid[k], beta[count]))
                else:
                    mfile.write("No SNPs found with rsid %s\n" % snp.rsid)
                count += 1
    mfile.close()

    gt_causal = gt_tot[:, gtmask]
    beta_sel = beta[betamask]

    # Get disease liabilities
    geno_eff = np.sum((beta_sel * gt_causal), axis=1)
    rand_var = np.sqrt(1 - opts.sigma_herited_sq)
    rand_eff = np.random.normal(0, rand_var, gt_causal.shape[0])
    liability = geno_eff + rand_eff

    # Select liabilities above the threshold of normal distribution truncating at K (disease prevalence)
    cutoff = scipy.stats.norm.ppf(1 - opts.prevalence)
    cases = np.where(liability >= cutoff)[0]
    #case_candidates = np.sum(liability >= cutoff)
    #print ("Case to be sampled from %i \n" % case_candidates)
    #if case_candidates < target_cases:
    #    raise ValueError("Target cases greater than cutoff, increase prevalence or decrease target cases")
    #else:
    #   case_candidates = np.where(liability >= cutoff)[0]
    #   cases = random.sample(list(case_candidates), target_cases)

    pheno = np.zeros(nsample, dtype=int)
    pheno[cases] = 1
    #for case in cases:
    #    pheno[case] = 1

    # Read sample file and write new
    counter = 0
    k = 0
    insample = os.path.join(gendir, '%s_QC.sample' % study)
    outsample = os.path.join(outdir, samplefilename)
    outfile = open(outsample, 'w')
    with open(insample, 'r') as samfile:
        for line in samfile:
            counter += 1
            if (counter <= 2):
               outfile.write(line)
            else:
               mline = line.split()
               newline = "%s %s %s %s %s %s %g\n" % (mline[0], mline[1], mline[2], mline[3], mline[4], mline[5], pheno[k])
               outfile.write(newline)
               k += 1

    print ("Artificial phenotypes generated.")
    print ("No. of cases:", np.sum(pheno))
    print ("No.of controls:", nsample - np.sum(pheno))
    print ()

# Create phenotype for each study

# Write the samples
