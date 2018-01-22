'''
Read the SNPs from all studies, find the common SNPs and write them for each study + combined

Requires the following files for each 'study':
    locidir/study/locusprefix.map
    locidir/study/locusprefix.gen

Outputs the following files:

 1) for each 'study':
      outdir/study/locusprefix.map [SNP annotation]
      outdir/study/locusprefix.gen [Oxford]
      outdir/study/locusprefix.bgen [Oxford]
      outdir/study/locusprefix.matgen [BIMBAM]
      outdir/study/locusprefix.bimbam_map [BIMBAM SNP annotation]
 2) for combined:
      outdir/combined/locusprefix.map [SNP annotation]
      outdir/combined/locusprefix.gen [Oxford]
      outdir/combined/locusprefix.bgen [Oxford]
      outdir/combined/locusprefix.matgen [BIMBAM]
      outdir/combined/locusprefix.bimbam_map [BIMBAM SNP annotation]
'''

import argparse
import os
import collections
import subprocess
import numpy as np

def parse_args():

    parser = argparse.ArgumentParser(description='Pick random SNPs for simulation')

    parser.add_argument('-dl', '--locidir',
                        type=str,
                        dest='locidir',
                        metavar='DIR',
                        help='path of the loci dosages directory')

    parser.add_argument('-fl', '--locusprefix',
                        type=str,
                        dest='locusprefix',
                        metavar='FILE PREFIX',
                        help='prefix of locusname')

    parser.add_argument('-st', '--studies',
                        nargs='*',
                        type=str,
                        dest='studies',
                        metavar='STR',
                        help='list of study names')

    parser.add_argument('-do', '--outdir',
                        type=str,
                        dest='outdir',
                        metavar='STR',
                        help='path of the output directory')

    parser.add_argument('--qctool', 
                        type=str,
                        dest='qctool',
                        metavar='STR',
                        help='path of the qctool')

    opts = parser.parse_args()
    return opts


SNPINFO_FIELDS = ['rsid', 'bp_location', 'alt_allele', 'ref_allele', 'studies']
class SnpInfo(collections.namedtuple('_SnpInfo', SNPINFO_FIELDS)):
    __slots__ = ()

def read_snps(locidir, studies, locusprefix):
    nstudy = len(studies)
    snpinfo = list()
    for study in studies:
        snpstudy = list()
        mapfile = os.path.join(locidir, study, locusprefix + '.map')
        with open(mapfile, 'r') as mfile:
            for line in mfile:
                mline = line.split()
                this_snp = SnpInfo(rsid = mline[1],
                                   bp_location = int(mline[3]),
                                   alt_allele = mline[4],
                                   ref_allele = mline[5],
                                   studies = [])
                snpstudy.append(this_snp)
        snpinfo.append(snpstudy)
    return snpinfo


def read_genotype(srcfile, common_snps):
    dosage = [None for x in common_snps]
    common_rsids = [x.rsid for x in common_snps]
    with open(srcfile, 'r') as infile:
        for line in infile:
            linesplit = line.split()
            rsid = linesplit[1].strip()
            if rsid in common_rsids:
                mindex = common_rsids.index(rsid)
                snp = common_snps[mindex]
                alt_allele = linesplit[3].strip()
                ref_allele = linesplit[4].strip()
                snpdosage = linesplit[5:]
                if alt_allele == snp.alt_allele and ref_allele == snp.ref_allele:
                    dosage[mindex] = snpdosage
                elif alt_allele == snp.ref_allele and ref_allele == snp.alt_allele:
                    print ("WARNING: Flipping alleles")
                    snpdosage = [2 - x for x in snpdosage]
                    dosage[mindex] = snpdosage
                else:
                    print ("WARNING: Alleles don't match for {0}. Original {1} {2}. Found {3} {4}".format(rsid, snp.alt_allele, snp.ref_allele, alt_allele, ref_allele))
                    dosage[mindex] = snpdosage

    dosage = np.array(dosage, dtype=np.float64)
    return dosage


def create_gt_matrix(dosage):
    nsnps = dosage.shape[0]
    nsample = int(dosage.shape[1] / 3)
    genotype = np.empty([nsample, nsnps], dtype=np.float64)
    for j in range(nsample):
        maj_ind = j * 3            # Allele AA, where A is the alternate allele
        het_ind = maj_ind + 1      # Allele AB, where B is the reference allele
        min_ind = het_ind + 1      # Allele BB
        genotype[j,:] = 2 * dosage[:, min_ind] + dosage[:, het_ind] # [AA, AB, BB] := [0, 1, 2]
    return genotype.T


def write_genotype(destdir, locusprefix, dosage, common_snps, qctool):
    if not os.path.exists(destdir):
        os.makedirs(destdir)
    genfile = os.path.join(destdir, locusprefix + '.gen')
    mapfile = os.path.join(destdir, locusprefix + '.map')
    bgenfile = os.path.join(destdir, locusprefix + '.bgen')
    with open(genfile, 'w') as genout, open(mapfile, 'w') as mapout:
        for i, snp in enumerate(common_snps):
            snpdosage = dosage[i]
            thisline  = "{0} {1} {2} {3} {4} ".format('---', snp.rsid, snp.bp_location, snp.alt_allele, snp.ref_allele)
            thisline += " ".join(['{:g}'.format(x) for x in snpdosage])
            thisline += "\n"
            genout.write(thisline)
            mapout.write("{0} {1} 0 {2} {3} {4}\n".format('---', snp.rsid, snp.bp_location, snp.alt_allele, snp.ref_allele))
    qctool_command = '{:s} -g {:s} -og {:s}'.format(qctool, genfile, bgenfile)
    subprocess.call([qctool_command], shell=True)
    return

def write_matgen(destdir, locusprefix, genmat, common_snps):
    if not os.path.exists(destdir):
        os.makedirs(destdir)
    genfile = os.path.join(destdir, locusprefix + '.matgen')
    mapfile = os.path.join(destdir, locusprefix + '.matmap')
    with open(genfile, 'w') as genout, open(mapfile, 'w') as mapout:
        for i, snp in enumerate(common_snps):
            snpdosage = genmat[i]
            thisline  = "{0}, {1}, {2}, ".format(snp.rsid, snp.ref_allele, snp.alt_allele)
            thisline += ", ".join(['{:g}'.format(x) for x in snpdosage])
            thisline += "\n"
            genout.write(thisline)
            mapout.write("{0}, {1}, 1\n".format(snp.rsid, snp.bp_location, 1))
    return



opts = parse_args()
studynames = opts.studies
locusprefix = opts.locusprefix
locidir = os.path.realpath(opts.locidir)
outdir = os.path.realpath(opts.outdir)
qctool = opts.qctool

nstudy = len(studynames)

## Read the SNPs in all studies
snpinfo = read_snps(locidir, studynames, locusprefix)

## Find the common SNPs
rsid_sets = list()
for i in range(nstudy):
    study_set = set([x.rsid for x in snpinfo[i]])
    rsid_sets.append(study_set)

common_rsids = rsid_sets[0]
for i in range(1, nstudy):
    updated_common_rsids = common_rsids & rsid_sets[i]
    common_rsids = updated_common_rsids

common_snps = [x for x in snpinfo[-1] if x.rsid in common_rsids] # this works because common SNPs must be present in the last study as well.

dosage = list()
gtmatrix = list()

## Write common SNPs for each study
for i, study in enumerate(studynames):
    print (study)
    srcfile  = os.path.join(locidir, study, locusprefix + '.gen')
    study_dosage = read_genotype(srcfile, common_snps)
    study_gtmatrix = create_gt_matrix(study_dosage)
    destdir = os.path.join(outdir, study)
    write_genotype(destdir, locusprefix, study_dosage, common_snps, qctool)
    write_matgen(destdir, locusprefix, study_gtmatrix, common_snps)
    dosage.append(study_dosage)
    gtmatrix.append(study_gtmatrix)
    

combined_dosage = np.hstack(dosage)
combined_gtmatrix = np.hstack(gtmatrix)
destdir = os.path.join(outdir, 'combined')
write_genotype(destdir, locusprefix, combined_dosage, common_snps, qctool)
write_matgen(destdir, locusprefix, combined_gtmatrix, common_snps)
