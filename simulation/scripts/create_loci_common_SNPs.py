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


def copy_genotype(srcfile, genfile, mapfile, common_snps):
    common_rsids = [x.rsid for x in common_snps]
    dosage = list()
    with open(srcfile, 'r') as infile, open(genfile, 'w') as genout, open(mapfile, 'w') as mapout:
        for line in infile:
            linesplit = line.split()
            rsid = linesplit[1].strip()
            snpdosage = linesplit[5:]
            if rsid in common_rsids:
                snp = common_snps[common_rsids.index(rsid)]
                alt_allele = linesplit[3].strip()
                ref_allele = linesplit[4].strip()
                if alt_allele == snp.alt_allele and ref_allele == snp.ref_allele:
                    genout.write(line)
                elif alt_allele == snp.ref_allele and ref_allele == snp.alt_allele:
                    print ("WARNING: Flipping alleles")
                    snpdosage = [2 - x for x in snpdosage]
                    newline  = " ".join(linesplit[:5])
                    newline += " ".join(['{:g}'.format(x) for x in snpdosage])
                    newline += "\n"
                    genout.write(newline)
                else:
                    print ("WARNING: Alleles don't match for rsid {0}. Original {1} {2}. Found {3} {4}".format(rsid, snp.alt_allele, snp.ref_allele, alt_allele, ref_allele))
                    genout.write(line)
                mapout.write("{0} {1} 0 {2} {3} {4}\n".format(line.split()[0].strip(), snp.rsid, snp.bp_location, snp.alt_allele, snp.ref_allele))
                dosage.append(snpdosage)
    dosage = np.array(dosage)
    return dosage

def write_genotype(destdir, locusprefix, dosage, common_snps):
    if not os.path.exists(destdir):
        os.makedirs(destdir)
    genfile = os.path.join(destdir, locusprefix + '.gen')
    mapfile = os.path.join(destdir, locusprefix + '.map')
    with open(genfile, 'w') as genout, open(mapfile, 'w') as mapout:
        for i, snp in enumerate(common_snps):
            snpdosage = dosage[i]
            thisline  = "{0} {1} {2} {3} {4} ".format('---', snp.rsid, snp.bp_location, snp.alt_allele, snp.ref_allele)
            thisline += " ".join(['{:g}'.format(x) for x in snpdosage])
            thisline += "\n"
            genout.write(thisline)
            mapout.write("{0} {1} 0 {2} {3} {4}\n".format('---', snp.rsid, snp.bp_location, snp.alt_allele, snp.ref_allele))
    return


opts = parse_args()
studynames = opts.studies
locusprefix = opts.locusprefix
locidir = os.path.realpath(opts.locidir)
outdir = os.path.realpath(opts.outdir)

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

## Write common SNPs for each study
for i, study in enumerate(studynames):
    print (study)
    srcfile  = os.path.join(locidir, study, locusprefix + '.gen')
    study_dosage = read_genotype(srcfile, common_snps)
    destdir = os.path.join(outdir, study)
    write_genotype(destdir, locusprefix, study_dosage, common_snps)
    dosage.append(study_dosage)

combined_dosage = np.hstack(dosage)
destdir = os.path.join(outdir, 'combined')
write_genotype(destdir, locusprefix, combined_dosage, common_snps)
