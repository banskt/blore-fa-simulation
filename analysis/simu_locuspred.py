import numpy as np
import benchmark_plot_stylesheet as mplt
import receiving_operating_characteristic as roc
import collections
from sklearn import metrics

nloci = 200
nsim = 5
simstart = 1
basedir='/mnt/storage/saikat/work/multivariate-gwas-bayesian-logistic-regression/artificial_phenotype/metaprod/simulations_nold_feat/'
#basedir='/mnt/storage/saikat/work/multivariate-gwas-bayesian-logistic-regression/artificial_phenotype/metaprod/simulations_withld_feat/'
#basedir='/mnt/storage/saikat/work/multivariate-gwas-bayesian-logistic-regression/artificial_phenotype/metaprod/simulations_nold_hg0.25_vnormal/'
prcplot = 'simu_locus_precision_recall_nold.pdf'
#prcplot = 'simu_locus_precision_recall_withld.pdf'
mysimrange = range(simstart, nsim + simstart)

INFO_FIELDS = ['stat', 'causality']
class LocusResult(collections.namedtuple('_LocusResult', INFO_FIELDS)):
    __slots__ = ()

def xyplotvals(data):
    estimated = np.array([x.stat for x in data])
    actual = np.array([x.causality for x in data])
    _fpr, _tpr, _precision, _recall, _nselected, _auc = roc.roc_values(actual, estimated, -20000)
    xvals = np.array(_recall)
    yvals = np.array(_precision)
    ystd = np.zeros(xvals.shape)
    return (xvals, yvals, ystd)

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


bammgwas = list()
bammgwas_feat = list()
caviarbf_c1 = list()
caviarbf_c2 = list()
paintor_comb = list()
paintor_indv = list()
paintor_fa = list()
finemap = list()
snptest = list()
bimbam = list()
metacca = list()

for sim in mysimrange:
    simdir = basedir + 'sim%03i/' % (sim + 1)

    causal_file = simdir + 'samples/causal.snplist'

    causal_loci = list()
    with open(causal_file, 'r') as mfile:
        for mline in mfile:
            if mline.startswith('Locus'):
                mline_split = mline.split()
                causal_loci.append(int(mline_split[3]))

    # CAVIARBF c1
    #for locus in range(nloci):
    #    resfile = simdir + 'caviarbf/c1_weighted/Locus.%03i.prior0.statistics' % (locus + 1)
    #    with open(resfile, 'r') as mfile:
    #        next(mfile)
    #        next(mfile)
    #        mline = next(mfile).strip()
    #        mline_split = mline.split()
    #        prob = float(mline_split[0])
    #    causality = causal_loci[locus]
    #    mres = LocusResult(stat = prob, causality = causality)
    #    caviarbf_c1.append(mres)

    # CAVIARBF c2
    #for locus in range(nloci):
    #    resfile = simdir + 'caviarbf/c2_weighted/Locus.%03i.prior0.statistics' % (locus + 1)
    #    with open(resfile, 'r') as mfile:
    #        next(mfile)
    #        next(mfile)
    #        mline = next(mfile).strip()
    #        mline_split = mline.split()
    #        prob = float(mline_split[0])
    #    causality = causal_loci[locus]
    #    mres = LocusResult(stat = prob, causality = causality)
    #    caviarbf_c2.append(mres)

    ## PAINTOR combined LD
    #for locus in range(nloci):
    #    resfile = simdir + 'paintor/combined_LD/output/Locus.%03i.results' % (locus + 1)
    #    problist = list()
    #    with open(resfile, 'r') as mfile:
    #        next(mfile)
    #        for mline in mfile:
    #            mlinesplit = mline.split()
    #            rsid = mlinesplit[2]
    #            problist.append(float(mlinesplit[4]))
    #    prob = max(problist)
    #    causality = causal_loci[locus]
    #    mres = LocusResult(stat = prob, causality = causality)
    #    paintor_comb.append(mres)

    ## PAINTOR combined LD with annotation
    #for locus in range(nloci):
    #    resfile = simdir + 'paintor/combined_LD_FA/output_fa/Locus.%03i.results' % (locus + 1)
    #    problist = list()
    #    with open(resfile, 'r') as mfile:
    #        next(mfile)
    #        for mline in mfile:
    #            mlinesplit = mline.split()
    #            rsid = mlinesplit[2]
    #            problist.append(float(mlinesplit[4]))
    #    prob = max(problist)
    #    causality = causal_loci[locus]
    #    mres = LocusResult(stat = prob, causality = causality)
    #    paintor_fa.append(mres)
    #    
    ## PAINTOR individual LD
    #for locus in range(nloci):
    #    resfile = simdir + 'paintor/individual_LD/output/Locus.%03i.results' % (locus + 1)
    #    problist = list()
    #    with open(resfile, 'r') as mfile:
    #        next(mfile)
    #        for mline in mfile:
    #            mlinesplit = mline.split()
    #            rsid = mlinesplit[1]
    #            problist.append(float(mlinesplit[7]))
    #    prob = max(problist)
    #    causality = causal_loci[locus]
    #    mres = LocusResult(stat = prob, causality = causality)
    #    paintor_indv.append(mres)

            
    ## FINEMAP
    #for locus in range(nloci):
    #    resfile = simdir + 'finemap/c4/Locus.%03i.config' % (locus + 1)
    #    with open(resfile, 'r') as mfile:
    #        next(mfile)
    #        mline = mfile.readline()
    #        mlinesplit = mline.split()
    #        prob = float(mlinesplit[-2])
    #        causality = causal_loci[locus]
    #        mres = LocusResult(stat = prob, causality = causality)
    #        finemap.append(mres)
                
    # PVALUES
    for locus in range(nloci):
        resfile = simdir + 'snptest/meta/Locus.%03i.meta.out' % (locus + 1)
        problist = list()
        with open(resfile, 'r') as mfile:
            for mline in mfile:
                #if not mline.startswith('#') and not mline.startswith('alternate_ids'):
                if not mline.startswith('#') and not mline.startswith('chr'):
                    mlinesplit = mline.split()
                    rsid = mlinesplit[1]
                    pval = float(mlinesplit[5])
                    problist.append(-np.log10(pval))
        prob = max(problist)
        causality = causal_loci[locus]
        mres = LocusResult(stat = prob, causality = causality)
        snptest.append(mres)

    # BIMBAM
    for locus in range(nloci):
        resfile = simdir + 'bimbam/meta/output/Locus.%03i.meta.ssd-bf.txt' % (locus + 1)
        nsnps = file_len(resfile) - 1
        problist = list()
        with open(resfile, 'r') as mfile:
            next(mfile)
            for mline in mfile:
                mlinesplit = mline.split()
                rsid = mlinesplit[0]
                log10bf = float(mlinesplit[1])
                pi = 2 / nsnps
                postodd = (np.power(10, log10bf) * pi) / (1 - pi)
                ppa = postodd / (1 + postodd)
                problist.append(ppa)
        prob = max(problist)
        causality = causal_loci[locus]
        mres = LocusResult(stat = prob, causality = causality)
        bimbam.append(mres)

    # metaCCA 10 SNPs
    resfile = simdir + 'metaCCA/res_10.txt'
    with open(resfile, 'r') as mfile:
        for mline in mfile:
            mlinesplit = mline.split('\"')
            info = mlinesplit[1].split(' ')
            prob = float(info[1])
            locus = int(info[0].split('.')[1]) - 1
            causality = causal_loci[locus]
            mres = LocusResult(stat = prob, causality = causality)
            metacca.append(mres)

    # BammGWAS
    for locus in range(nloci):
        outfile = simdir + 'blore/zmax4_mu0_pi0.01_sig0.01/blore_meta_res/Locus.%03i.gen.res' % (locus + 1)
        with open(outfile, 'r') as mfile:
            mline = next(mfile).strip()
            mline_split = mline.split()
            prob = float(mline_split[2])
        causality = causal_loci[locus]
        mres = LocusResult(stat = prob, causality = causality)
        bammgwas.append(mres)
        
    # BammGWAS with functional annotations
    for locus in range(nloci):
        outfile = simdir + 'blore/zmax4_mu0_pi0.01_sig0.01_feature/blore_meta_res/Locus.%03i.gen.res' % (locus + 1)
        with open(outfile, 'r') as mfile:
            mline = next(mfile).strip()
            mline_split = mline.split()
            prob = float(mline_split[2])
            #if prob < 0.1:
            #    prob = 0.1 * random.uniform(0, 1)
        causality = causal_loci[locus]
        mres = LocusResult(stat = prob, causality = causality)
        bammgwas_feat.append(mres)
        

datadict = collections.defaultdict(lambda:0)
datadict['snptest'] = snptest
datadict['bimbam'] = bimbam
datadict['metacca'] = metacca
datadict['bammgwas'] = bammgwas
datadict['bammgwas_fa'] = bammgwas_feat
datadict['paintor_comb'] = paintor_comb 
datadict['paintor_indv'] = paintor_indv 
datadict['paintor_fa'] = paintor_fa
datadict['caviarbf_c1'] = caviarbf_c1
datadict['caviarbf_c2'] = caviarbf_c2
datadict['finemap'] = finemap

# Plot the PR curve
xydatadict = collections.defaultdict(lambda:0)
aucdict = collections.defaultdict(lambda:0)
for key, val in datadict.items():
    if len(val) > 0:
        xydatadict[key] = xyplotvals(val)
        estimated = np.array([x.stat for x in val])
        actual = np.array([x.causality for x in val])
        aucdict[key] = metrics.roc_auc_score(actual, estimated)
        #auc = metrics.roc_auc_score(actual, estimated)
        #print ("{:s}:\t{:g}".format(key, auc))
        
mxlabel = r'Recall'
mylabel = r'Precision'
labels = (mxlabel, mylabel)
limits = (-0.05, 1.05, 0.48, 1.02)
legendloc = 'lower left'
#legendloc = 'lower center'
bboxpos = (0.05, 0.05)
#bboxpos = (0.5, 0.05)
legendcol = 1
#xydatadict['metacca'][1][2] = 0.8
mplt.saveplot(prcplot, xydatadict, labels, limits, legendloc, bboxpos, legendcol, locuspred=True, aucdict = aucdict)
