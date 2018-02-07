import numpy as np
import benchmark_plot_stylesheet as mplt
import precisionld_recall as roc
import collections

nloci = 200
nsim = 50
startsim = 200
lddir = '/mnt/storage/saikat/work/multivariate-gwas-bayesian-logistic-regression/artificial_phenotype/metaprod/ldmap_nold/combined/'
#basedir = '/mnt/storage/saikat/work/multivariate-gwas-bayesian-logistic-regression/artificial_phenotype/metaprod/simulations_nold_hg0.25_vnormal/'
basedir='/mnt/storage/saikat/work/multivariate-gwas-bayesian-logistic-regression/artificial_phenotype/metaprod/simulations_nold_feat/'
recallplot = 'simu_finemap_recall_nold_30SNPs_tmp.pdf'
prcplot = 'simu_finemap_precision_recall_nold_tmp.pdf'
fdrplot = 'simu_finemap_altprecision_recall_nold_tmp.pdf'
mysimrange = range(startsim, nsim + startsim)
ldcut = 0.9

INFO_FIELDS = ['locus', 'rsid', 'stat', 'ld', 'causality']
class LDResult(collections.namedtuple('_LDResult', INFO_FIELDS)):
    __slots__ = ()

def xyplotvals_recall(data, ldcut, nloci):
    xvals = np.linspace(0, 200, 200)
    yests = list()
    for sim in range(nsim):
        _fpr, _tpr, _precision, _recall, _nselected, _precisionld, _fdr = roc.precisionld_recall_curve(data[sim], ldcut)
        _nselected = np.array(_nselected) / nloci
        yest = np.interp(xvals, _nselected, _recall)
        yests.append(yest)
    yests = np.array(yests)
    yvals = yests.mean(axis=0)
    ystd = yests.std(axis=0)
    return (xvals, yvals, ystd)


def xyplotvals_prc(data, ldcut, nloci):
    xvals = np.linspace(0, 1, 200)
    yests = list()
    for sim in range(nsim):
        _fpr, _tpr, _precision, _recall, _nselected, _precisionld, _fdr = roc.precisionld_recall_curve(data[sim], ldcut)
        yest = np.interp(xvals, _recall[1:], _precision[1:])
        yests.append(yest)
    yests = np.array(yests)
    yvals = yests.mean(axis=0)
    ystd = yests.std(axis=0)
    xvals = np.insert(xvals, 0, 0)
    yvals = np.insert(yvals, 0, 1)
    ystd = np.insert(ystd, 0, 0)
    return (xvals, yvals, ystd)

def xyplotvals_fdr(data, ldcut, nloci):
    xvals = np.linspace(0, 1, 200)
    yests = list()
    for sim in range(nsim):
        _fpr, _tpr, _precision, _recall, _nselected, _precisionld, _fdr = roc.precisionld_recall_curve(data[sim], ldcut)
        yest = np.interp(xvals, _recall[1:], _fdr[1:])
        yests.append(yest)
    yests = np.array(yests)
    yvals = yests.mean(axis=0)
    ystd = yests.std(axis=0)
    xvals = np.insert(xvals, 0, 0)
    yvals = np.insert(yvals, 0, 0)
    ystd = np.insert(ystd, 0, 0)
    return (xvals, yvals, ystd)



ldmatrix = [[] for x in range(nloci)]
for locus in range(nloci):
    rsidfile = lddir + 'Locus.%03i.LD.rsid' % (locus + 1)
    ldfile = lddir + 'Locus.%03i.LD' % (locus + 1)
    rsidlist = list()
    ldcol = collections.defaultdict(lambda:0)
    with open(rsidfile, 'r') as mfile:
        for mline in mfile:
            rsidlist.append(mline.split()[0])
    this_ldmat = np.loadtxt(ldfile)
    nsnps = len(rsidlist)
    for snp1 in range(nsnps):
        for snp2 in range(nsnps):
            rsid1 = rsidlist[snp1]
            rsid2 = rsidlist[snp2]
            ldcol[rsid1, rsid2] = np.square(this_ldmat[snp1, snp2])
    ldmatrix[locus] = ldcol


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
paintor_fa = list()
finemap = list()
snptest = list()
bimbam = list()

for sim in mysimrange:
    simdir = basedir + 'sim%03i/' % (sim + 1)

    causal_file = simdir + 'samples/causal.snplist'

    causal_rsids = [[] for x in range(nloci)]
    #causal_rsids = [[] for x in range(200)]
    counter = 0
    with open(causal_file, 'r') as mfile:
        for mline in mfile:
            if not mline.startswith('#'):
                mline_split =  mline.split()
                if mline_split[0] == 'Locus':
                    counter += 1
                else:
                    causal_rsids[counter - 1].append(mline_split[0])
    causal_loci = [min(len(x), 1) for x in causal_rsids]

    # BammGWAS
    thisres = list()
    for locus in range(nloci):
        outfile = simdir + 'blore/zmax4_muvar_pi0.01_sig0.01/blore_meta_res/Locus.%03i.gen.res' % (locus + 1)
        causals = causal_rsids[locus]
        ldrsq = ldmatrix[locus]
        with open(outfile, 'r') as mfile:
            for mline in mfile:
                mline_split = mline.split()
                if not mline_split[0] == 'Causal':
                    rsid = mline_split[0]
                    prob = float(mline_split[4])
                    #causality = 1 if rsid in causals else 0
                    ld = max([ldrsq[rsid, x] for x in causals]) if len(causals) > 0 else 0
                    causality = 1 if ld > ldcut else 0
                    mres = LDResult(locus = locus + 1,
                                    rsid = rsid,
                                    stat = prob,
                                    ld = ld,
                                    causality = causality)
                    thisres.append(mres)
    bammgwas.append(thisres)
    
    # BammGWAS with functional annotations
    thisres = list()
    for locus in range(nloci):
        outfile = simdir + 'blore/zmax4_muvar_pi0.01_sig0.01_feature/blore_meta_res/Locus.%03i.gen.res' % (locus + 1)
        causals = causal_rsids[locus]
        ldrsq = ldmatrix[locus]
        with open(outfile, 'r') as mfile:
            for mline in mfile:
                mline_split = mline.split()
                if not mline_split[0] == 'Causal':
                    rsid = mline_split[0]
                    prob = float(mline_split[4])
                    #causality = 1 if rsid in causals else 0
                    ld = max([ldrsq[rsid, x] for x in causals]) if len(causals) > 0 else 0
                    causality = 1 if ld > ldcut else 0
                    mres = LDResult(locus = locus + 1,
                                    rsid = rsid,
                                    stat = prob,
                                    ld = ld,
                                    causality = causality)
                    thisres.append(mres)
    bammgwas_feat.append(thisres)

    # CAVIARBF c1
    #thisres = list()
    #for locus in range(nloci):
    #    rsidfile = simdir + 'caviarbf/c1_weighted/Locus.%03i' % (locus + 1)
    #    resfile = simdir + 'caviarbf/c1_weighted/Locus.%03i.prior0.marginal' % (locus + 1)
    #    causals = causal_rsids[locus]
    #    ldrsq = ldmatrix[locus]
    #    rsidlist = list()
    #    with open(rsidfile, 'r') as mfile:
    #        for mline in mfile:
    #            rsidlist.append(mline.split()[0])
    #    with open(resfile, 'r') as mfile:
    #        for mline in mfile:
    #            mlinesplit = mline.split()
    #            index = int(mlinesplit[0]) - 1
    #            rsid = rsidlist[index]
    #            prob = float(mlinesplit[1])
    #            causality = 1 if rsid in causals else 0
    #            ld = max([ldrsq[rsid, x] for x in causals]) if len(causals) > 0 else 0
    #            mres = LDResult(locus = locus + 1,
    #                            rsid = rsid,
    #                            stat = prob,
    #                            ld = ld,
    #                            causality = causality)
    #            thisres.append(mres)
    #caviarbf_c1.append(thisres)
    
    # CAVIARBF c2
    thisres = list()
    for locus in range(nloci):
        rsidfile = simdir + 'caviarbf/c2_weighted/Locus.%03i' % (locus + 1)
        resfile = simdir + 'caviarbf/c2_weighted/Locus.%03i.prior0.marginal' % (locus + 1)
        causals = causal_rsids[locus]
        ldrsq = ldmatrix[locus]
        rsidlist = list()
        with open(rsidfile, 'r') as mfile:
            for mline in mfile:
                rsidlist.append(mline.split()[0])
        with open(resfile, 'r') as mfile:
            for mline in mfile:
                mlinesplit = mline.split()
                index = int(mlinesplit[0]) - 1
                rsid = rsidlist[index]
                prob = float(mlinesplit[1])
                #causality = 1 if rsid in causals else 0
                ld = max([ldrsq[rsid, x] for x in causals]) if len(causals) > 0 else 0
                causality = 1 if ld > ldcut else 0
                mres = LDResult(locus = locus + 1,
                                rsid = rsid,
                                stat = prob,
                                ld = ld,
                                causality = causality)
                thisres.append(mres)
    caviarbf_c2.append(thisres)

    # PAINTOR combined (= weighted) LD
    thisres = list()
    for locus in range(nloci):
        resfile = simdir + 'paintor/combined_LD/output/Locus.%03i.results' % (locus + 1)
        causals = causal_rsids[locus]
        ldrsq = ldmatrix[locus]
        with open(resfile, 'r') as mfile:
            next(mfile)
            for mline in mfile:
                mlinesplit = mline.split()
                rsid = mlinesplit[2]
                prob = float(mlinesplit[4])
                #causality = 1 if rsid in causals else 0
                ld = max([ldrsq[rsid, x] for x in causals]) if len(causals) > 0 else 0
                causality = 1 if ld > ldcut else 0
                mres = LDResult(locus = locus + 1,
                                rsid = rsid,
                                stat = prob,
                                ld = ld,
                                causality = causality)
                thisres.append(mres)
    paintor_comb.append(thisres)
    
    # PAINTOR individual LD
    #thisres = list()
    #for locus in range(nloci):
    #    resfile = simdir + 'paintor/individual_LD/output/Locus.%03i.results' % (locus + 1)
    #    causals = causal_rsids[locus]
    #    ldrsq = ldmatrix[locus]
    #    with open(resfile, 'r') as mfile:
    #        next(mfile)
    #        for mline in mfile:
    #            mlinesplit = mline.split()
    #            rsid = mlinesplit[1]
    #            prob = float(mlinesplit[7])
    #            causality = 1 if rsid in causals else 0
    #            ld = max([ldrsq[rsid, x] for x in causals]) if len(causals) > 0 else 0
    #            mres = LDResult(locus = locus + 1,
    #                            rsid = rsid,
    #                            stat = prob,
    #                            ld = ld,
    #                            causality = causality)
    #            thisres.append(mres)
    #paintor_indv.append(thisres)
    
    # PAINTOR combined (= weighted) LD + functional annotations
    thisres = list()
    for locus in range(nloci):
        resfile = simdir + 'paintor/combined_LD_FA/output_fa/Locus.%03i.results' % (locus + 1)
        causals = causal_rsids[locus]
        ldrsq = ldmatrix[locus]
        with open(resfile, 'r') as mfile:
            next(mfile)
            for mline in mfile:
                mlinesplit = mline.split()
                rsid = mlinesplit[2]
                prob = float(mlinesplit[4])
                #causality = 1 if rsid in causals else 0
                ld = max([ldrsq[rsid, x] for x in causals]) if len(causals) > 0 else 0
                causality = 1 if ld > ldcut else 0
                mres = LDResult(locus = locus + 1,
                                rsid = rsid,
                                stat = prob,
                                ld = ld,
                                causality = causality)
                thisres.append(mres)
    paintor_fa.append(thisres)
    
    # FINEMAP weighted LD
    thisres = list()
    for locus in range(nloci):
        resfile = simdir + 'finemap/c4/Locus.%03i.snp' % (locus + 1)
        causals = causal_rsids[locus]
        ldrsq = ldmatrix[locus]
        with open(resfile, 'r') as mfile:
            next(mfile)
            for mline in mfile:
                mlinesplit = mline.split()
                rsid = mlinesplit[1]
                prob = float(mlinesplit[3])
                #causality = 1 if rsid in causals else 0
                ld = max([ldrsq[rsid, x] for x in causals]) if len(causals) > 0 else 0
                causality = 1 if ld > ldcut else 0
                mres = LDResult(locus = locus + 1,
                                rsid = rsid,
                                stat = prob,
                                ld = ld,
                                causality = causality)
                thisres.append(mres)
    finemap.append(thisres)   
    
    # PVALUES
    thisres = list()
    for locus in range(nloci):
        resfile = simdir + 'snptest/meta/Locus.%03i.meta.out' % (locus + 1)
        causals = causal_rsids[locus]
        ldrsq = ldmatrix[locus]
        with open(resfile, 'r') as mfile:
            for mline in mfile:
                if not mline.startswith('#') and not mline.startswith('chr'):
                    mlinesplit = mline.split()
                    rsid = mlinesplit[1]
                    pval = float(mlinesplit[5])
                    #causality = 1 if rsid in causals else 0
                    ld = max([ldrsq[rsid, x] for x in causals]) if len(causals) > 0 else 0
                    causality = 1 if ld > ldcut else 0
                    mres = LDResult(locus = locus + 1,
                                    rsid = rsid,
                                    stat = -np.log(pval),
                                    ld = ld,
                                    causality = causality)
                    thisres.append(mres)                    
    snptest.append(thisres)


    # BIMBAM
    thisres = list()
    for locus in range(nloci):
        resfile = simdir + 'bimbam/meta/output/Locus.%03i.meta.ssd-bf.txt' % (locus + 1)
        nsnps = file_len(resfile) - 1
        causals = causal_rsids[locus]
        ldrsq = ldmatrix[locus]
        with open(resfile, 'r') as mfile:
            next(mfile)
            for mline in mfile:
                mlinesplit = mline.split()
                rsid = mlinesplit[0]
                log10bf = float(mlinesplit[1])
                pi = 2 / nsnps
                postodd = (np.power(10, log10bf) * pi) / (1 - pi)
                ppa = postodd / (1 + postodd)
                #causality = 1 if rsid in causals else 0
                ld = max([ldrsq[rsid, x] for x in causals]) if len(causals) > 0 else 0
                causality = 1 if ld > ldcut else 0
                mres = LDResult(locus = locus + 1,
                                rsid = rsid,
                                stat = ppa,
                                ld = ld,
                                causality = causality)
                thisres.append(mres)                   
    bimbam.append(thisres)

datadict = collections.defaultdict(lambda:0)
datadict['paintor_comb'] = paintor_comb 
datadict['paintor_fa'] = paintor_fa
datadict['finemap'] = finemap
#datadict['caviarbf_c1'] = caviarbf_c1
datadict['caviarbf_c2'] = caviarbf_c2
datadict['snptest'] = snptest
datadict['bimbam'] = bimbam
datadict['bammgwas'] = bammgwas
datadict['bammgwas_fa'] = bammgwas_feat

# Plot the recall
xydatadict = collections.defaultdict(lambda:0)
for key, val in datadict.items():
    if len(val) > 0:
        xydatadict[key] = xyplotvals_recall(val, ldcut, nloci)
mxlabel = r'Average number of SNPs selected per locus'
mylabel = r'Recall'
labels = (mxlabel, mylabel)
#limits = (0, 200, -0.02, 1.02)
limits = (0, 30, -0.02, 0.801)
legendloc = 'lower right'
#bboxpos = (0.95, 0.05)
bboxpos = (0.98, 0.02)
legendcol = 1
mplt.saveplot(recallplot, xydatadict, labels, limits, legendloc, bboxpos, legendcol)

# Plot the PR curve
xydatadict = collections.defaultdict(lambda:0)
for key, val in datadict.items():
    if len(val) > 0:
        xydatadict[key] = xyplotvals_prc(val, ldcut, nloci)
mxlabel = r'Recall'
mylabel = r'Precision'
labels = (mxlabel, mylabel)
limits = (-0.05, 1.05, -0.05, 1.05)
legendloc = 'upper right'
bboxpos = (0.95, 0.95)
legendcol = 1
mplt.saveplot(prcplot, xydatadict, labels, limits, legendloc, bboxpos, legendcol)

# Plot the "good" False Positive
xydatadict = collections.defaultdict(lambda:0)
for key, val in datadict.items():
    if len(val) > 0:
        xydatadict[key] = xyplotvals_fdr(val, ldcut, nloci)
mxlabel = r'Recall'
mylabel = r'Fraction of false positives with $r_{\mathrm{c}}^2$ > 0.9'
labels = (mxlabel, mylabel)
limits = (-0.05, 1.05, -0.05, 1.05)
legendloc = 'upper right'
bboxpos = (1.0065, 1.0075)
legendcol = 2
mplt.saveplot(fdrplot, xydatadict, labels, limits, legendloc, bboxpos, legendcol)
