import numpy as np
import collections
import os
import matplotlib.pyplot as plt
from pip_comparison_functions import create_plot
from pip_comparison_functions import read_data

basedir='/mnt/storage/saikat/work/multivariate-gwas-bayesian-logistic-regression/quasi_laplace_gwas/results/mcmc_vs_qL'
mu_list = ['0', 'var']
bvsr_list = ['linear', 'probit']
phenotype_list = ['fixed', 'bimodal']
    
nloci = 200
nsim = 10
startsim = {'fixed': 0, 'bimodal': 100}

locusfile = os.path.join(basedir, 'LOCUSNAMES')
with open(locusfile, 'r') as mfile:
    locusprefixes = mfile.readlines()
locusprefixes = [x.strip() for x in locusprefixes]

for phenotype in phenotype_list:
    for bvsr in bvsr_list:
        for mu in mu_list:
            outdir = os.path.abspath('{:s}_effectsize'.format(phenotype))
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            filename = 'pip_{:s}_beta_{:s}_mcmc_qLmu{:s}.png'.format(phenotype, bvsr, mu)
            outfile = os.path.join(outdir, filename)
            xlabel = 'piMASS ({:s} BVSR model)'.format(bvsr)
            ylabel = 'BLORE (quasi-Laplace logistic)'
            print ("Creating {:s}".format(filename))
            result = read_data(nsim, startsim[phenotype], basedir, locusprefixes, mu, bvsr)
            create_plot(result, outfile, xlabel, ylabel)
