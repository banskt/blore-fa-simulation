import numpy as np
import collections
import os
import matplotlib.pyplot as plt
import pip_comparison_functions

INFO_FIELDS = ['locus', 'rsid', 'qlpip', 'mcmcpip' , 'ld', 'causality']
class PIPResult(collections.namedtuple('_PIPResult', INFO_FIELDS)):
    __slots__ = ()

def create_plot(result, outfile, mxlabel, mylabel):
    blore_true = [np.log10(x.qlpip)   for x in result if x.causality == 1]
    bvsr_true  = [np.log10(x.mcmcpip) for x in result if x.causality == 1]
    blore_false = [np.log10(x.qlpip)   for x in result if x.causality == 0]
    bvsr_false  = [np.log10(x.mcmcpip) for x in result if x.causality == 0]
    
    bordercolor = '#2B2B2B'
    bordercolor = '#333333'
    borderwidth = 2
    colors = ['#A6BDD7', '#C10020']
    figsize = (12, 12)
    axis_font_size = 30
    label_font_size = 25
    legend_font_size = 25
    
    fig = plt.figure(figsize = figsize)
    ax1 = fig.add_subplot(111)
    
    mlabel = "Non-causal SNPs"
    ax1.scatter(bvsr_false, blore_false, color=colors[0], s = 10, alpha = 0.3, label = mlabel)
    ax1.plot([-6, 1], [-6, 1], lw=borderwidth, color=bordercolor)
    mlabel = "Causal SNPs"
    ax1.scatter(bvsr_true, blore_true, color=colors[1], s = 10, alpha = 1, label = mlabel)
    
    ax1.set_xlabel(mxlabel, {'size': axis_font_size, 'color': bordercolor}, labelpad = 60)
    ax1.set_ylabel(mylabel, {'size': axis_font_size, 'color': bordercolor}, labelpad = 60)
    
    for ax in [ax1]:
        ax.set_ylim(-2.8, 0.2)
        ax.set_xlim(-2.8, 0.2)
        mticks = [-2, -1, 0]
        ax.set_xticks(mticks)
        ax.set_yticks(mticks)
        tl = ["{:g}".format(np.power(10.,x)) for x in mticks]
        ax.set_xticklabels(tl)
        ax.set_yticklabels(tl)
    
        legendtitle = 'SNPs'
        legend = ax.legend(loc='upper left', bbox_to_anchor=(0.02, 0.98),
                           #scatterpoints = 100,
                           markerscale=5,
                           frameon = True, borderpad = 1.5, labelspacing = 1.5
                           #title = legendtitle
                          )
        for l in legend.legendHandles:
            l.set_alpha(1)
        lframe = legend.get_frame()
        lframe.set_edgecolor(bordercolor)
        lframe.set_linewidth(borderwidth)
    
        ax.tick_params(axis='both', which = 'major',
                       length = 10, width = borderwidth, pad=10,
                       labelsize = label_font_size,
                       color = bordercolor,
                       labelcolor = bordercolor,
                       bottom = True, top = False, left = True, right = False
                      )
        for side, border in ax.spines.items():
            border.set_linewidth(borderwidth)
            border.set_color(bordercolor)
        for fonts in ([legend.get_title()] + legend.texts):
            fonts.set_fontsize(legend_font_size)
            fonts.set_color(bordercolor)
        ax.grid(color='dimgray', lw=0.5, alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(outfile, bbox_inches='tight')
    return


def read_data(nsim, startsim, basedir, locusprefixes, muvar, bvsr_model):
    result = list()
    for sim in range(startsim, startsim + nsim):
        simdir = os.path.join( basedir, 'simulations', 'sim{:03d}'.format(sim + 1) )
        causal_file = os.path.join( simdir, 'samples/causal.snplist' )
        causal_rsids = list()
        with open(causal_file, 'r') as mfile:
            for mline in mfile:
                if not mline.startswith('#'):
                    mline_split =  mline.split()
                    if mline_split[0] != 'Locus':
                        causal_rsids.append(mline_split[0])
        for locus in locusprefixes:
            rsid_list = list()
            ql_pip = collections.defaultdict(lambda:0)
            mcmc_pip = collections.defaultdict(lambda:0)
            path = 'blore/meta_without_feature/zmax4_mu{:s}_pi0.01_sig0.01/blore_meta_res/{:s}.gen.res'
            outfile = os.path.join( simdir, path.format(muvar, locus) )
            with open(outfile, 'r') as mfile:
                for mline in mfile:
                    mline_split = mline.split()
                    if not mline_split[0] == 'Causal':
                        rsid = mline_split[0].strip()
                        ql_pip[rsid] = max(float(mline_split[4]), 1e-8)
                        if rsid not in rsid_list:
                            rsid_list.append(rsid)
            path = 'pimass/c4_5e6_{:s}/output/{:s}.mcmc.txt'
            outfile = os.path.join( simdir, path.format(bvsr_model, locus) )
            with open(outfile, 'r') as mfile:
                next(mfile)
                for mline in mfile:
                    mline_split = mline.split()
                    rsid = mline_split[0].strip()
                    mcmc_pip[rsid] = max(float(mline_split[3]), 1e-8)
                    if rsid not in rsid_list:
                        rsid_list.append(rsid)
            for rsid in rsid_list:
                causality = 1 if rsid in causal_rsids else 0
                ld = 0 # max([ldrsq[rsid, x] for x in causals]) if len(causals) > 0 else 0
                mres = PIPResult(locus = locus,
                                 rsid = rsid,
                                 qlpip = ql_pip[rsid],
                                 mcmcpip = mcmc_pip[rsid],
                                 ld = 0,
                                 causality = causality)
                result.append(mres)
    return result
