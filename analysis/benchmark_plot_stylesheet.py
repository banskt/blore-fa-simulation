import numpy as np
import collections
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def plot_benchmark(ax, data, color, myls, mlabel):
    xvals = data[0]
    yvals = data[1]
    ystd  = data[2]
    yupper = np.minimum(yvals + ystd, 1)
    ylower = yvals - ystd
    color = color
    linestyle = myls
    #ax.fill_between(xvals, ylower, yupper, color=color, alpha=0.15)

    coef = 4
    _dash = 1
    _dot = 0.3
    _dashspace = 0.8
    _dotspace = 0.5
    if myls == 'solid':
        mydash = []
    elif myls == 'dashed':
        mydash = [_dash, _dashspace]
    elif myls == 'dotted':
        mydash = [_dot, _dotspace]
    elif myls == 'dashdot':
        mydash = [_dash, _dashspace, _dot, _dashspace]
    elif myls == 'dashdotdot':
        mydash = [_dash, _dashspace, _dot, _dashspace, _dot, _dashspace]
    elif myls == 'dashdashdot':
        mydash = [_dash, _dashspace, _dashspace, _dashspace, _dot, _dashspace]
    mydash = [x * coef for x in mydash]

    ax.plot(xvals, yvals, color=color, dashes = mydash, lw=4, label=mlabel)
    return


def saveplot(filename, datadict, mlabels, mlim, legendloc, bboxpos, legendcol, locuspred=False, aucdict=None):
    ''' Use the same plot params for different benchmarks.
    '''

    kelly_colors_hex = [
        '#FFB300', # Vivid Yellow
        '#803E75', # Strong Purple
        '#FF6800', # Vivid Orange
        '#A6BDD7', # Very Light Blue
        '#C10020', # Vivid Red
        '#CEA262', # Grayish Yellow
        '#817066', # Medium Gray
    
        # The following don't work well for people with defective color vision
        '#007D34', # Vivid Green
        '#F6768E', # Strong Purplish Pink
        '#00538A', # Strong Blue
        '#FF7A5C', # Strong Yellowish Pink
        '#53377A', # Strong Violet
        '#FF8E00', # Vivid Orange Yellow
        '#B32851', # Strong Purplish Red
        '#F4C800', # Vivid Greenish Yellow
        '#7F180D', # Strong Reddish Brown
        '#93AA00', # Vivid Yellowish Green
        '#593315', # Deep Yellowish Brown
        '#F13A13', # Vivid Reddish Orange
        '#232C16', # Dark Olive Green
        ]
    
    colors = kelly_colors_hex
    marker_list = ['8', '>', 'd', '<', '*', 'p', '^', 's', 'h', 'v', 'D', r'$\clubsuit$']
    
    #bordercolor = '#2B2B2B'
    bordercolor = '#333333'
    borderwidth = 2
    figsize = (12,12)
    axis_font_size = 30
    label_font_size = 25
    legend_font_size = 25


    fig = plt.figure(figsize = figsize)
    ax1 = fig.add_subplot(111)

    mcolor = collections.defaultdict(lambda:0)
    mstyle = collections.defaultdict(lambda:0)
    mlegend = collections.defaultdict(lambda:0)

    mlegend['snptest']     = 'SNPTEST / META'
    mlegend['bimbam']      = 'BIMBAM'
    mlegend['metacca']     = 'metaCCA'
    mlegend['bammgwas']    = 'B-LORE'
    mlegend['bammgwas_fa'] = 'B-LORE FG'
    mlegend['caviarbf_c2'] = 'CAVIARBF'
    mlegend['finemap']     = 'FINEMAP'
    mlegend['paintor_comb']= 'PAINTOR'
    mlegend['paintor_fa']  = 'PAINTOR FG'

    mstyle['snptest']     = 'solid'
    mstyle['bimbam']      = 'dotted'
    mstyle['metacca']     = 'dashdot'
    mstyle['bammgwas']    = 'dashed'
    mstyle['bammgwas_fa'] = 'solid'
    mstyle['caviarbf_c2'] = 'dashdotdot'
    mstyle['finemap']     = 'dashdashdot'
    mstyle['paintor_comb']= 'dotted'
    mstyle['paintor_fa']  = 'solid'

    mcolor['snptest']     = colors[0]
    mcolor['bimbam']      = '#1f78b4' # colors[9] # blue
    mcolor['bammgwas']    = '#e31a1c' # colors[4] # red
    mcolor['bammgwas_fa'] = '#e31a1c' # colors[4]
    mcolor['metacca']     = colors[7] # green
    mcolor['finemap']     = '#33a02c' # colors[7] # colors[16]
    mcolor['caviarbf_c2'] = '#6a3d9a' # colors[18]
    mcolor['paintor_comb']= colors[6]
    mcolor['paintor_fa']  = colors[6]

    if (locuspred):

        #mlegend['snptest']     = 'SNPTEST / META (AUC = {:4.2f})'.format(aucdict['snptest'])
        #mlegend['bimbam']      = 'BIMBAM (AUC = {:4.2f})'.format(aucdict['bimbam'])
        #mlegend['metacca']     = 'metaCCA (AUC = {:4.2f})'.format(aucdict['metacca'])
        #mlegend['bammgwas']    = 'B-LORE (AUC = {:4.2f})'.format(aucdict['bammgwas'])
        #mlegend['bammgwas_fa'] = 'B-LORE FG (AUC = {:4.2f})'.format(aucdict['bammgwas_fa'])
        #mlegend['caviarbf_c2'] = 'CAVIARBF (AUC = {:4.2f})'.format(aucdict['caviarbf_c2'])
        #mlegend['finemap']     = 'FINEMAP (AUC = {:4.2f})'.format(aucdict['caviarbf_c2'])

        for key, val in datadict.items():
            if len(val[0]) > 0:
                plot_benchmark(ax1, val, mcolor[key], mstyle[key], mlegend[key])

    else:

        for key, val in datadict.items():
            if len(val[0]) > 0:
                plot_benchmark(ax1, val, mcolor[key], mstyle[key], mlegend[key])

    ax1.set_xlabel(mlabels[0], {'size': axis_font_size, 'color': bordercolor}, labelpad = 15)
    ax1.set_ylabel(mlabels[1], {'size': axis_font_size, 'color': bordercolor}, labelpad = 20)
    ax1.set_xlim(mlim[0], mlim[1])
    ax1.set_ylim(mlim[2], mlim[3])

    (lines, labels) = ax1.get_legend_handles_labels()
    legend = ax1.legend(lines, labels,
                        loc=legendloc,
                        bbox_to_anchor=bboxpos,
                        handlelength = 6.4,
                        frameon = True,
                        borderpad = 1.5,
                        labelspacing = 1.5,
                        ncol = legendcol)
    lframe = legend.get_frame()
    lframe.set_edgecolor(bordercolor)
    lframe.set_linewidth(borderwidth)
    lframe.set_facecolor('white')

    for fonts in ([legend.get_title()] + legend.texts):
        fonts.set_fontsize(legend_font_size)
        fonts.set_color(bordercolor)
    ax1.tick_params(axis='both', which = 'major',
                   length = 10, width = borderwidth, pad=10,
                   labelsize = label_font_size,
                   color = bordercolor,
                   labelcolor = bordercolor,
                   bottom = True, top = False, left = True, right = False
                  )
    for side, border in ax1.spines.items():
        border.set_linewidth(borderwidth)
        border.set_color(bordercolor)
    ax1.grid(color='dimgray', lw=0.5, alpha=0.5)

    plt.tight_layout()
    plt.savefig(filename, bbox_inches='tight')
    return
