 #!/usr/bin/env python

'''pmrm vs ANI plot.

This script plots pm/rm vs ANI from the pairwise all vs. all pmrm data.

It builds a figure similar to 01b_fastANI_scatter_pyGAM.py but with the
difference that pmrm is plotted on the y-axis against ANI.

This script requires the following input parameters:

    * the *_pmpr_pairwise_data.tsv file from 04b_pmrrm_analyses.py
    * Output file prefix

This script returns two publication ready figures in pdf format.

*NOTE change -ymin to negative values when using -lg True.

This script requires python 3.6+ and the following modules:

    * matplotlib
    * numpy
    * pandas
    * seaborn
    * scipy
    * pygam - https://pygam.readthedocs.io/
        - pip install pygam
        - conda install -c conda-forge pygam
        - conda install -c conda-forge scikit-sparse nose
        - goes way faster with more RAM. More RAM, more fast.
    * datashader for -z True density plot
        - pip install datashader
        - conda install datashader
        ** point size and alpha do not work with -z True

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Jan 2024
License :: GNU GPLv3
Copyright 2024 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse, sys, random
from collections import defaultdict
import numpy as np, pandas as pd; np.random.seed(0)
from scipy.stats import spearmanr as corr
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns; sns.set(style="white", color_codes=True)
import datashader as ds
from datashader.mpl_ext import dsshow
from pygam import LinearGAM, s, l

class TickRedrawer(matplotlib.artist.Artist):
    #https://stackoverflow.com/questions/19677963/
    #matplotlib-keep-grid-lines-behind-the-graph-but-the-y-and-x-axis-above
    """Artist to redraw ticks."""

    __name__ = "ticks"

    zorder = 10

    @matplotlib.artist.allow_rasterization
    def draw(self, renderer: matplotlib.backend_bases.RendererBase) -> None:
        """Draw the ticks."""
        if not self.get_visible():
            self.stale = False
            return

        renderer.open_group(self.__name__, gid=self.get_gid())

        for axis in (self.axes.xaxis, self.axes.yaxis):
            loc_min, loc_max = axis.get_view_interval()

            for tick in axis.get_major_ticks() + axis.get_minor_ticks():
                if tick.get_visible() and loc_min <= tick.get_loc() <= loc_max:
                    for artist in (tick.tick1line, tick.tick2line):
                        artist.draw(renderer)

        renderer.close_group(self.__name__)
        self.stale = False


def parse_pmrm_file(pmrm_file, lg):
    """Read the file and grabs the data"""
    # initialize dict to store data
    data = {'gpair': [], 'pmrm': [], 'ani': [], 'f100': []}
    # read through the file and populate the dict
    with open(pmrm_file, 'r') as file:
        # grab the header
        header = file.readline()
        for line in file:
            # split each line by tabs
            X = line.rstrip().split('\t')
            # get the genome pair names, pmrm, ANI, and F100
            gpair, pmrm, ani, f100 = X[0], float(X[1]), float(X[2]), float(X[3])
            if pmrm == 0: continue
            data['gpair'].append(gpair)
            pmrm = np.log(pmrm) if lg else pmrm
            data['pmrm'].append(pmrm)
            data['ani'].append(ani)
            data['f100'].append(f100)

    # convert data dict to pandas data frame
    df = pd.DataFrame(data)

    return df


def gather_stats(xs, ys):
    """Computes correlation, mean, and median on df columns xs and ys """

    # Compute Pearson Correlation Coefficient
    print("\nCalculating statistics.")
    pcorr = corr(xs, ys)

    # Compute ANI mean and median
    ani_mean = np.mean(xs)
    ani_median = np.median(xs)
    frag_mean = np.mean(ys)
    frag_median = np.median(ys)

    # Compile dictionairy
    df_stats = {
        'pcorr': pcorr,
        'ani_mean': ani_mean,
        'ani_median': ani_median,
        'frag_mean': frag_mean,
        'frag_median': frag_median
        }

    return df_stats


def get_GAM(df):

    # Trendline with pyGAM
    print('\nCalculating trendline with pyGAM.')
    X = df["ani"].to_numpy()
    X = X[:, np.newaxis]
    y = df["pmrm"].to_list()

    # Build and fit the GAM model using gridsearch to tune param
    lam = np.logspace(-5, 5, 10)
    # A
    splineA = s(0, n_splines=4, constraints='convex')
    gamA = LinearGAM(splineA).gridsearch(X, y, lam=lam)
    r2A = gamA.statistics_['pseudo_r2']['explained_deviance']
    # B
    splineB = s(0, n_splines=4, constraints='concave')
    gamB = LinearGAM(splineB).gridsearch(X, y, lam=lam)
    r2B = gamB.statistics_['pseudo_r2']['explained_deviance']
    # C
    splineC = s(0, n_splines=4, constraints='monotonic_inc')
    gamC = LinearGAM(splineC).gridsearch(X, y, lam=lam)
    r2C = gamC.statistics_['pseudo_r2']['explained_deviance']

    # AA
    splineAA = s(0, n_splines=6, constraints='convex')
    gamAA = LinearGAM(splineAA).gridsearch(X, y, lam=lam)
    r2AA = gamAA.statistics_['pseudo_r2']['explained_deviance']
    # BB
    splineBB = s(0, n_splines=6, constraints='concave')
    gamBB = LinearGAM(splineBB).gridsearch(X, y, lam=lam)
    r2BB = gamBB.statistics_['pseudo_r2']['explained_deviance']
    # CC
    splineCC = s(0, n_splines=6, constraints='monotonic_inc')
    gamCC = LinearGAM(splineCC).gridsearch(X, y, lam=lam)
    r2CC = gamCC.statistics_['pseudo_r2']['explained_deviance']

    # AA
    splineAAA = s(0, n_splines=8, constraints='convex')
    gamAAA = LinearGAM(splineAAA).gridsearch(X, y, lam=lam)
    r2AAA = gamAAA.statistics_['pseudo_r2']['explained_deviance']
    # BB
    splineBBB = s(0, n_splines=8, constraints='concave')
    gamBBB = LinearGAM(splineBB).gridsearch(X, y, lam=lam)
    r2BBB = gamBBB.statistics_['pseudo_r2']['explained_deviance']
    # CC
    splineCCC = s(0, n_splines=8, constraints='monotonic_inc')
    gamCCC = LinearGAM(splineCCC).gridsearch(X, y, lam=lam)
    r2CCC = gamCCC.statistics_['pseudo_r2']['explained_deviance']

    # D
    gamD = LinearGAM(l(0)).gridsearch(X, y)
    r2D = gamD.statistics_['pseudo_r2']['explained_deviance']

    # choose the best
    gams = [gamA, gamB, gamC, gamAA, gamBB, gamCC, gamAAA, gamBBB, gamCCC, gamD]
    r2s = [r2A, r2B, r2C, r2AA, r2BB, r2CC, r2AAA, r2BBB, r2CCC, r2D]

    r2 = max(r2s)
    i = r2s.index(r2)
    gam = gams[i]
    XX = gam.generate_X_grid(term=0, n=500)

    return gam, XX


def pmrm_ani_plot(df, xmin, xmax, xt, ymin, ymax, yt, z, p, a, gam, mm, lg, op, s):

    # Gather Stats
    df_stats = gather_stats(df['ani'], df['pmrm'])

    stats_line = (
        f"Spearman rho: {round(df_stats['pcorr'][0], 2)}\n"
        f"p value: {round(df_stats['pcorr'][1], 2)}"
        )

    # Set Colors and markers
    main_color = '#933b41'
    second_color = '#737373'
    vline_color = '#000000'
    color = '#252525'
    marker = '.' #'o'

    # build plot
    g = sns.JointGrid(x="ani", y="pmrm", data=df)
    print('\nComputing KDEs for marginal plots.')
    # x margin kde plot
    sns.kdeplot(
            x=df["ani"],
            ax=g.ax_marg_x,
            legend=False,
            color=color
            )
    # y margin kde plot
    sns.kdeplot(
            y=df['pmrm'],
            ax=g.ax_marg_y,
            legend=False,
            color=color
            )
    # main panel scatter plot
    if z: # density scatter plot with datashader
        print('\nComputing plot densities.')
        import datashader as ds
        from datashader.mpl_ext import dsshow
        dsartist = dsshow(
                        df,
                        ds.Point("ani", "pmrm"),
                        ds.count(),
                        norm="log",
                        aspect="auto",
                        ax=g.ax_joint,
                        width_scale=3.,
                        height_scale=3.
                        )

    else: # regular scatter plot
        print('\nPlotting data.')
        n = len(df)
        rast = True if n >= 5000 else False
        g.ax_joint.plot(
                df["ani"],
                df["pmrm"],
                marker,
                ms=p,
                alpha=a,
                color=color,
                rasterized=rast,
                )

    if gam:
        # Trendline with pyGAM

        gam, XX = get_GAM(df)

        g.ax_joint.plot(
                    XX,
                    gam.predict(XX),
                    color='#2171b5',
                    linestyle='--',
                    linewidth=1.0,
                    )
        g.ax_joint.plot(
                    XX,
                    gam.prediction_intervals(XX, width=0.95),
                    color='#9ecae1',
                    linestyle='--',
                    linewidth=1.0
                    )
        r2 = gam.statistics_['pseudo_r2']['explained_deviance']
        GAM_line = f"GAM Pseudo R-Squared: {r2:.4f}"
        g.ax_joint.text(
            0.75, 0.1, GAM_line,
            fontsize=10, color=second_color,
            verticalalignment='top', horizontalalignment='right',
            transform=g.ax_joint.transAxes
            )

    # plot title, labels, text
    ptitle = f'r/m vs. ANI'
    g.ax_marg_x.set_title(ptitle, fontsize=18, y=1.02)

    g.ax_marg_x.text(
        -0.05, 1.02, s, style='italic',
        fontsize=18, #color=second_color,
        verticalalignment='bottom', horizontalalignment='left',
        transform=g.ax_marg_x.transAxes
        )

    g.ax_joint.set_xlabel(
        'Average nucleotide identity (%)',
        fontsize=12, y=-0.02
        )
    ylab = 'ln(r/m)' if lg else 'r/m'
    g.ax_joint.set_ylabel(
        ylab,
        fontsize=12, x=-0.02
        )
    g.ax_joint.text(
        0.28, 0.98, stats_line,
        fontsize=10, color=second_color,
        verticalalignment='top', horizontalalignment='right',
        transform=g.ax_joint.transAxes
        )

    # set the axis parameters / style
    hstep = xt/10
    g.ax_joint.set_xticks(np.arange(xmin, xmax+hstep, xt))
    g.ax_joint.set_xlim(left=xmin-hstep, right=xmax+hstep)


    istep = yt/10
    g.ax_joint.set_yticks(np.arange(ymin, ymax+istep, yt))
    g.ax_joint.set_ylim(bottom=ymin-istep, top=ymax+istep)
    g.ax_joint.tick_params(axis='both', labelsize=12)
    g.ax_joint.tick_params(
        axis='both', which='major', direction='inout', color='k',
        width=2, length=6, bottom=True, left=True, zorder=3
        )

    # set grid style
    g.ax_joint.yaxis.grid(
        which="major", color='#d9d9d9', linestyle='--', linewidth=1
        )
    g.ax_joint.xaxis.grid(
        which="major", color='#d9d9d9', linestyle='--', linewidth=1
        )
    g.ax_joint.set_axisbelow(True)
    g.ax_joint.add_artist(TickRedrawer())

    if mm:
        # Plot mean and median
        _ = g.ax_joint.axvline(
            x=df_stats['ani_mean'], ymin=0, ymax=1,
            color=vline_color, linewidth=2, linestyle='--',
            label='Mean'
            )
        _ = g.ax_joint.axhline(
            y=df_stats['frag_mean'], xmin=0, xmax=1,
            color=vline_color, linewidth=2, linestyle='--',
            )
        _ = g.ax_joint.axvline(
            x=df_stats['ani_median'], ymin=0, ymax=1,
            color=vline_color, linewidth=2, linestyle=':',
            label='Mean'
            )
        _ = g.ax_joint.axhline(
            y=df_stats['frag_median'], xmin=0, xmax=1,
            color=vline_color, linewidth=2, linestyle=':',
            )

        # Build legend for mean and median
        g.ax_joint.legend(
            loc='lower left',
            fontsize=12,
            markerscale=1.5,
            numpoints=1,
            frameon=False,
            ncol=2
            )

    # adjust layout, save, and close
    g.fig.set_figwidth(7)
    g.fig.set_figheight(5)
    g.savefig(f'{op}_pmrm_ani.pdf')
    plt.close()


    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--pmrm_input_file',
        help='Please specify the pmrm input file!',
        metavar='.',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file_prefix',
        help='Please specify the output file prefix!',
        metavar='.',
        type=str,
        required=True
        )
    parser.add_argument(
        '-xmin', '--xaxis_minimum',
        help='OPTIONAL: Minimum value to plot on x-axis. (Default=95.0)',
        metavar='.',
        type=float,
        default=80.0,
        required=False
        )
    parser.add_argument(
        '-xmax', '--xaxis_maximum',
        help='OPTIONAL: Maximum value to plot on x-axis. (Default=100.0)',
        metavar='.',
        type=float,
        default=100.0,
        required=False
        )
    parser.add_argument(
        '-xt', '--xaxis_step_size',
        help='OPTIONAL: x-axis ticks step increment. (Default=1.0)',
        metavar='.',
        type=float,
        default=5.0,
        required=False
        )
    parser.add_argument(
        '-ymin', '--yaxis_minimum',
        help='OPTIONAL: Minimum value to plot on y-axis. (Default=0.0)',
        metavar='.',
        type=float,
        default=0,
        required=False
        )
    parser.add_argument(
        '-ymax', '--yaxis_maximum',
        help='OPTIONAL: Maximum value to plot on y-axis. (Default=5.0)',
        metavar='.',
        type=float,
        default=240,
        required=False
        )
    parser.add_argument(
        '-yt', '--yaxis_step_size',
        help='OPTIONAL: y-axis ticks step increment. (Default=1.0)',
        metavar='.',
        type=float,
        default=40.0,
        required=False
        )
    parser.add_argument(
        '-z', '--add_density_layer',
        help='OPTIONAL: Input -z True to add density layer (Default=None).',
        metavar='',
        type=str,
        default=None,
        required=False
        )
    parser.add_argument(
        '-p', '--point_size',
        help='OPTIONAL: Plot point size (Default=4.0) [Not for -z True]',
        metavar='.',
        type=float,
        default=4.0,
        required=False
        )
    parser.add_argument(
        '-a', '--point_alpha',
        help='OPTIONAL: Plot point alpha value (Default=0.10) [Not for -z True]',
        metavar='.',
        type=float,
        default=0.10,
        required=False
        )
    parser.add_argument(
        '-g', '--generate_GAM_trendline',
        help='OPTIONAL: Input -g True to add trendline with GAM (Default=None).',
        metavar='',
        type=str,
        default=None,
        required=False
        )
    parser.add_argument(
        '-m', '--plot_mean_median',
        help='OPTIONAL: Input -m True to plot mean and median (Default=None).',
        metavar='',
        type=str,
        default=None,
        required=False
        )
    parser.add_argument(
        '-lg', '--log_yaxis',
        help='OPTIONAL: Set -lg True to plot ln(yaxis) (Default=None).',
        metavar='',
        type=str,
        default=None,
        required=False
        )
    parser.add_argument(
        '-s', '--species',
        help='OPTIONAL: Species label for the figure (Default: My Species)!',
        metavar='',
        type=str,
        nargs='+',
        default=['My', 'Species'],
        required=False
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n')

    # define parameters
    pmrm_file = args['pmrm_input_file']
    op = args['output_file_prefix']
    xmin = args['xaxis_minimum']
    xmax = args['xaxis_maximum']
    xt = args['xaxis_step_size']
    ymin = args['yaxis_minimum']
    ymax = args['yaxis_maximum']
    yt = args['yaxis_step_size']
    z = args['add_density_layer']
    p = args['point_size']
    a = args['point_alpha']
    gam = args['generate_GAM_trendline']
    mm = args['plot_mean_median']
    lg = args['log_yaxis']
    s = ' '.join(args['species'])

    # read in pmrm
    df = parse_pmrm_file(pmrm_file, lg)

    # build plots
    _ = pmrm_ani_plot(df, xmin, xmax, xt, ymin, ymax, yt, z, p, a, gam, mm, lg, op, s)
                               
    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')

if __name__ == "__main__":
    main()
