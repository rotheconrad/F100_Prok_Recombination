 #!/usr/bin/env python

'''blastANI (ani.rb) shared genome fraction vs ANI correlation.

This script reads blastANI output of one vs many or all vs all and builds
a scatter plot with the mean, median and correlation. The y-axis is the
shared genome fraction and the x-axis is the ANI. The blastANI output
includes the count of bidirectional fragment mappings and the total query
fragments along with the ANI estimate. The shared genome fraction is the
bidirectional fragment mapping count / the total fragments.

mean, meadian and correlation of ANI values (x-axis) and Shared / Total
Fragments (y-axis). Returns a scatter plot of the data as a pdf file.

This script takes the following input parameters:

    * all vs all ani file from blastANI from 06b_ANIrb_to_TSV.py (str)
    * organism name for the plot title (str)
    * An output file name. Use .pdf (str)

This script returns a publication ready figure in pdf format.

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
Date Created :: April 2022
License :: GNU GPLv3
Copyright 2022 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse, random
from collections import defaultdict
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np, pandas as pd; np.random.seed(0)
import seaborn as sns; sns.set(style="white", color_codes=True)
from scipy.stats import pearsonr as corr
from pygam import LinearGAM


"""
# Set the default sans-serif font to Helvetica
matplotlib.rcParams['font.sans-serif'] = "Helvetica"
# Set default font to sans-serif
matplotlib.rcParams['font.family'] = "sans-serif"
"""

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


def gather_data(ani_file, xmin, xmax):
    """Reads in the tsv ALLvsAll ANI """

    print("\nReading data.")
    data_dict = {'xs': [], 'ys': []}

    with open(ani_file, 'r') as f:
        # grab header
        header = f.readline()
        for l in f:
            # split each line by tabs
            X = l.rstrip().split('\t')
            # for the case where allv ani.rb doesn't full finish
            if len(X) != 10:
                print(l)
                continue
            # fill params for each line
            g1 = X[0]
            g2 = X[1]
            ani1 = float(X[2])
            ani2 = float(X[3])
            ani3 = float(X[4])
            frag1 = int(X[5])
            frag2 = int(X[6])
            frag3 = int(X[7])
            tfrag1 = int(X[8])
            tfrag2 = int(X[9])
            # chose which ani and frags to use
            ani = ani3
            ratio = frag3 / tfrag1
            # store data
            data_dict['xs'].append(ani)
            data_dict['ys'].append(ratio)

    # convert to dataframe
    df = pd.DataFrame(data_dict)
    df = df[df['xs'] <= xmax] 
    df = df[df['xs'] >= xmin]
    n = len(df)

    return df, n


def gather_stats(df):
    """Computes correlation, mean, and median on df columns xs and ys """

    # Compute Pearson Correlation Coefficient
    print("\nCalculating statistics.")
    pcorr = corr(df['xs'], df['ys'])

    # Compute ANI mean and median
    ani_mean = np.mean(df['xs'])
    ani_median = np.median(df['xs'])
    frag_mean = np.mean(df['ys'])
    frag_median = np.median(df['ys'])

    # Compile dictionairy
    df_stats = {
        'pcorr': pcorr,
        'ani_mean': ani_mean,
        'ani_median': ani_median,
        'frag_mean': frag_mean,
        'frag_median': frag_median
        }

    print(f"\nANI mean: {ani_mean:.2f}\nANI median: {ani_median:.2f}")
    print(f"\nFrag mean: {frag_mean:.2f}\nFrag median: {frag_median:.2f}")

    return df_stats


def fastANI_scatter_plot(
                df, n, species, outfile, xmin, xmax, xstep, p , a, z, g
                ):
    """Takes the data and builds the plot"""

    # Gather Stats
    df_stats = gather_stats(df)

    stats_line = (
        f"Pearson r: {round(df_stats['pcorr'][0], 2)}\n"
        f"p value: {round(df_stats['pcorr'][1], 2)}"
        )

    # Set Colors and markers
    main_color = '#933b41'
    second_color = '#737373'
    vline_color = '#000000'
    color = '#252525'
    marker = '.' #'o'

    # build plot
    g = sns.JointGrid(x="xs", y="ys", data=df)
    print('\nComputing KDEs for marginal plots.')
    # x margin kde plot
    sns.kdeplot(
            x=df["xs"],
            ax=g.ax_marg_x,
            legend=False,
            color=color
            )
    # y margin kde plot
    sns.kdeplot(
            y=df["ys"],
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
                        ds.Point("xs", "ys"),
                        ds.count(),
                        norm="log",
                        aspect="auto",
                        ax=g.ax_joint,
                        width_scale=3.,
                        height_scale=3.
                        )

    else: # regular scatter plot
        print('\nPlotting data.')
        rast = True if n >= 5000 else False
        g.ax_joint.plot(
                df["xs"],
                df["ys"],
                marker,
                ms=p,
                alpha=a,
                color=color,
                rasterized=rast,
                )

    if g:
        # Trendline with pyGAM
        print('\nCalculating trendline with pyGAM.')
        X = df["xs"].to_numpy()
        X = X[:, np.newaxis]
        y = df["ys"].to_list()

        gam = LinearGAM().gridsearch(X, y)
        XX = gam.generate_X_grid(term=0, n=500)

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
    species_name = ' '.join(species.split('_'))
    ptitle = f'{species_name} (n={n})'
    g.ax_marg_x.set_title(ptitle, fontsize=18, y=1.02)

    g.ax_joint.set_xlabel(
        'Average nucleotide identity (%)',
        fontsize=12, y=-0.02
        )
    g.ax_joint.set_ylabel(
        'Shared / total fragments',
        fontsize=12, x=-0.02
        )
    g.ax_joint.text(
        0.25, 0.99, stats_line,
        fontsize=10, color=second_color,
        verticalalignment='top', horizontalalignment='right',
        transform=g.ax_joint.transAxes
        )

    # set the axis parameters / style
    hstep = xstep/10
    g.ax_joint.set_xticks(np.arange(xmin, xmax+hstep, xstep))
    g.ax_joint.set_xlim(left=xmin-hstep, right=xmax+hstep)

    ymin = round(df['ys'].min(),2)
    ystep = round((1-ymin)/10, 2)

    g.ax_joint.set_yticks(np.arange(ymin, 1.1, ystep))
    g.ax_joint.set_ylim(bottom=ymin-ystep, top=1+ystep/2)
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

    """
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
    """

    # adjust layout, save, and close
    g.fig.set_figwidth(7)
    g.fig.set_figheight(5)
    g.savefig(outfile)
    plt.close()


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_file',
        help='Please specify the input file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-s', '--species_name',
        help='Species name for plot title. ex: Escherichia_coli',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file',
        help='Please specify the output file.',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-xmin', '--xaxis_minimum',
        help='OPTIONAL: Minimum value to plot on x-axis. (Default=95.0)',
        metavar='',
        type=float,
        default=95.0,
        required=False
        )
    parser.add_argument(
        '-xmax', '--xaxis_maximum',
        help='OPTIONAL: Maximum value to plot on x-axis. (Default=100.0)',
        metavar='',
        type=float,
        default=100.0,
        required=False
        )
    parser.add_argument(
        '-t', '--xaxis_step_size',
        help='OPTIONAL: X-axis ticks step increment. (Default=1.0)',
        metavar='',
        type=float,
        default=1.0,
        required=False
        )
    parser.add_argument(
        '-p', '--point_size',
        help='OPTIONAL: Size of the plotted points (Default=2.0)',
        metavar='',
        type=float,
        default=2.0,
        required=False
        )
    parser.add_argument(
        '-a', '--point_alpha',
        help='OPTIONAL: Alpha value of the plotted points (Default=0.10)',
        metavar='',
        type=float,
        default=0.10,
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
        '-g', '--generate_GAM_trendline',
        help='OPTIONAL: Input -g True to add trendline with GAM (Default=None).',
        metavar='',
        type=str,
        default=None,
        required=False
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n')

    # define parameters
    infile = args['input_file']
    outfile = args['output_file']
    species = args['species_name']
    xmin = args['xaxis_minimum']
    xmax = args['xaxis_maximum']
    xstep = args['xaxis_step_size']
    p = args['point_size']
    a = args['point_alpha']
    z = args['add_density_layer']
    g = args['generate_GAM_trendline']

    # read in the data
    df, n = gather_data(infile, xmin, xmax)

    # build the plot
    _ = fastANI_scatter_plot(
                df, n, species, outfile, xmin, xmax, xstep, p, a, z, g
                )

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')

if __name__ == "__main__":
    main()

