#!/usr/bin/env python

''' Compare fastANI vs blastANI results (calculates difference)

Reads and filters equivalent all vs all fastANI and BlastN tabular 
output files. Calculates the differences between the scores of pairwise
comparisons.

Reads in files, stores ANI and shared frac by gene pair names and then
computes the difference between. Creates scatter plot of differences.

The blastANI file is a tsv file output by 06b_ANIrb_to_TSV.py
ani.rb from enveomics is used to calculate blast based ANI. Output from
ani.rb is saved to a text file and parsed to a tsv w/ 06b_ANIrb_to_TSV.py

The fastANI file is output straight from fastANI.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: September 2022
License :: GNU GPLv3
Copyright 2022 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
import matplotlib
import matplotlib.pyplot as plt
import numpy as np, pandas as pd; np.random.seed(0)
import seaborn as sns; sns.set(style="white", color_codes=True)
from scipy.stats import pearsonr as corr


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


def parse_fastani(fastin, min_ani, min_frac):

    # Reads fastani file, filters results, returns dict

    # {read-pair: [ani, shared_frac]}
    fast_dict = {}

    with open(fastin, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            qry = X[0].split('/')[-1].split('.')[0]
            ref = X[1].split('/')[-1].split('.')[0]
            # check for self matches and discard
            if qry == ref: continue
            ani = float(X[2])
            shared_frags = int(X[3])
            total_frags = int(X[4])
            shared_frac = shared_frags / total_frags

            if shared_frac >= min_frac and ani >= min_ani:
                # check for reciprocal matches and discard lower score
                sorted_pair = '-'.join(sorted([qry, ref]))
                if sorted_pair in fast_dict:
                    old_ani = fast_dict[sorted_pair][0]
                    if ani > old_ani:
                        fast_dict[sorted_pair] = [ani, shared_frac]
                else:
                    fast_dict[sorted_pair] = [ani, shared_frac]

    return fast_dict


def compute_difference(fastANI, blastANI):

    # compute difference per genome pair between fastANI and blastANI
    # stores data in dict of {"ani": [], "frac": []}
    # returns pandas dataframe

    data = {"ani": [], "frac": []}

    for gpair, values in fastANI.items():
        if gpair in blastANI:
            # retrieve values
            fani = values[0]
            ffrac = values[1]
            bani = blastANI[gpair][0]
            bfrac = blastANI[gpair][1]
            # compute difference and store
            x = fani - bani
            y = ffrac - bfrac
            data["ani"].append(x)
            data["frac"].append(y)
        else:
            print(f"{gpair} not in blastANI data")

    for gpair in blastANI.keys():
        if gpair not in fastANI:
            print(f"{gpair} not in fastANI data")

    df = pd.DataFrame(data)

    return df


def gather_stats(df):

    # computes correlation, mean, and median

    # Compute Pearson Correlation Coefficient
    print("\nCalculating statistics.")
    pcorr = corr(df['ani'], df['frac'])

    # Compute ANI mean and median
    ani_mean = np.mean(df['ani'])
    ani_median = np.median(df['ani'])
    frag_mean = np.mean(df['frac'])
    frag_median = np.median(df['frac'])

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


def build_plots(df, title, outfile, p, a, z):

    # builds scatter plot of difference between fastani - blastani
    # shared genome fraction difference on y axis
    # ani difference on x axis

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
    g = sns.JointGrid(x="ani", y="frac", data=df)
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
            y=df["frac"],
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
                        ds.Point("ani", "frac"),
                        ds.count(),
                        norm="log",
                        aspect="auto",
                        ax=g.ax_joint,
                        width_scale=3.,
                        height_scale=3.
                        )

    else: # regular scatter plot
        print('\nPlotting data.')
        rast = True if len(df) >= 5000 else False
        g.ax_joint.plot(
                df["ani"],
                df["frac"],
                marker,
                ms=p,
                alpha=a,
                color=color,
                rasterized=rast,
                )

    # plot title, labels, text
    g.ax_marg_x.set_title(title, fontsize=18, y=1.02)

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
    #hstep = anitep/10
    #g.ax_joint.set_xticks(np.arange(xmin, xmax+hstep, xstep))
    #g.ax_joint.set_xlim(left=xmin-hstep, right=xmax+hstep)

    #ymin = round(df['frac'].min(), 2)
    #ystep = round((1-ymin)/10, 2)

    #g.ax_joint.set_yticks(np.arange(ymin, 1.1, ystep))
    #g.ax_joint.set_ylim(bottom=ymin-ystep, top=1+ystep/2)
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
    g.savefig(outfile)
    plt.close()

    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-f1', '--input_fastani_file1',
        help='Please specify the input fastANI file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-f2', '--input_fastani_file2',
        help='Please specify the input fastANI file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file_name',
        help='Please specify the output file name (Use .pdf)!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-t', '--plot_title',
        help='(OPTIONAL) Please specify the plot title!',
        metavar='',
        type=str,
        required=False
        )
    parser.add_argument(
        '-min_ani', '--minimum_ani_value',
        help='(OPTIONAL) Specify the minimum ani value (Default = 0).',
        metavar='',
        type=int,
        required=False,
        default=0
        )
    parser.add_argument(
        '-min_frac', '--minimum_aligned_frac',
        help='(OPTIONAL) Specify the minimum aligned_frac (Default = 0.0).',
        metavar='',
        type=float,
        required=False,
        default=0.0
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
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n')

    # define parameters
    fastin1 = args['input_fastani_file1']
    fastin2 = args['input_fastani_file2']
    outfile = args['output_file_name']
    title = args['plot_title']
    min_ani = args['minimum_ani_value']
    min_frac = args['minimum_aligned_frac']
    p = args['point_size']
    a = args['point_alpha']
    z = args['add_density_layer']

    if not title:
        f1 = fastin1.split('/')[-1].split('.')[0]
        f2 = fastin2.split('/')[-1].split('.')[0]
        title = f"{f1} - {f2}"

    # parse the fastANI file - returns dict {read_pair: [ANI, shared_fraction]}
    print('\nParsing fastANI file ...\n')
    fastANI1 = parse_fastani(fastin1, min_ani, min_frac)
    fastANI2 = parse_fastani(fastin2, min_ani, min_frac)

    # compare the two
    df = compute_difference(fastANI1, fastANI2)

    # build some plots
    print('\nBuilding plots ...')
    _ = build_plots(df, title, outfile, p, a, z)

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')


if __name__ == "__main__":
    main()
