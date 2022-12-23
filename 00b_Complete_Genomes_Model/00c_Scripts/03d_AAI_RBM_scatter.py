 #!/usr/bin/env python

'''AAI.rb RBM F100 vs ANI scatter plot.

F100 = (count of 100% sequence similarity) / (total RBMs)

This script makes a scatter plot of the tsv from 03c_AAI_RBM_F100.py
with F100 (y-axis) and ANI (x-axis).

This script returns a publication ready figure in pdf format.

This script requires python 3.6+ and the following modules:

    * matplotlib
    * numpy
    * pandas
    * seaborn
    * scipy
    * datashader for -z True density plot
        - pip install datashader
        - conda install datashader

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

import argparse
#import os
import matplotlib
matplotlib.use('Agg')
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


def gather_data(infile, xmin, xmax, D):

    data_dict = {'xs': [], 'ys': []}

    with open(infile, 'r') as f:
        for l in f:
            X = l.rstrip().split('\t')
            gpair = X[0]
            ani = float(X[1])

            # convert ANI to D (100 - ANI)
            if D:
                ani = 1 - (ani/100)

            f100 = float(X[2])

            data_dict['xs'].append(ani)
            data_dict['ys'].append(f100)

    df = pd.DataFrame(data_dict)

    # Remove matches below xmin and above xmax ANI
    if D:
        df = df[df['xs'] <= 1 - (xmin/100)]
        df = df[df['xs'] >= 1 - (xmax/100)]
        print(f"\nMinimum Distance in data: {df['xs'].min():.2f}")
        print(f"Maximum Distance in data: {df['xs'].max():.2f}")
    else:
        df = df[df['xs'] <= xmax]
        df = df[df['xs'] >= xmin]
        print(f"\nMinimum ANI in data: {df['xs'].min():.2f}")
        print(f"Maximum ANI in data: {df['xs'].max():.2f}")

    gcount = len(df)

    return df, gcount
    

def f100_scatter_plot(
            df, n, species, outfile, xmin, xmax, xstep, log, D, p, a, z
            ):
    """Takes the data and builds the plot"""

    # Set Colors and markers
    main_color = '#933b41'
    second_color = '#737373'
    vline_color = '#000000'
    color = '#252525'
    marker = '.' #'o'

    # setup plot
    g = sns.JointGrid(x="xs", y="ys", data=df)
    g.fig.set_figwidth(7)
    g.fig.set_figheight(5)

    # add data
    sns.kdeplot(
            x=df["xs"],
            ax=g.ax_marg_x,
            legend=False,
            color=color
            )
    sns.kdeplot(
            y=df["ys"],
            ax=g.ax_marg_y,
            #vertical=True,
            legend=False,
            color=color
            )

    if z:
        import datashader as ds
        from datashader.mpl_ext import dsshow
        
        dsartist = dsshow(
                        df,
                        ds.Point("xs", "ys"),
                        ds.count(),
                        #vmin=0,
                        #vmax=100,
                        norm="log",
                        aspect="auto",
                        ax=g.ax_joint
                        )
        g.ax_joint.plot(
                df["xs"],
                df["ys"],
                marker,
                ms=0,
                alpha=0,
                rasterized=True,
                )

    else:
        rast = True if n >= 1000 else False
        g.ax_joint.plot(
                df["xs"],
                df["ys"],
                marker,
                ms=p,
                alpha=a,
                color=color,
                rasterized=rast,
                #markeredgecolor='#000000'
                )

    # plot title, labels, text
    species_name = ' '.join(species.split('_'))
    ptitle = f'{species_name} (n={n})'
    g.ax_marg_x.set_title(ptitle, fontsize=18, y=1.02)

    if D and log:
        xlab = 'Distance (1 - (ANI/100))'
        xpos, ypos = 0.25, 0.6
    elif D:
        xlab = 'Distance (1 - ANI/100))'
        xpos, ypos = 0.99, 0.99
    elif log:
        xlab = 'Average nucleotide identity (%)'
        xpos, ypos = 0.25, 0.6
    else:
        xlab = 'Average nucleotide identity (%)'
        xpos, ypos = 0.25, 0.99

    g.ax_joint.set_xlabel(xlab, fontsize=12, y=-0.02)
    g.ax_joint.set_ylabel(
        "F100 (100% RBMs) / (total RBMs)",
        fontsize=12, x=-0.02
        )

    # set the axis parameters / style    
    if D and log:
        g.ax_joint.set_xscale('log')
    elif D:
        hstep = xstep/10
        g.ax_joint.set_xticks(np.arange(-hstep, xmax+hstep, hstep), minor=True)
        g.ax_joint.set_xticks(np.arange(0, xmax, xstep))
        g.ax_joint.set_xlim(left=-hstep, right=xmax+hstep)
    elif log:
        g.ax_joint.set_xscale('log')
    else:
        hstep = xstep/10
        g.ax_joint.set_xticks(np.arange(xmin-hstep, xmax+hstep, hstep), minor=True)
        g.ax_joint.set_xticks(np.arange(xmin, xmax+hstep, xstep))
        g.ax_joint.set_xlim(left=xmin-hstep, right=xmax+hstep)

    g.ax_joint.set_yticks(np.arange(0.0, 1.1, 0.1))
    g.ax_joint.set_ylim(bottom=-0.02, top=1.02)
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

    # save, and close
    g.savefig(outfile, dpi=300)
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
        help='Please specify the output file?',
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
        '-l', '--xaxis_log_scale',
        help='OPTIONAL: Input -l True to set x-axis to log scale (Default=None).',
        metavar='',
        type=str,
        default=None,
        required=False
        )
    parser.add_argument(
        '-d', '--convert_ani_distance',
        help='OPTIONAL: Input -d True for ANI = D (1-ANI/100) (Default=None).',
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
    parser.add_argument(
        '-z', '--add_density_layer',
        help='OPTIONAL: Input -z True to add density layer (Default=None).',
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
    log = args['xaxis_log_scale']
    D = args['convert_ani_distance']
    p = args['point_size']
    a = args['point_alpha']
    z = args['add_density_layer']

    # read in the data
    df, n = gather_data(infile, xmin, xmax, D)

    # build the plot
    _ = f100_scatter_plot(
                df, n, species, outfile, xmin, xmax, xstep, log, D, p, a, z
                )

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')

if __name__ == "__main__":
    main()
