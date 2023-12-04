 #!/usr/bin/env python

'''F100 vs ANI GAM.

This script reads the tsv output from 03c_AAI_RBM_F100.py and fits a GAM
model to the data using the pyGAM package.

F100 is the ratio of (100% identical RBMs) / (Total RBMs).
RBMs are reciprocal best blast matches.
GAM is generalized linear model. 

This script takes the following input parameters:

    * f100 - all vs all f100 tsv file from 03c_AAI_RBM_F100.py (str)
    * title - Plot title (str)
    * op - An output file prefix (str)

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
    * datashader for density scatter plot if > 500 genome pairs.
        - pip install datashader
        - conda install datashader

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: June 2022
License :: GNU GPLv3
Copyright 2022 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse, sys, random
from collections import defaultdict
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np, pandas as pd; np.random.seed(0)
import seaborn as sns; sns.set(style="white", color_codes=True)
import datashader as ds
from datashader.mpl_ext import dsshow
#from scipy.stats import pearsonr as corr
from pygam import LinearGAM, s


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


def read_f100_file(f100_file, labelled):
    """Read the file and grabs the data"""
    # initialize dict for data {genome1-genome2: [data]}
    name_dict = {}
    # read through the file and populate the dict
    with open(f100_file, 'r') as file:
        # grab the header
        header = file.readline()
        for line in file:
            # split each line by tabs
            X = line.rstrip().split('\t')
            # get the genome pair names, ANI, and F100
            qry_genome = X[0]
            ref_genome = X[1]
            if qry_genome == ref_genome: continue
            ANI = float(X[2])
            F100 = float(X[3])
            # grouping meta data for plotting and stats
            clade1, clade2 = 0, 0
            if labelled:
                clade1 = int(X[6])
                clade2 = int(X[7])
                """
                if clade1 == clade2:
                    with open(f'Salruber_Clade-0{clade1}.tsv', 'a') as out:
                        out.write(line)
                """
            # sort genome file names and combine
            names = [qry_genome, ref_genome]
            names.sort()
            gname = '-'.join(names)

            # remove F100 == 0 and F100 = 1
            if F100 == 0 or F100 == 1: continue
            # Transform Y to linearize
            F100 = np.log(F100)

            # store data
            name_dict[gname] = [ANI, F100, clade1, clade2]

    return name_dict


def gather_data(f100_file, labelled, xmin, xmax):
    """Reads the file and organized the data in a dict"""

    print("\nReading data.")
    name_dict = read_f100_file(f100_file, labelled)
    data_dict = {'gpair': [], 'xs': [], 'ys': [], 'clade1': [], 'clade2': []}

    # write data to arrays
    for gpair, metrics in name_dict.items():
        ani = metrics[0]
        f100 = metrics[1]
        clade1 = metrics[2]
        clade2 = metrics[3]
        data_dict['xs'].append(ani)
        data_dict['ys'].append(f100)
        data_dict['clade1'].append(clade1)
        data_dict['clade2'].append(clade2)
        data_dict['gpair'].append(gpair)
    # convert to dataframe
    df = pd.DataFrame(data_dict)
    df = df[df['xs'] <= xmax] 
    df = df[df['xs'] >= xmin]

    print(df)
    return df

def f100_scatter_plot(
                df, df2, labelled, title, outfile, xmin, xmax, xstep, p , a, e
                ):
    """Takes the data and builds the plot"""

    # open file to write out genome pairs outside the confidence interval.
    sig_pair_out = open(f'{outfile}_sig-pairs.tsv', 'w')
    # write a file header.
    sig_pair_out.write(
                'Genome pair\tClade 1\tClade 2\tPrediction\t'
                'p value\t2.5% c.i.\t97.5% c.i.\n'
                )

    # Set Colors and markers
    grid_color = '#d9d9d9'
    main_color = '#933b41'
    second_color = '#737373'
    vline_color = '#000000'
    color = '#252525'
    marker = '.' #'o'

    # build plot
    gg = sns.JointGrid(x="xs", y="ys", data=df)

    # x margin hist plot
    sns.histplot(
            x=df["xs"],
            ax=gg.ax_marg_x,
            legend=False,
            color=color,
            stat='probability'
            )
    # y margin hist plot
    sns.histplot(
            y=df["ys"],
            ax=gg.ax_marg_y,
            legend=False,
            color=color,
            stat='probability'
            )
    # main panel scatter plot
    print('\nComputing plot densities.')
    dsartist = dsshow(
                    df,
                    ds.Point("xs", "ys"),
                    ds.count(),
                    norm="log",
                    aspect="auto",
                    ax=gg.ax_joint,
                    width_scale=3.,
                    height_scale=3.
                    )
    dsartist.zorder = 2.5

    # Trendline with pyGAM
    print('\nCalculating trendline with pyGAM.')
    X = df["xs"].to_numpy()
    X = X[:, np.newaxis]
    y = df["ys"].to_list()

    # Build and fit the GAM model using gridsearch to tune param
    lam = np.logspace(-5, 5, 10)
    spline = s(0, n_splines=5, constraints='convex')
    gam = LinearGAM(spline).gridsearch(X, y, lam=lam)
    # constraints ['convex', 'concave', 'monotonic_inc', 'monotonic_dec']
    #gam.gridsearch(X, y, lam=lam)

    print('\n', 'lam: ', gam.lam, '\n')
    #gam.summary()
    XX = gam.generate_X_grid(term=0)
    
    XX = np.linspace(e, 100, 500)

    gg.ax_joint.plot(
                XX,
                gam.predict(XX),
                color='#2171b5',
                linestyle='--',
                linewidth=1.0,
                zorder=2.8
                )
    gg.ax_joint.plot(
                XX,
                gam.prediction_intervals(XX, width=0.95),
                color='#9ecae1',
                linestyle='--',
                linewidth=1.0,
                zorder=2.8
                )
    r2 = gam.statistics_['pseudo_r2']['explained_deviance']
    GAM_line = f"GAM Pseudo R-Squared: {r2:.4f}"
    gg.ax_joint.text(
        0.1, 0.9, GAM_line,
        fontsize=10, color=second_color,
        verticalalignment='top', horizontalalignment='left',
        transform=gg.ax_joint.transAxes
        )

    # plot title, labels, text
    gg.ax_marg_x.set_title(title, fontsize=18, y=1.02)

    gg.ax_joint.set_xlabel(
        'Average nucleotide identity (%)',
        fontsize=12, y=-0.02
        )
    gg.ax_joint.set_ylabel(
        #'F100 [(100% identical RBMs) / (total RBMs)]',
        'ln(F100)',
        fontsize=12, x=-0.02
        )

    # plot second file overlay
    if df2 is not None:

        # define the marker set        
        marker_set = {
                0: 'o', 1: 'v', 2: '^', 3: '<', 4: '>', 5: 's',
                6: 'p', 7: '*', 8: 'X', 9: 'D', 10: '1', 11: '2',
                12: '3', 13: '4', 
                }
        # keep track of clades
        clades = {}
        # hold legend elements
        legend_elements = []

        # get the data in an easy to iterate fashion
        overlay = list(
                    zip(
                        df2["xs"],
                        df2["ys"],
                        df2["clade1"],
                        df2["clade2"],
                        df2["gpair"]
                        )
                    )
        # iterate throught the data
        for ox, oy, c1, c2, gpair in overlay:

            # check and color points by confidence interval
            ci = gam.prediction_intervals(ox, width=0.95)[0]

            recombines = None
            # if below the lower c.i. Non-recombining
            if oy < ci[0]:
                recombines = "Non-recombining"
            # if above the upper c.i. Recombing
            elif oy > ci[1]:
                recombines = "Recombining"

            if recombines:

                # write the genome pair to file
                sig_pair_out.write(
                            f'{gpair}\t{c1}\t{c2}\t{recombines}\t'
                            f'{oy}\t{ci[0]}\t{ci[1]}\n'
                            )

                if labelled:
                    c, marker, zorder = get_label(c1, c2, marker_set)
                    clades[c1] = ''
                else:
                    c = '#00ffbc'
                    marker = 's'
                    zorder = 2.7
                gg.ax_joint.scatter(
                            ox, oy, facecolors='none', edgecolors=c,
                            linewidth=0.5, marker=marker, s=10,
                            alpha=0.5, zorder=zorder
                            )
            else:
                if labelled:
                    c, marker, zorder = get_label(c1, c2, marker_set)
                    clades[c1] = ''
                else:
                    c = '#c493ff'
                    marker = 'o'
                    zorder = 2.7
                gg.ax_joint.scatter(
                            ox, oy, facecolors='none', edgecolors=c,
                            linewidth=0.5, marker=marker, s=10,
                            alpha=0.5, zorder=zorder
                            )

    # set the axis parameters / style
    hstep = xstep/10
    gg.ax_joint.set_xticks(np.arange(e, xmax+hstep, xstep))
    gg.ax_joint.set_xlim(left=e-hstep, right=xmax+hstep)

    gg.ax_joint.set_ylim(top=0.1)

    # set grid style
    gg.ax_joint.yaxis.grid(
        which="major", color='#d9d9d9', linestyle='--', linewidth=1
        )
    gg.ax_joint.xaxis.grid(
        which="major", color='#d9d9d9', linestyle='--', linewidth=1
        )
    gg.ax_joint.set_axisbelow(True)
    gg.ax_joint.add_artist(TickRedrawer())

    # adjust layout, save, and close
    gg.fig.set_figwidth(7)
    gg.fig.set_figheight(5)
    gg.savefig(f'{outfile}_GAMplot.pdf')
    plt.close()


def get_label(c1, c2, marker_set):
    # if the data is labelled, check within or between clade
    # within clade labelled with square marker
    # between clade keeps circle marker

    if c1 == c2:
        # within clade
        c = '#00ffbc'
        marker = 's' #marker_set[c1]
        zorder = 2.8
    else:
        # between clades
        c = '#c493ff'
        marker = 'o' #marker_set[c1]
        zorder = 2.7

    return c, marker, zorder


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
        '-o', '--output_file_prefix',
        help='Please specify the output file prefix!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-t', '--plot_title',
        help='OPTIONAL: Please specify the plot title (defalut: F100 vs ANI)!',
        metavar='',
        type=str,
        nargs='+',
        default=['F100', 'vs', 'ANI'],
        required=False
        )
    parser.add_argument(
        '-xmin', '--ANI_minimum',
        help='OPTIONAL: Minimum ANI value to include. (Default=70.0)',
        metavar='',
        type=float,
        default=70.0,
        required=False
        )
    parser.add_argument(
        '-xmax', '--ANI_maximum',
        help='OPTIONAL: Maximum ANI value to include. (Default=100.0)',
        metavar='',
        type=float,
        default=100.0,
        required=False
        )
    parser.add_argument(
        '-i2', '--input_file2',
        help='OPTIONAL: Second file data to plot on top of first file data.',
        metavar='',
        type=str,
        required=False,
        default=None
        )
    parser.add_argument(
        '-l', '--input_file2_labelled',
        help='OPTIONAL: Second input file has clade labels (default: False).',
        metavar='',
        type=str,
        required=False,
        default=None
        )
    parser.add_argument(
        '-s', '--xaxis_step_size',
        help='OPTIONAL: X-axis ticks step increment. (Default=5.0)',
        metavar='',
        type=float,
        default=5.0,
        required=False
        )
    parser.add_argument(
        '-p', '--point_size',
        help='OPTIONAL: Size of the plotted points (Default=4.0)',
        metavar='',
        type=float,
        default=4.0,
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
        '-e', '--extrapolate_min_ANI',
        help='OPTIONAL: Minimum ANI for extrapolation (Default=70)',
        metavar='',
        type=int,
        default=70,
        required=False
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n')

    # define parameters
    infile = args['input_file']
    infile2 = args['input_file2']
    labelled = args['input_file2_labelled']
    outfile = args['output_file_prefix']
    title = ' '.join(args['plot_title'])
    xmin = args['ANI_minimum']
    xmax = args['ANI_maximum']
    xstep = args['xaxis_step_size']
    p = args['point_size']
    a = args['point_alpha']
    e = args['extrapolate_min_ANI']


    # read in the data
    df = gather_data(infile, False, xmin, xmax)
    # read in second data
    if infile2:
        df2 = gather_data(infile2, labelled, xmin, xmax)
    else:
        df2 = None

    # build the plot
    _ = f100_scatter_plot(
                df, df2, labelled, title, outfile, xmin, xmax, xstep, p, a, e
                )
                               
    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')

if __name__ == "__main__":
    main()
