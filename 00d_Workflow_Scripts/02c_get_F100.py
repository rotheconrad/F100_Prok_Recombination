 #!/usr/bin/env python

'''Calculate F100 and ANI from concatenated RBMs

This script reads through a concatenated RBM file and calculates F100 and
ANI for each genome pair. It writes a new tsv file with columns:
Genome 1, Genome 2, ANI, F100, F99.5, F99

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
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns; sns.set(style="white", color_codes=True)

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


def fastANI_visual_dist_plot(xs, title, outfile, stat):
    """Takes the data and builds the plot"""

    # Set Colors and markers
    main_color = '#933b41'
    second_color = '#737373'
    vline_color = '#000000'
    color = '#252525'

    # build plot
    g = sns.displot(x=xs, kde=True, color=color)

    # plot title, labels, text
    g.ax.set_title(title, fontsize=18, y=1.02)

    g.ax.set_xlabel(
        'RBM sequence similarity (%)',
        fontsize=12, y=-0.02
        )
    g.ax.set_ylabel(
        'Count',
        fontsize=12, x=-0.02
        )

    # set the axis parameters / style
    g.ax.tick_params(axis='both', labelsize=12)
    g.ax.tick_params(
        axis='both', which='major', direction='inout', color='k',
        width=2, length=6, bottom=True, left=True, zorder=3
        )

    # set grid style
    g.ax.grid(
        which="major", color='#d9d9d9', linestyle='--', linewidth=1
        )
    g.ax.grid(
        which="major", color='#d9d9d9', linestyle='--', linewidth=1
        )
    g.ax.set_axisbelow(True)
    g.ax.add_artist(TickRedrawer())

    # add ANI line
    ani, f100, f995, f990 = stat[0], stat[1], stat[2], stat[3]
    stext = (
            f'ANI: {ani:.2f}\nF100: {f100:.2f}\n'
            f'F99.5: {f995:.2f}\nF99: {f990:.2f}'
            )
    g.ax.text(
        0.25, 0.99, stext,
        fontsize=10, color=second_color,
        verticalalignment='top', horizontalalignment='right',
        transform=g.ax.transAxes
        )  
    _ = g.ax.axvline(
        x=ani, ymin=0, ymax=1,
        color=vline_color, linewidth=2, linestyle='--',
        )

    g.fig.set_figwidth(7)
    g.fig.set_figheight(5)

    # adjust layout, save, and close
    g.savefig(outfile, dpi=300)
    plt.close()


def process_concatenated_rbm(infile, outpre):
    """ Reads the file does the work """

    #initialize dictionAAary to store data
    data = defaultdict(list)
    stats = {}

    # read through file
    with open(infile, 'r') as ifile:
        for line in ifile:
            # split up each line
            X = line.rstrip().split('\t')
            g1 = '_'.join(X[0].split('_')[:-2]) # grab genome 1 name
            g2 = '_'.join(X[1].split('_')[:-2]) # grab genome 2 name
            pid = float(X[2]) # grab sequence similarity
            # if g1 = g2 skip these lines
            if g1 == g2: continue
            # create list of genome names and sort 
            names = [g1, g2]
            names.sort()
            # store genome pair in data
            data['-'.join(names)].append(pid)

    # open the output file
    with open(f"{outpre}_F100.tsv", 'w') as ofile:
        header = "Genome 1\tGenome 2\tANI\tF100\tF99.5\tF99\n"
        ofile.write(header)
        # read through data dict
        for gpair, values in data.items():
            # calculate metrics
            rbms = len(values)
            ANI = sum(values) / rbms
            F100 = len([i for i in values if i >= 100]) / rbms
            F995 = len([i for i in values if i >= 99.5]) / rbms
            F990 = len([i for i in values if i >= 99]) / rbms
            # write out the line
            genomes = gpair.split('-')
            genome1 = genomes[0]
            genome2 = genomes[1]
            lineout = f'{genome1}\t{genome2}\t{ANI}\t{F100}\t{F995}\t{F990}\n'
            ofile.write(lineout)
            # store stats
            stats[gpair] = [ANI, F100, F995, F990]

    return data, stats


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
        help='Please specify the output file prefix?',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-p', '--plot_histograms',
        help='(OPTIONAL) set -p True to plot RBM histograms (default False)',
        metavar='',
        type=str,
        required=False,
        default=None
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n')

    # define parameters
    infile = args['input_file']
    outpre = args['output_file_prefix']
    plots = args['plot_histograms']

    # run it
    data, stats = process_concatenated_rbm(infile, outpre)

    # build the plots
    if plots:
        for gpair, xs in data.items():
            outfile = f'{outpre}_{gpair}.pdf'
            stat = stats[gpair]
            _ = fastANI_visual_dist_plot(xs, gpair, outfile, stat)

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')

if __name__ == "__main__":
    main()
