 #!/usr/bin/env python

'''fastANI visual hist kde plot.

This script reads the visual output from fastANI and plots the
distribution of fragment sequence identity values.

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
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np, pandas as pd; np.random.seed(0)
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


def gather_data(ani):
    """Reads in the tsv ALLvsAll ANI """

    xs = []

    with open(ani, 'r') as f:
        for l in f:
            X = l.rstrip().split('\t')
            ani = float(X[2])
            xs.append(ani)

    return xs


def fastANI_visual_dist_plot(xs, species, outfile):
    """Takes the data and builds the plot"""

    # Set Colors and markers
    main_color = '#933b41'
    second_color = '#737373'
    vline_color = '#000000'
    color = '#252525'

    # build plot
    g = sns.displot(x=xs, kde=True, color=color)

    # plot title, labels, text
    species_name = ' '.join(species.split('_'))
    g.ax.set_title(species_name, fontsize=18, y=1.02)

    g.ax.set_xlabel(
        'Fragment sequence similarity (%)',
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

    g.fig.set_figwidth(7)
    g.fig.set_figheight(5)

    # adjust layout, save, and close
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
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n')

    # define parameters
    infile = args['input_file']
    outfile = args['output_file']
    species = args['species_name']

    # read in the data
    xs = gather_data(infile)

    # build the plot
    _ = fastANI_visual_dist_plot(xs, species, outfile)

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')

if __name__ == "__main__":
    main()
