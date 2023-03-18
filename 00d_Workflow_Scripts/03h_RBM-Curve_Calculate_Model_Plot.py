#!/usr/bin/env python

''' Calculate and Plot a rRBMs Rarefaction Curve

Takes a tab separated binary RBM presence/absence matrix with genomes
as the columns and RBMs as the rows. Genome IDs should be on the 
first line (row), and RBM IDs should in the first column. All other
values should be either a 0 for rRBM absent in the genome (column) or
1 for rRBM is present in the genome.

rRBMs modelling based on:
Zhang et. al. 2018, 10.3389/fmicb.2018.00577
Tettelin et. al. 2008, https://doi.org/10.1016/j.mib.2008.09.006

This script requires the lmfit package.

https://lmfit.github.io/lmfit-py/
conda install -c conda-forge lmfit

-------------------------------------------
Author :: Roth Conrad & Carlos Ruiz
Email :: rotheconrad@gatech.edu, cruizperez3@gatech.edu
GitHub :: https://github.com/rotheconrad, https://github.com/cruizperez/
Date Created :: January 24th, 2020
License :: GNU GPLv3
Copyright 2020 Roth Conrad & Carlos Ruiz
All rights reserved
-------------------------------------------
'''

import argparse
from collections import defaultdict
import numpy as np
import pandas as pd
from lmfit.models import PowerLawModel, ExpressionModel
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def plot_rRBMs_curve(
                        dfout, PLM_totl, PLM_spec,
                        omega_core, omega_new,
                        prm, pc, out, ymax, ystep
                        ):
    ''' This function builds a plot of the rRBMs calculations'''

    print('\n\nBuilding a plot of the data ...')

    # Set x-axis values
    x = dfout['n'].unique()

    # Calculate mean, median, and quartiles
    dfmean = dfout.groupby('n').mean()
    dfmedian = dfout.groupby('n').median()
    df025 = dfout.groupby('n').quantile(q=0.025)
    df975 = dfout.groupby('n').quantile(q=0.975)

    # Set the colors
    H1 = '#933b41'
    totl = '#933b41'
    core = '#0868ac'
    spec = '#b35806'
    new = '#542788'
    other = '#bdbdbd'
    # alpha value to use for IQRs
    a = 0.2
    # model marker, alpha, and size
    mdl = 'd'
    m = 0.3
    md = 5

    # Build the plot
    fig, (ax1, ax2) = plt.subplots(
        2, 1,
        gridspec_kw={'height_ratios': [2, 1]},
        figsize=(20, 14),
        sharex=True,
        sharey=False
        )

    # plot title, labels, and text
    ax1.set_title(
        'Recombinant RBM (rRBMs) Curves',
        color=H1, fontsize=50, y=1.02
        )
    ax2.set_title(
        'New rRBMs per genome addition',
        color=new, fontsize=32, y=1.02
        )
    ax1.set_ylabel('rRBMs count', fontsize=28)
    ax2.set_ylabel('New RBMs (%)', fontsize=28)
    ax2.set_xlabel('Genome count', fontsize=28)

    # Emperical data text
    total_genomes=len(x)+1

    stext = (
        f"Total rRBMs: {dfmean.rRBMs.tolist()[-1]}  |  "
        f"Core rRBMs: {dfmean.CorerRBMs.tolist()[-1]}  |  "
        f"Non-recombinant RBMs: {dfmean.NonRecombinant.tolist()[-1]}"
        )
    ax1.text(
        0.5, 0.02, stext,
        fontsize=18, color='#737373',
        horizontalalignment='center', transform=ax1.transAxes
        )
    ttext = (
        f"Number of Genomes: {total_genomes}  |  "
        f"Number of Permutations: {prm}"
        )
    ax1.text(
        0.5, 0.94, ttext,
        fontsize=22, color='#737373',
        horizontalalignment='center', transform=ax1.transAxes
        )
    mtext = (
        f"Mean rRBMs per Genome: {int(dfmean.totalrRBMs.mean())}  |  "
        f"New rRBMs per Genome: {int(dfmean.NewrRBMs.mean())}  |  "
        f"New rRBMs at n={total_genomes}: {int(dfmean.NewrRBMs.tolist()[-1])}"
        )
    ax2.text(
        0.5, 0.90, mtext,
        fontsize=18, color='#737373',
        horizontalalignment='center', transform=ax2.transAxes
        )

    # Modelled data text
    totlgamma = PLM_totl['Gamma']
    specgamma = PLM_spec['Gamma']

    if totlgamma < 0: totllabel = 'Closed'
    elif totlgamma <= 1: totllabel = 'Open'
    else: print('Error in rRBMs Model Gamma Paramerter')
    if specgamma < 0: speclabel = 'Closed'
    elif specgamma <= 1: speclabel = 'Open'
    else: print('Error in Non-recombinant Model Gamma Paramerter')

    modeltext = (
        f'Total rRBMs Model \u03B3 = {totlgamma:.2f}, {totllabel}  |  '
        f'Non-recombinant Model \u03B3 = {specgamma:.2f}, {speclabel}  |  '
        f'Core rRBMs Model \u03A9 = {omega_core:.2f}'
        )

    ax1.text(
        0.5, -0.05, modeltext,
        fontsize=18, color='#000000',
        horizontalalignment='center', transform=ax1.transAxes
        )
    rtext = (
        f'Ratio at n={total_genomes}: {dfmean.NewrRBMsRatio.tolist()[-1]:.2f}%\n'
        f'New rRBMs Ratio Model \u03A9 = {omega_new:.2f}%'
        )
    ax2.text(
        len(x)-2, dfmean.NewrRBMsRatio.mean() + 5, rtext,
        fontsize=18, color=new, horizontalalignment='right'
        )


    # Plot rRBMs Medians, Means and IQR
    ax1.plot(x, dfmedian.rRBMs, color=totl, linestyle='--', lw=1)
    ax1.plot(x, dfmean.rRBMs, color=totl, linestyle=':', lw=2)
    ax1.fill_between(x, df025.rRBMs, df975.rRBMs, color=totl, alpha=a)
    ax1.plot(
        x, dfmean.rRBMs_PLM,
        color=totl, marker=mdl, linestyle='', markersize=md, alpha=m
        )

    # Plot Non-recombinant Medians, Means and IQR
    ax1.plot(x, dfmedian.NonRecombinant, color=spec, linestyle='--', lw=1)
    ax1.plot(x, dfmean.NonRecombinant, color=spec, linestyle=':', lw=2)
    ax1.fill_between(
                    x, df025.NonRecombinant, df975.NonRecombinant,
                    color=spec, alpha=a
                    )
    ax1.plot(x, dfmean.NonRecombinant_PLM,
        color=spec, marker=mdl, linestyle='', markersize=md, alpha=m
        )

    # Plot Core Medians, Means and IQR
    ax1.plot(x, dfmedian.CorerRBMs, color=core, linestyle='--', lw=1)
    ax1.plot(x, dfmean.CorerRBMs, color=core, linestyle=':', lw=2)
    ax1.fill_between(x, df025.CorerRBMs, df975.CorerRBMs, color=core, alpha=a)
    ax1.plot(x, dfmean.CorerRBMs_EDM,
        color=core, marker=mdl, linestyle='', markersize=md, alpha=m
        )

    # Build ax1 Plot Legend
    l_totl = Line2D(
        [0],[0], color='w', label='Total rRBMs',
        markerfacecolor=totl, marker='o', markersize=18, alpha=m
        )
    l_core = Line2D(
        [0],[0], color='w', label='Core rRBMs',
        markerfacecolor=core, marker='o', markersize=18, alpha=m
        )
    l_specific = Line2D(
        [0],[0], color='w', label='Non-recombinant RBMs',
        markerfacecolor=spec, marker='o', markersize=18, alpha=m
        )
    l_IQR = Line2D(
        [0],[0], color='w', label='95% IQR',
        markerfacecolor=other, marker='s', markersize=20
        )
    l_mean = Line2D(
        [0],[0], color=other, linestyle=':', lw=4, label='Mean'
        )
    l_median = Line2D(
        [0],[0], color=other, linestyle='--', lw=4, label='Median'
        )
    l_model = Line2D(
        [0],[0], color='w', label='Model Fit',
        markerfacecolor=other, marker=mdl, markersize=20
        )

    legend_elements = [
                        l_totl,
                        l_core,
                        l_specific,
                        l_IQR,
                        l_mean,
                        l_median,
                        l_model
                        ]

    ax1.legend(
        handles=legend_elements,
        loc='upper left',
        fontsize=18,
        fancybox=True,
        framealpha=0.0,
        frameon=False
        )

    # Plot rRBMs Ratio Plot
    ax2.plot(x, dfmedian.NewrRBMsRatio, color=new, linestyle='--', lw=1)
    ax2.plot(x, dfmean.NewrRBMsRatio, color=new, linestyle=':', lw=2)
    ax2.fill_between(
                x, df025.NewrRBMsRatio, df975.NewrRBMsRatio,
                color=new, alpha=a
                )
    ax2.plot(x, dfmean.NewrRBMsRatio_EDM,
        color=new, marker=mdl, linestyle='', markersize=md, alpha=m
        )

    # set the axis parameters / style
    ax1.yaxis.grid(which="both", color='#d9d9d9', linestyle='--', linewidth=1)

    # y axis minimum
    ymin = -200
    if ymax and ystep:
        ax1.set_yticks(range(0, ymax+1, ystep))
    elif ymax and not ystep:
        ax1.set_ylim(ymin, ymax)
    elif ystep and not ymax:
        start, end = ax1.get_ylim()
        ax1.set_yticks(np.arange(start, end, ystep))
    else:
        _, ymax = ax1.get_ylim()
        ax1.set_ylim(ymin, ymax)

    ax1.minorticks_on()
    ax1.set_xlabel('')
    ax1.tick_params(labelsize=22)
    for spine in ax1.spines.values(): spine.set_linewidth(2)
    ax1.set_axisbelow(True)

    ax2.yaxis.grid(which="both", color='#d9d9d9', linestyle='--', linewidth=1)
    ax2.minorticks_on()
    ax2.tick_params(labelsize=22)
    ax2.set_xticks(range(0, len(x)+2, 10))
    ax2.set_xlim(-1, len(x)+2)
    for spine in ax2.spines.values(): spine.set_linewidth(2)
    ax2.set_axisbelow(True)

    plt.subplots_adjust(
        left = 0.09,
        right = 0.98,
        bottom = 0.07,
        top = 0.87,
        hspace = 0.2
        )

    plt.savefig(f'{out}.pdf')
    plt.close()


def model_decay_curve(dfout, Column):
    ''' This function models the decay curve for core, and new rRBMs'''

    # Model the Core rRBMs curve using an Exponential Decay Function:
    # Fc = Kc*exp(-N/τc) + Ω
    print(f'\n\nFitting Exponential Decay function to {Column}')
    print('Using Exponential Decay Function: K*exp-(N/\u03C4) + \u03A9')
    # Initialize model
    Custom_Exponential = ExpressionModel(
                                    'A * exp(-x/tau) + omega',
                                    independent_vars=['x']
                                    )
    # Initialize custom parameters
    Expression_Params = Custom_Exponential.make_params()
    # add params with tuples: (NAME VALUE VARY MIN  MAX  EXPR  BRUTE_STEP)
    Expression_Params.add_many(
                            ('A', 5, True, 0, None, None, 0.1),
                            ('tau', 5, True, 0, None, None, 0.1),
                            ('omega', 5, True, 0, None, None, 0.1)
                            )
    EDM_Fit = Custom_Exponential.fit(
                                    dfout[Column],
                                    Expression_Params,
                                    x=dfout['n']
                                    )
    dfout[f'{Column}_EDM'] = EDM_Fit.best_fit
    omega = EDM_Fit.best_values['omega']

    return dfout, omega


def model_rRBMs_curve(dfout, Column):
    ''' This function models the total, core, and new rRBMs curves'''

    # Initialize dictionary to store model results
    params = {} # PowerLawModel for rRBMs curve
    # Model the rRBMs curve using a Powerlaw function: Ps = κn^γ
    print('\n\nFitting PowerLaw function to {Column}')
    print('Using Power Law Function K*N^\u03B3 ...')
    PLM = PowerLawModel() # Initialize the model
    # Guess starting values from the data
    PLM_Guess_Start = PLM.guess(
                            dfout[Column],
                            x=dfout['n'],
                            )
    # Fit model to the data (optimize parameters)
    PLM_Fit = PLM.fit(
                    dfout[Column],
                    PLM_Guess_Start,
                    x=dfout['n']
                    )
    # Add modeled data to dfout data frame
    dfout[f'{Column}_PLM'] = PLM_Fit.best_fit
    # Store K and Gamma parameters 
    params['K'] = float(PLM_Fit.best_values['amplitude'])
    params['Gamma'] = float(PLM_Fit.best_values['exponent'])

    return dfout, params

def calculate_rRBMs_curve(binary_matrix, prm, c, out):
    # Read in the binary matrix with first column as RBM names (index).
    df = pd.read_csv(binary_matrix, sep='\t') #, index_col=0)
    # Define columns for output data frame
    colout = [
            'Trial', 'n', 'n/N', 'rRBMs', 'CorerRBMs',
            'NonRecombinant', 'NewrRBMs', 'NewrRBMsRatio', 'totalrRBMs'
            ]
    # Initialize output data frame
    dfout = pd.DataFrame(columns=colout)
    # Initialize row counter for ouput dataframe
    rowcount = 0
    # Set the total number of genomes to n
    N = df.shape[1]
    # Calculate the total number of rRBMs per genome
    rRBMs_per_genome = df.sum(axis=0)

    for j in range(prm):

        print('Running Permutation:', j+1)
        # Here I initialize(empty) a dataframe for each permutation.
        # Then I randomly shuffle the gneome order for each permutation.

        # Initialize dataframe to keep track of permutation
        dfprm = pd.DataFrame()
        # Random shuffle the genomes for each permutation
        dfrand = df.sample(frac=1, axis=1, random_state=j*42)
        # Set convenient list of genome names (column names)
        genomes = dfrand.columns.tolist()

        # Add first genome to dfprm
        gnm = genomes[0]
        dfprm[gnm] = dfrand[gnm]
        # Start the rRBM curve count
        rRBMcurve = dfprm.sum(axis=1).values
        rRBMprev = (rRBMcurve != 0).sum()

        for n in range(1,N):
            # Then for each permutation we step through the columns
            # And calculate values for each genome addition.

            # Set new value for current number of genomes
            nN = 1 / (n+1)
            # select current genome name from prm genomes list
            gnm = genomes[n]
            # add current genome column to dfprm to step the iteration
            dfprm[gnm] = dfrand[gnm]
            # calculate the n/N rRBMs values at this step
            rRBMcurve = dfprm.sum(axis=1).values / (n+1)

            # calculate total rRBMs size
            totl = (rRBMcurve != 0).sum()
            # calculate core rRBMs size
            core = (rRBMcurve >= c).sum()
            # calculate Non-recombinant RBMs
            nonrec = (rRBMcurve == 0).sum()
            # Grab number of rRBMs in current genome
            totalrRBMs = rRBMs_per_genome[gnm]
            # Calculate new rRBMs per genome / the number of rRBMs in genome
            nrRBMs = totl - rRBMprev
            nratio = nrRBMs / totalrRBMs * 100
            rRBMprev = totl

            # Define new row for dfout dataframe.
            z = [j+1, n+1, nN, totl, core, nonrec, nrRBMs, nratio, totalrRBMs]
            # Add new row to dfout dataframe
            dfout.loc[rowcount] = z
            # increment the rowcount
            rowcount += 1

    dfout.to_csv(f'{out}_data.tsv', sep='\t')

    return dfout

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-b', '--binary_matrix',
        help='Please specify the binary.tsv input file!',
        #metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-p', '--number_permutations',
        help='Optional: Number of permutations to run (default: 100)',
        #metavar='',
        type=int,
        default=100,
        required=False
        )
    parser.add_argument(
        '-c', '--percent_core',
        help='Optional: Genome fraction with RBM for core (default: 1.0).',
        #metavar='',
        type=float,
        default=1.0,
        required=False
        )
    parser.add_argument(
        '-y', '--yaxis_max',
        help='Optional: Set range of y-axis rRBMs size (eg: 20000)',
        #metavar='',
        type=int,
        required=False
        )
    parser.add_argument(
        '-s', '--yaxis_step',
        help='Optional: Set y-axis step increment (eg: 2000)',
        #metavar='',
        type=int,
        required=False
        )
    parser.add_argument(
        '-o', '--output_prefix',
        help='What do you want to name the output files?',
        #metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Run this scripts main function
    print('\n\nCalculating the rarefaction curve ...\n')

    dfout = calculate_rRBMs_curve(
                            args['binary_matrix'],
                            args['number_permutations'],
                            args['percent_core'],
                            args['output_prefix']
                            )

    # Fit a Power Law Model to the rRBMs Curve
    dfout, PLM_totl = model_rRBMs_curve(dfout, 'rRBMs')
    # Fit a Power Law Model to the Non-recombinant RBMs
    dfout, PLM_spec = model_rRBMs_curve(dfout, 'NonRecombinant')
    # Fit an Exponential Decay Model to the Core Genome
    dfout, omega_core = model_decay_curve(dfout, 'CorerRBMs')
    # Fit an Exponential Decay Model to the New rRBMs per genome count
    dfout, omega_new = model_decay_curve(dfout, 'NewrRBMsRatio')
    # Plot the results
    plot_rRBMs_curve(
                        dfout,
                        PLM_totl,
                        PLM_spec,
                        omega_core,
                        omega_new,
                        args['number_permutations'],
                        args['percent_core'],
                        args['output_prefix'],
                        args['yaxis_max'],
                        args['yaxis_step']
                        )

    print(
        '\n\nCongratulations!! '
        'The script seems to have finished successfully.\n\n'
        )

if __name__ == "__main__":
    main()