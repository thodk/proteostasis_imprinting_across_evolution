#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
import numpy
import pandas
import itertools
import sys
sys.path.append('../')
from core_functions import taxonomy_colors
from core_functions import taxonomy_pairs_colors
from scipy.stats import gaussian_kde
from matplotlib import pyplot as plt
import matplotlib.lines
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerLine2D
import matplotlib
matplotlib.rcParams['font.sans-serif'] = ['FreeSans', ]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['axes.titlepad'] = 2
matplotlib.rcParams['axes.labelpad'] = 2


class SymHandler(HandlerLine2D):
    
    def create_artists(self, legend, orig_handle,xdescent, ydescent, 
                       width, height, fontsize, trans):
        xx= 0.6*height
        return super(SymHandler, self).create_artists(legend, orig_handle,xdescent, 
                     xx, width, height, fontsize, trans)


def intra_distances_worker(fin, taxonomy):
    df = pandas.read_csv(fin, index_col=0, sep='\t')
    species = species_df.loc[species_df.taxonomy == taxonomy, :].species_name
    df = df.loc[species, :]
    df = df.loc[:, species]
    matrix = df.values
    triu_indices = numpy.triu_indices(matrix.shape[0], 1, matrix.shape[1])
    pair_indices = zip(triu_indices[0], triu_indices[1])
    values = []
    for pair in pair_indices:
        values.append(matrix[pair[0], pair[1]])
    kernel = gaussian_kde(values, bw_method=0.25)
    X = numpy.linspace(0,1,S)
    Y = kernel.pdf(X) / sum(kernel.pdf(X)) * (S/100.)
    return X, Y


def inter_distances_worker(fin, taxonomy1, taxonomy2):
    df = pandas.read_csv(fin, index_col=0, sep='\t')
    species1 = species_df.loc[species_df.taxonomy == taxonomy1, :].species_name
    species2 = species_df.loc[species_df.taxonomy == taxonomy2, :].species_name  
    df = df.loc[species1, :]
    df = df.loc[:, species2]
    values = df.values.flatten()
    kernel = gaussian_kde(values, bw_method=0.25)
    X = numpy.linspace(0,1,S)
    Y = kernel.pdf(X) / sum(kernel.pdf(X)) * (S/100.)
    return X, Y


def worker(fin, fout):
    
    # Get intra-distances (distances between the members of each taxonomy)
    M1 = numpy.zeros((S,4)) # a 1000x4 matrix to store the KDE frequencies
    colnames = ['X']
    M1[:,0] = numpy.linspace(0,1,S) # first column are the x values
    for i, taxonomy in enumerate(taxonomy_colors.keys()):
        X, Y = intra_distances_worker(fin, taxonomy)
        M1[:,i+1] = Y
        colnames.append(taxonomy)
    M1 = pandas.DataFrame(M1, columns=colnames) # the matrix of intra-distances densities
    
    # Get inter-distances (distances between the members of each taxonomic class
    # and those of the other classes)
    M2 = numpy.zeros((S,4)) # a 1000x4 matrix to store the KDE frequencies
    colnames = ['X']
    M2[:,0] = numpy.linspace(0,1,S)
    for i, pair in enumerate(taxonomy_pairs):
        X, Y = inter_distances_worker(fin, pair[0], pair[1])
        colnames.append(','.join(pair))
        M2[:,i+1] = Y
    M2 = pandas.DataFrame(M2, columns=colnames) # the matrix of intra-distances densities
    
    # Make the plots
    fig = plt.figure(1, figsize=(8,6))
    grids = GridSpec(2,1, hspace=0.3) # two rows, one column
    xlabel1 = 'Intra-distances'
    xlabel2 = 'Inter-distances'
    ylabel = 'Frequency\n'

    
    # Plot for intra-distances
    ax = fig.add_subplot(grids[0,0])
    patches = []
    for taxonomy in M1.columns.values[1:]:
        ax.plot(M1.loc[:,'X'], M1.loc[:,taxonomy], linewidth=2, 
                c=taxonomy_colors[taxonomy], label=taxonomy.capitalize())
        p = Line2D([0], [0], marker='o', color=taxonomy_colors[taxonomy], 
                   markerfacecolor=taxonomy_colors[taxonomy], 
                   markersize=10, label=taxonomy.capitalize())
        patches.append(p)
    
    ax.grid(True, linestyle='--', alpha=0.25, linewidth=0.75)
    ax.legend(loc=9, handles=patches, ncol=3, fancybox=True, framealpha=0.5, 
              facecolor='#cccccc', borderpad=0, handleheight=2, frameon=False,
              handlelength=0.25, fontsize=10, columnspacing=1, handletextpad=0.5,
              handler_map={matplotlib.lines.Line2D: SymHandler()})
    ax.xaxis.set_tick_params(pad=1)
    ax.yaxis.set_tick_params(pad=1)
    ax.set_xticks(numpy.arange(0, 1.01, 0.2))
    ax.set_xlabel(xlabel1, fontsize=10)
    max_value = M1.iloc[:,1:].max().max()
    step = int(round(max_value+0.01, 2) *100)/4.
    up_limit = numpy.math.ceil(step)*4
    up_limit = up_limit/100.
    yticks = numpy.linspace(0, up_limit, 5)
    ax.set_yticks(yticks)
    ax.set_ylabel(ylabel, fontsize=10, labelpad=0.2)
    ax.tick_params(axis='both', which='major', labelsize=8)
    # Plot for inter-disctances
    ax = fig.add_subplot(grids[1,0])
    patches = []
    for pair_str in M2.columns.values[1:]:
        pair = tuple(pair_str.split(','))
        label = [p.capitalize() for p in pair]
        label = ' & '.join(label)
        ax.plot(M2.loc[:,'X'], M2.loc[:,pair_str], linewidth=2,
                c=taxonomy_pairs_colors[pair], label=label)
        p = Line2D([0], [0], marker='o', color=taxonomy_pairs_colors[pair],
                   markerfacecolor=taxonomy_pairs_colors[pair], markersize=10,
                   label=label)
        patches.append(p)   
    ax.grid(True, linestyle='--', alpha=0.25, linewidth=0.75)
    ax.legend(loc=9, handles=patches, ncol=3, fancybox=True, framealpha=0.5, 
              facecolor='#cccccc', borderpad=0, handleheight=2, frameon=False,
              handlelength=0.25, fontsize=10, columnspacing=1, handletextpad=0.5,
              handler_map={matplotlib.lines.Line2D: SymHandler()})
    ax.xaxis.set_tick_params(pad=1)
    ax.yaxis.set_tick_params(pad=1)    
    ax.set_xticks(numpy.arange(0, 1.01, 0.2))
    ax.set_xlabel(xlabel2, fontsize=10)
    max_value = M2.iloc[:,1:].max().max()
    step = int(round(max_value+0.01, 2) *100)/4.
    up_limit = numpy.math.ceil(step)*4
    up_limit = up_limit/100.
    yticks = numpy.linspace(0, up_limit, 5)
    ax.set_yticks(yticks)
    ax.set_ylabel(ylabel, fontsize=10, labelpad=0.2)
    ax.tick_params(axis='both', which='major', labelsize=8)  
    
    #plt.savefig(fout+'.tiff', dpi=900, format='png', bbox_inches = 'tight', 
    #            pad_inches = 0.05)
    plt.savefig(fout+'.png', dpi=600, format='png', bbox_inches = 'tight', 
                pad_inches = 0.05)
    plt.close()



if __name__ == '__main__':
    
    '''
    DESCRIPTION:
    
    This script creates the density plots for the inter- and intra-distances of
    taxonomic classes, for the four phylogenetic criteria. Given the calculated 
    values from the distance matrices, a KDE function is constructed. Then 
    1000 values from 0 to 1 (X values) are used to calculate the densities 
    (Y values). A figure with two subplots is created for each phylogenetic
    criterion. The first one depicts the KDEs of intra-distances and the second 
    one those of inter-distances for the three taxonomic kingdoms.
    '''

    species_df = pandas.read_csv('./files/species.tsv', sep='\t', index_col=0)
    species_dict = dict(zip(species_df.species_name, species_df.taxonomy))
    taxonomy_pairs = list(itertools.combinations(taxonomy_colors.keys(), r=2))
    taxonomy_pairs = [sorted(p) for p in taxonomy_pairs]
    S = 1000 # variable used as input in the distances KDEs, to calculate frequencies 
    
    filenames = {
        'rRNA':'./rRNAs_analysis/final_distance_matrix.tsv',
        'hsp40':'./HSPs_analysis/hsp40/final_distance_matrix.tsv',
        'hsp70':'./HSPs_analysis/hsp70/final_distance_matrix.tsv',
        'proteostasis_network':'./PN_analysis/final_distance_matrix.tsv'
    }
    
    main_dir = './comparison_of_clustergrams/'
    if not os.path.exists(main_dir):
        os.makedirs(main_dir)
    
    for criterion, fin in filenames.items():
        fout = main_dir+criterion+'_distribution_of_distances'
        worker(fin, fout)



