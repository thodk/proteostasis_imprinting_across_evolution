# -*- coding: utf-8 -*-

import os
import numpy
import pandas
import itertools
from matplotlib import pyplot as plt
import matplotlib.lines
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerLine2D
from scipy.stats import gaussian_kde
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


def intra_distances_worker(filename, category):
    df = pandas.read_csv(filename, index_col=0, sep='\t')
    species = species_df.loc[species_df.taxonomy == category, :].species_name
    df = df.loc[species, :]
    df = df.loc[:, species]
    print df.shape
    
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


def inter_distances_worker(filename, category1, category2):
    df = pandas.read_csv(filename, index_col=0, sep='\t')
    species1 = species_df.loc[species_df.taxonomy == category1, :].species_name
    species2 = species_df.loc[species_df.taxonomy == category2, :].species_name  
    df = df.loc[species1, :]
    df = df.loc[:, species2]
    print category1,  category2, df.shape
    
    values = df.values.flatten()
        
    kernel = gaussian_kde(values, bw_method=0.25)
    X = numpy.linspace(0,1,S)
    Y = kernel.pdf(X) / sum(kernel.pdf(X)) * (S/100.)
    return X, Y


def worker(filename, output):
    
    # Intra Distances
    M1 = numpy.zeros((S,4))
    colnames = ['X']
    M1[:,0] = numpy.linspace(0,1,S)
    for i, category in enumerate(colors_dict.keys()):
        X, Y = intra_distances_worker(filename, category)
        M1[:,i+1] = Y
        colnames.append(category)
    M1 = pandas.DataFrame(M1, columns=colnames)
    
    # Inter Distances
    M2 = numpy.zeros((S,4))
    colnames = ['X']
    M2[:,0] = numpy.linspace(0,1,S)
    for i, pair in enumerate(categories_pairs_dict.keys()):
        X, Y = inter_distances_worker(filename, pair[0], pair[1])
        colnames.append(','.join(pair))
        M2[:,i+1] = Y
    M2 = pandas.DataFrame(M2, columns=colnames)
    
    # PLOT
    fig = plt.figure(1, figsize=(4.5,4.5))
    grids = GridSpec(2,1, hspace=0.25)
    xlabel1 = 'Intra-distances'
    xlabel2 = 'Inter-distances'   
    xlabel_size = 8
    xticks_size = 5
    ylabel = 'Frequency\n'
    ylabel_size = 8
    yticks_size = 5
    
    # Intra - distances plot
    ax = fig.add_subplot(grids[0,0])
    
    patches = []
    for category in M1.columns.values[1:]:
        color = colors_dict[category]
        ax.plot(M1.loc[:,'X'], M1.loc[:,category], linewidth=1, c=color, 
                label=category.capitalize())
        p = Line2D([0], [0], marker='o', color=color, markerfacecolor=color, 
                   markersize=5, label=category.capitalize())
        patches.append(p)
    
    ax.grid(True, linestyle='--', alpha=0.75, linewidth=0.75)
    ax.legend(loc=9, handles=patches, ncol=3, fancybox=True, framealpha=0.5, 
              facecolor='#cccccc', borderpad=0, handleheight=2, frameon=False,
              handlelength=0.5, fontsize=6, columnspacing=2, handletextpad=1,
              handler_map={matplotlib.lines.Line2D: SymHandler()})
              
    ax.xaxis.set_tick_params(pad=1)
    ax.yaxis.set_tick_params(pad=1)
    
    plt.xticks(fontsize=xticks_size)
    plt.xlabel(xlabel1, fontsize=xlabel_size)
    
    max_value = M1.iloc[:,1:].max().max()
    up_limit = int(round(max_value+0.01, 2) *100)
    step = up_limit/4.
    up_limit = numpy.math.ceil(step)*4
    up_limit = up_limit/100.
    yticks = numpy.linspace(0, up_limit, 5)
    plt.yticks(yticks, fontsize=yticks_size)
    plt.ylabel(ylabel, fontsize=ylabel_size, labelpad=0.2)
    
    # Inter - distances plot
    ax = fig.add_subplot(grids[1,0])
    
    patches = []
    for pair in M2.columns.values[1:]:
        color = categories_pairs_dict[tuple(pair.split(','))]
        label = [i.capitalize() for i in pair.split(',')]
        label = ' & '.join(label)
        ax.plot(M2.loc[:,'X'], M2.loc[:,pair], linewidth=1, c=color, label=label)
        p = Line2D([0], [0], marker='o', color=color, markerfacecolor=color, 
                   markersize=5, label=label)
        patches.append(p)
        
    ax.grid(True, linestyle='--', alpha=0.75, linewidth=0.75)
    ax.legend(loc=9, handles=patches, ncol=3, fancybox=True, framealpha=0.5, 
              facecolor='#cccccc', borderpad=0, handleheight=2, handletextpad=1,
              handlelength=0.5, fontsize=6, columnspacing=2, frameon=False,
              handler_map={matplotlib.lines.Line2D: SymHandler()})
              
    ax.xaxis.set_tick_params(pad=1)
    ax.yaxis.set_tick_params(pad=1)    
    
    plt.xticks(fontsize=xticks_size)
    plt.xlabel(xlabel2, fontsize=xlabel_size)
    
    max_value = M2.iloc[:,1:].max().max()
    up_limit = int(round(max_value+0.01, 2) *100)
    step = up_limit/4.
    up_limit = numpy.math.ceil(step)*4
    up_limit = up_limit/100.
    yticks = numpy.linspace(0, up_limit, 5)
    plt.yticks(yticks, fontsize=yticks_size)
    plt.ylabel(ylabel, fontsize=ylabel_size, labelpad=0.2)
    
    #plt.savefig(output+'.tiff', dpi=800,  format='tiff', bbox_inches = 'tight',
                #pad_inches = 0.05)
    plt.savefig(output+'.png', dpi=800, format='png', bbox_inches = 'tight', 
                pad_inches = 0.05)

    plt.close()



if __name__ == '__main__':
    species_df = pandas.read_csv('./files/species.tsv', sep='\t', index_col=0)
    species_dict = dict(zip(species_df.species_name, species_df.taxonomy))
    
    colors_df = pandas.read_csv('./files/colors.csv', index_col=0)
    colors_dict = dict(zip(colors_df.taxonomy, colors_df.color))
    categories_pairs = [p for p in itertools.combinations(colors_dict.keys(),2)]
    
    categories_pairs_dict = {('archaea', 'bacteria'): '#E95C20',
                             ('archaea', 'eukaryotes'): '#006747',
                             ('bacteria', 'eukaryotes'): '#4F2C1D'}
    
    S = 1000
    
    filenames = {
        'rRNA':'./rRNAs_analysis/final_distance_matrix.tsv',
        'hsp40':'./HSPs_analysis/hsp40/final_distance_matrix.tsv',
        'hsp70':'./HSPs_analysis/hsp70/final_distance_matrix.tsv',
        'proteostasis_network':'./PN_analysis/final_distance_matrix.tsv'
    }
    
    main_dir = './comparison_of_phylogenetic_trees/'
    if not os.path.exists(main_dir):
        os.makedirs(main_dir)
    
    for method, filename in filenames.items():
        worker(filename, main_dir+method+'_distances_distributions')



