#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy
import pandas   
from sklearn.manifold import MDS
from matplotlib import pyplot as plt
import sys
sys.path.append('../')
from core_functions import taxonomy_colors



if __name__ == '__main__':
    
    '''
    DESCRIPTION:

    This script creates the second part of the 3rd figure of the manuscript. 
    It depicts the species of each kingdom on a 2-dimensional space, based on
    their semantic distances given 1) the whole PN profiles, 2) the conserved
    PN profile and 3) the differential PN profile.
    '''
    
    filenames = {
        'whole':'./PN_analysis/final_distance_matrix.tsv',
        'common':'./PN_components/common_and_differential/common_components_distance_matrix.tsv',
        'diff':'./PN_components/common_and_differential/differential_components_distance_matrix.tsv',
    }    
    titles = {
        'whole': 'Whole Proteostasis Network',
        'common': 'Conserved Core Network',
        'diff': 'Differential Network'
    }
    
    df = pandas.read_csv('./files/species.tsv', sep='\t', index_col=0)
    species_taxonomy_mapping = dict(zip(df.abbreviation, df.taxonomy))
    species_names_mapping = dict(zip(df.species_name, df.abbreviation))
    
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(8,2.5), frameon=False)
    
    for i, key in enumerate(['whole', 'common', 'diff']):
        
        matrix = pandas.read_csv(filenames[key], sep='\t', index_col=0)
        matrix.index = [species_names_mapping[s] for s in matrix.index]
        mds_obj = MDS(dissimilarity='precomputed', n_components=2, n_jobs=4, 
                      random_state=100, eps=1e-6, n_init=10)
        coords = mds_obj.fit_transform(matrix.values)
        coords = pandas.DataFrame(coords, index=matrix.index)
        print(key, coords.shape)
        colors = [taxonomy_colors[species_taxonomy_mapping[i]] for i in matrix.index]  
    
        lim = 0.7
        ticks = numpy.arange(-0.8, 0.8, 0.2)
    
        axes[i].scatter(coords.loc[:,0], coords.loc[:,1], c=colors, s=5, 
                        linewidths=.0, edgecolor='black',  alpha=0.85)
        axes[i].set_xlim([-lim,lim])
        axes[i].set_ylim([-lim,lim])
        axes[i].set_xticks(ticks, size=0)
        axes[i].set_yticks(ticks, size=0)
        axes[i].set_title(titles[key], size=12)
        axes[i].set_xticks(ticks, size=6)
        axes[i].set_xlabel('MDS Component 1', size=10, labelpad=1)
        axes[i].set_yticks(ticks, size=6)
        axes[i].set_ylabel('MDS Component 2', size=10, labelpad=1)
        axes[i].tick_params(labelsize=6, pad=0.75, width=0.75, length=2, direction='out')
        #plt.setp(ax.get_xticklabels(), visible=False)
        #plt.setp(ax.get_yticklabels(), visible=False)
        axes[i].grid(True, linestyle='--', alpha=0.5, linewidth=0.3)
    fig.tight_layout()
    plt.savefig('./PN_components/common_and_differential/scatter_plots.png', 
                dpi=800, format='png')
    #plt.savefig('./PN_components/common_and_differential/scatter_plots.tiff', 
                #dpi=600, format='tiff')    
    

 
    
    
    