#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import pandas
from scipy.spatial import distance
from scipy.cluster import hierarchy
from sklearn.metrics import silhouette_score
from sklearn.metrics import silhouette_samples
from sklearn.metrics import homogeneity_score
import numpy
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerLine2D
import matplotlib
matplotlib.rcParams['font.sans-serif'] = ['FreeSans', ]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['axes.titlepad'] = 4


class SymHandler(HandlerLine2D):
    def create_artists(self, legend, orig_handle,xdescent, ydescent, width, 
                       height, fontsize, trans):
        xx= self.y_align*height
        return super(SymHandler, self).create_artists(legend, orig_handle,xdescent, 
                     xx, width, height, fontsize, trans)
SH = SymHandler()
SH.y_align = 0.45


def weighted_silhouette_score(matrix, ref_labels, labels):
    silhouettes = silhouette_samples(matrix, labels=labels, metric='precomputed')
    
    weights = {}
    for tax in list(set(ref_labels)):
        ratio = 1./(ref_labels.count(tax)*3)
        weights.update({tax:ratio})
    
    s = 0
    for i in range(len(ref_labels)):
        label = ref_labels[i]
        weight = weights[label]
        s = s + silhouettes[i]*weight
    
    return round(s,3)



if __name__ == '__main__':
    
    '''
    DESCRIPTION:

    This script creates the 4rth plot of tha manuscript, in order to compare 
    the silhouette and homogeneity scores between the 20 conserved mechanisms,
    PN and ribosomal sequences, when they are used as criteria for the separation
    of the 3 taxonomic kingdoms.
    '''
    
    species_df = pandas.read_csv('./files/species.tsv', sep='\t', index_col=0)
    taxonomy_mapping = dict(zip(species_df.species_name.tolist(), species_df.taxonomy.tolist()))
    main_dir = './other_terms_analysis/'
    other_terms = os.listdir(main_dir+'pathway_analysis_inputs/')
    
    data = []

    # Calculate Homogeneity and Silhouette scores for the 20 conserved mechanisms
    # with and without PN components
    for other_term in other_terms:
        
        tmp_data = [other_term]
        
        tmp_dir = './other_terms_analysis/semantic_analysis/'+other_term+'/'
        matrix_df = pandas.read_csv(tmp_dir+'final_distance_matrix.tsv', sep='\t', index_col=0)
        condensed_matrix = distance.squareform(matrix_df.values, force='to_vector')
        Z = hierarchy.linkage(condensed_matrix, method='ward')
        labels = list(hierarchy.fcluster(Z, criterion='maxclust',t=3))
        silhouette = round(silhouette_score(matrix_df.values, labels=labels, metric='precomputed'), 2)
        ref_labels = [taxonomy_mapping[s] for s in matrix_df.index]
        weighted_silhouette = weighted_silhouette_score(matrix_df.values, ref_labels, labels=labels)    
        homogeneity = round(homogeneity_score(ref_labels, labels),2)
        tmp_data.append(homogeneity)
        tmp_data.append(silhouette)
        tmp_data.append(weighted_silhouette)
        
        tmp_dir = './other_terms_analysis/semantic_analysis_without_PN_components/'+other_term+'/'
        matrix_df = pandas.read_csv(tmp_dir+'final_distance_matrix.tsv', sep='\t', index_col=0)
        condensed_matrix = distance.squareform(matrix_df.values, force='to_vector')
        Z = hierarchy.linkage(condensed_matrix, method='ward')
        labels = list(hierarchy.fcluster(Z, criterion='maxclust',t=3))
        silhouette = round(silhouette_score(matrix_df.values, labels=labels, metric='precomputed'), 2)
        ref_labels = [taxonomy_mapping[s] for s in matrix_df.index]
        weighted_silhouette = weighted_silhouette_score(matrix_df.values, ref_labels, labels=labels)    
        homogeneity = round(homogeneity_score(ref_labels, labels),2)
        tmp_data.append(homogeneity)
        tmp_data.append(silhouette)
        tmp_data.append(weighted_silhouette)
        
        data.append(tmp_data)
        
    columns = ['term', 'homogeneity_score1', 'silhouette_score1', 'weighted_silhouette_score1',
               'homogeneity_score2', 'silhouette_score2', 'weighted_silhouette_score2'] 
    results_df = pandas.DataFrame(data, columns=columns)
    results_df.mechanism = results_df.term.apply(lambda x: x.replace('_', ' ').capitalize())
    results_df.sort_values(by='homogeneity_score1', ascending=False, inplace=True)
    results_df.reset_index(drop=True, inplace=True)
    
    
    # Calculate Homogeneity and Silhouette scores for the ribosomal sequences
    matrix_df = pandas.read_csv('./rRNAs_analysis/final_distance_matrix.tsv', sep='\t', 
                                index_col=0)
    condensed_matrix = distance.squareform(matrix_df.values, force='to_vector')
    Z = hierarchy.linkage(condensed_matrix, method='ward')
    labels = list(hierarchy.fcluster(Z, criterion='maxclust',t=3))
    ref_labels = [taxonomy_mapping[s] for s in matrix_df.index]
    rRNA_homogeneity = round(homogeneity_score(ref_labels, labels), 3)
    rRNA_weighted_silhouette = weighted_silhouette_score(matrix_df.values, ref_labels,
                                                         labels=labels)  
    
    # Calculate Homogeneity and Silhouette scores for the PN profiles   
    matrix_df = pandas.read_csv('./PN_analysis/final_distance_matrix.tsv', sep='\t',
                                index_col=0)
    condensed_matrix = distance.squareform(matrix_df.values, force='to_vector')
    Z = hierarchy.linkage(condensed_matrix, method='ward')
    labels = list(hierarchy.fcluster(Z, criterion='maxclust',t=3))
    ref_labels = [taxonomy_mapping[s] for s in matrix_df.index]
    proteostasis_homogeneity = round(homogeneity_score(ref_labels, labels), 3)
    proteostasis_weighted_silhouette = weighted_silhouette_score(matrix_df.values, 
                                                                 ref_labels, 
                                                                 labels=labels)  
    
    # Create the two-panels bar plot
    fig = plt.figure(1, figsize=(7,6), frameon=False)
    grids = GridSpec(1,2, hspace=0, wspace=0.1)
    
    # Add bar plot for the Homogeneity score
    ax = fig.add_subplot(grids[0,0])
    step = 1
    barwidth = 0.5
    x0 = 2
    x1 = [x0+i*step for i in range(len(results_df))]
    x2 = [x0+i*step+barwidth for i in range(len(results_df))]
    ylabels = []
    s = {'Atp': 'ATP', 'dna':'DNA', 'Dna': 'DNA', 'Trna': 'tRNA'}
    for term in results_df.term.tolist():
        term = term.capitalize()
        term = term.replace('_', ' ')
        for k,v in s.items():
            term = term.replace(k,v)
        ylabels.append(term)
    i_bars = ax.barh([0, 1], [rRNA_homogeneity, proteostasis_homogeneity],
                     height=barwidth, color='#650E16', linewidth=2, edgecolor='#A55960')  
    ax.plot(results_df.homogeneity_score2, x1, linewidth=2, color='#650E16', zorder=2)
    ax.scatter(results_df.homogeneity_score2, x1, marker='o', s=16, c='#650E16', zorder=3)
    hs_bars = ax.barh(x1, results_df.homogeneity_score1, height=barwidth, color='#A55960', zorder=1)
    
    y_ticks = [-0.03, 0.87] + [x0+i*step-0.03 for i in range(len(results_df))]
    ylabels.insert(0, 'Ribosomal Sequences')
    ylabels.insert(1, 'Proteostasis')
    plt.yticks(y_ticks, ylabels, rotation=0, fontsize=12)
    plt.xticks(numpy.arange(0, 1.05, 0.2), fontsize=12)
    ax.invert_yaxis()
    ax.set_axisbelow(True)
    ax.xaxis.grid(True, linestyle='--', alpha=0.5, linewidth=0.7)
    plt.xlim([0,1.05])
    
    
    # Add bar plot for the Silhouette score
    ax = fig.add_subplot(grids[0,1])
    step = 1
    barwidth = 0.5
    x1 = [x0+i*step for i in range(len(results_df))]
    x2 = [x0+i*step+barwidth for i in range(len(results_df))]
    i_bars = ax.barh([0, 1], [rRNA_weighted_silhouette, proteostasis_weighted_silhouette], 
                     height=barwidth, color='#0e655d', linewidth=2, edgecolor='#59a59e')  
    hs_bars = ax.barh(x1, results_df.weighted_silhouette_score1, height=barwidth, color='#59a59e')
    ax.plot(results_df.weighted_silhouette_score2, x1, linewidth=2,  color='#0e655d', zorder=2)
    ax.scatter(results_df.weighted_silhouette_score2, x1, marker='o', s=16, c='#0e655d', zorder=3)
    plt.yticks([], [])
    plt.xticks(numpy.arange(0, 1.05, 0.2), fontsize=12)
    ax.invert_yaxis()
    ax.set_axisbelow(True)
    ax.xaxis.grid(True, linestyle='--', alpha=0.5, linewidth=0.7)
    plt.xlim([0,1.05])
    
    # Legend
    patches = []
    patches_labels = []
    p = Line2D([0], [0], marker='o', color='#A55960', markerfacecolor='#A55960',
               markersize=14, label='Homogeneity Score (entire profile)')
    patches.append(p)
    patches_labels.append('Homogeneity Score (entire profile)')
    p = Line2D([0], [0], marker='o', color='#650E16', markerfacecolor='#650E16',
               markersize=14, label='Homogeneity Score (without PN components)')
    patches.append(p)
    patches_labels.append('Homogeneity Score (without PN components)')
    p = Line2D([0], [0], marker='o', color='#59a59e', markerfacecolor='#59a59e',
               markersize=14, label='Silhouette Score (entire profile)')
    patches.append(p)
    patches_labels.append('Silhouette Score (entire profile)')
    p = Line2D([0], [0], marker='o', color='#0e655d', markerfacecolor='#0e655d',
               markersize=14, label='Silhouette Score (without PN components)')
    patches.append(p)
    patches_labels.append('Silhouette Score (without PN components)')
    plt.legend(handles=patches, 
               labels=patches_labels, 
               bbox_to_anchor=(0.4, -0.55, 0.55, 0.5),
               fontsize=14, ncol=2, fancybox=False, framealpha=0, 
               handlelength=0.5, handletextpad=0.8,
               handleheight=2,
               labelspacing=0.8,
               handler_map={matplotlib.lines.Line2D: SH})
         
    main_dir = './other_terms_analysis/'
    plt.savefig(main_dir+'bars_plot_figure_4.png', dpi=900, format='png', bbox_inches='tight', pad_inches=0.01)
    #plt.savefig(main_dir+'bars_plot_01.tiff', dpi=800, format='tiff', bbox_inches='tight', pad_inches=0.01)
    plt.close()

