#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas
import numpy
import sys
import os
from sklearn.decomposition import NMF
from sklearn.metrics import homogeneity_score
from sklearn.metrics import silhouette_score
from sklearn.cluster import KMeans
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib import pyplot as plt
import warnings
warnings.filterwarnings("ignore")
sys.path.append('../')
from core_functions import taxonomy_colors
import heatmap


if __name__ == '__main__':
    
    '''
    DESCRIPTION:

    The non-negative matrix factorization is applied on the association_matrix
    (derived in the previous step), using 3 of components. The whole workflow
    is described in the section 2.5.2 of the manuscript. The outputs are two 
    matrices, which indicate the relevance of species and GO-BP clusters to 
    the derived NMF-based components. Also these matrices are illustrated in 
    the respective heatmaps, using the heatmap.py script from the upper directory.
    '''

    if not os.path.exists('./PN_components/nmf/'):
        os.makedirs('./PN_components/nmf/')
    
    df = pandas.read_csv('./files/species.tsv', sep='\t', index_col=0)
    species_taxonomy_mapping = dict(zip(df.abbreviation, df.taxonomy))
    species_names_mapping = dict(zip(df.abbreviation, df.species_name))
    
    
    association_df = pandas.read_csv('./PN_components/association_matrix.tsv',
                                 sep='\t', index_col=0)
    association_df = association_df.transpose()
    association_matrix = association_df.values
    true_labels = [species_taxonomy_mapping[s] for s in association_df.index]
    
    # Run Non Negative Matrix Decomposition
    nmf_results = {}
    for k in numpy.arange(1,16):
        nmf_model = NMF(n_components=k, tol=1e-4, max_iter=1600)
        W = nmf_model.fit_transform(association_matrix)
        H = nmf_model.components_
        error = nmf_model.reconstruction_err_
        nmf_results.update({k:(W,H,error)})
    
    # Caclulate homogeneity_score & silhouette_score
    scores = []
    for k in numpy.arange(1,16):
        kmeans_model = KMeans(n_clusters=3, random_state=1234)
        predictions = kmeans_model.fit_predict(nmf_results[k][0])
        hm = homogeneity_score(labels_true=true_labels, labels_pred=predictions).round(2)
        sil = silhouette_score(nmf_results[k][0], labels=predictions).round(2)
        scores.append([k, hm, sil])
    scores = numpy.array(scores)
    
    #
    # PLOT 1: homogeneity_score & silhouette_score
    #
    fig, ax = plt.subplots(1, 1, figsize=(8,4))
    ax.scatter(scores[:,0], scores[:,1], marker='o', c='#7f121c')
    ax.plot(scores[:,0], scores[:,1], label='Homogeneity', linewidth=2, color='#7f121c')
    ax.scatter(scores[:,0], scores[:,2], marker='o', c='#127f75')
    ax.plot(scores[:,0], scores[:,2], label='Silhouette', linewidth=2, color='#127f75')
    ax.set_xticks(numpy.arange(1,16))
    ax.set_xlabel('Number of Components', size=16)
    ax.set_yticks(numpy.arange(0, 1.01, 0.1))
    ax.set_ylabel('Score', size=16)
    ax.legend(frameon=False, loc=4)
    ax.grid(True, alpha=0.25)
    plt.tight_layout()
    plt.savefig('./PN_components/nmf/NMF_K_selection.png', dpi=600, format='png')
    #plt.savefig('./PN_components/nmf/NMF_K_selection.tiff', dpi=600, format='tiff')    
    
    # Components matrix
    H = nmf_results[3][1]
    H = H / numpy.linalg.norm(H, ord=numpy.inf, axis=1).reshape(-1,1)
    HT = H.transpose()
    HT_df = pandas.DataFrame(HT, index=association_df.columns, columns=range(1, HT.shape[1]+1))
    HT_df.to_csv('./PN_components/nmf/components_matrix.tsv', sep='\t')
    HT[HT == 0] = -1e-8    
    
    common_terms = HT_df.loc[(HT_df > 0).all(axis=1),].index.tolist()

    #
    # PLOT 2: components matrix (H)
    #
    row_clustering_data={'data':HT, 'method':'ward', 'metric':'euclidean', 'ratio':0.15}
    col_clustering_data=None
    
    legend_ticks = numpy.arange(0, 1.01, 0.1)
    legend_labels = ["{:.1f}".format(f) for f in legend_ticks]
    legend_title = 'Association\nCoefficient\n'
    legend_data = {'x':0.0, 'y':0.03, 'w':0.06, 'h':0.7, 'labels':legend_labels, 
                   'labels_size':12, 'cbar_kws':{'ticks':legend_ticks}, 
                   'title':legend_title, 'title_size':16}
    
    x_axis_data = {'labels': ['Component '+str(i) for i in range(1, HT.shape[1]+1)], 
                   'fontdict':{'size':13.0, 'rotation':0}}
    y_labels = [s.capitalize() for s in association_df.columns]
    y_specific_labels_format = {}
    for i, label in enumerate(y_labels):
        if label.lower() in common_terms:
            tmp_dict = {'color': '#000000', 'weight':900}
        else:
            tmp_dict = {'color': '#000000', 'weight':300}
        y_specific_labels_format.update({label:tmp_dict})
    y_axis_data = {'labels': y_labels, 'fontdict':{'size':14}, 
                   'specific_labels_format':y_specific_labels_format}
       
    # Define some necessary paramters for the 'heatmap' module
    heatmap_data={'data':HT, 'type':'features', 'x_ratio':0.8, 'y_ratio':1}  
    #cmap = LinearSegmentedColormap.from_list("my_colormap", ('#cdcdcd', '#000000'), N=100)
    cmap = LinearSegmentedColormap.from_list("my_colormap", ('#e0eaf1', '#537895'), N=100)
    cmap.set_under('#ffffff')
    
    # Call the function to initiate the clustergram construction
    c = heatmap.Clustergram(heatmap_data, figsize=(6,15), cmap=cmap, 
                            y_axis_data=y_axis_data,
                            x_axis_data=x_axis_data,
                            row_clustering_data=row_clustering_data,
                            col_clustering_data=None,
                            row_colors_data = None,
                            col_colors_data = None,
                            colors_legend_data = None, vmin=-1e-8,
                            norm=Normalize(0,1), legend_data=legend_data, 
                            linecolor='#e0e0e0', linewidth=0.01)
    
    # construction
    c.construction()
    c.set_coords()
    c.set_labels()
    c.clustergrid.ax_heatmap.set_frame_on(False)
    c.clustergrid.ax_heatmap.tick_params(pad=1, width=0.05, length=0.3, direction='out')
    c.clustergrid.cax.set_title(label='Normalized\nCoefficient', pad=10, size=16)
    c.clustergrid.savefig('./PN_components/nmf/components.png', dpi=600, format='png', 
                          pad_inches=0.02)
    #c.clustergrid.savefig('./PN_components/nmf/components.tiff', dpi=600, format='tiff', 
                          #pad_inches=0.02)

    # Transformed matrix
    W = nmf_results[3][0]
    W = W / numpy.linalg.norm(W, ord=numpy.inf, axis=1).reshape(-1,1)
    WT = W.transpose()
    W_df = pandas.DataFrame(W, index=association_df.index, columns=range(1, W.shape[1]+1))
    W_df.to_csv('./PN_components/nmf/transformed_matrix.tsv', sep='\t')
    WT[WT == 0] = -1e-8

    #
    # PLOT 3: tranformed matrix (W)
    #
    col_clustering_data={'data':W, 'method':'ward', 'metric':'euclidean', 'ratio':0.3}
    row_clustering_data=None
    
    legend_ticks = numpy.arange(0, 1.01, 0.2)
    legend_labels = ["{:.1f}".format(f) for f in legend_ticks]
    legend_title = 'Normalized\nCoefficient'
    legend_data = {'x':0.11, 'y':0.25, 'w':0.02, 'h':0.4, 'labels':legend_labels, 
                   'labels_size':8, 'cbar_kws':{'ticks':legend_ticks}, 
                   'title':legend_title, 'title_size':12}
    
    y_axis_data = {'labels': ['Component '+str(i) for i in range(1, WT.shape[0]+1)], 'fontdict':{'size':12}}
    
    x_labels = [i for i in list(association_df.index)]
    x_labels_for_plot = [species_names_mapping[i] for i in x_labels]
    x_labels_dict = dict(zip(x_labels, x_labels_for_plot))
    x_specific_labels_format = {}
    for label in x_labels:
        tmp_dict = {'color': taxonomy_colors[species_taxonomy_mapping[label]],
                    'weight':600}
        x_specific_labels_format.update({x_labels_dict[label]:tmp_dict})
    x_axis_data = {'labels': x_labels_for_plot,
                   'specific_labels_format':x_specific_labels_format, 
                   'fontdict':{'size':1, 'rotation':90}}
    
    colors = [taxonomy_colors[species_taxonomy_mapping[i]] for i in x_labels]    
    colors_data = {'colors':colors, 'ratio':0.02}
    colors_legend_data = {}
    colors_legend_data.update({'patches': [[], []]})
    for key, value in taxonomy_colors.items():
        colors_legend_data['patches'][0].append(key.capitalize())
        p = Line2D([0], [0], marker='o', color=value, markerfacecolor=value,
                   markersize=10, label=key.capitalize())
        colors_legend_data['patches'][1].append(p)
    
    colors_legend_data.update({'title':'Taxonomy'})
    colors_legend_data.update({'bbox_anchor':(0.2,0.2), 'fontsize':10, 
                               'handlelength':0.75, 'handletextpad':0.25,
                               'handleheight':1.5, 'title_size':12})
    
    # Define some necessary paramters for the 'heatmap' module
    heatmap_data={'data':WT, 'type':'features', 'x_ratio':0.8, 'y_ratio':1}  
    #cmap = LinearSegmentedColormap.from_list("my_colormap", ('#cdcdcd', '#000000'), N=100)
    cmap = LinearSegmentedColormap.from_list("my_colormap", ('#f1e8e0', '#957053'), N=100)
    cmap.set_under('#ffffff')
    
    # Call the function to initiate the clustergram construction
    c = heatmap.Clustergram(heatmap_data, figsize=(8,3), cmap=cmap, 
                            y_axis_data=y_axis_data,
                            x_axis_data=x_axis_data,
                            row_clustering_data=row_clustering_data,
                            col_clustering_data=col_clustering_data,
                            row_colors_data = None,
                            col_colors_data = colors_data,
                            colors_legend_data = colors_legend_data, vmin=-1e-8,
                            norm=Normalize(0,1), legend_data=legend_data, 
                            linecolor='#e0e0e0', linewidth=0)
    
    # construction
    c.construction()
    c.set_coords()
    c.set_labels()
    c.clustergrid.ax_heatmap.set_frame_on(False)
    c.clustergrid.ax_heatmap.tick_params(pad=1, width=0.05, length=0.3, direction='out')
    c.clustergrid.cax.set_title(label='Normalized\nCoefficient', pad=8, size=12)
    c.clustergrid.savefig('./PN_components/nmf/transformed_matrix.png', dpi=600, 
                          format='png', pad_inches=0.02)
    #c.clustergrid.savefig('./PN_components/nmf/transformed_matrix.tiff', dpi=600, 
                          #format='tiff', pad_inches=0.02)