#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas
import sys
import numpy
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib import pyplot as plt
sys.path.append('../')
import heatmap
from core_functions import taxonomy_colors



if __name__ == '__main__':
    
    '''
    DESCRIPTION:

    This script creates a part of the 3rd figure of the manuscript. This subplot
    illustrates in a heatmap the association of each species with the derived 
    semantic clusters. 
    '''
    
    # Load the names of species and a mapping for their taxonomy
    df = pandas.read_csv('./files/species.tsv', sep='\t', index_col=0)
    species_taxonomy_mapping = dict(zip(df.abbreviation, df.taxonomy))
    species_names_mapping = dict(zip(df.abbreviation, df.species_name))
    
    # Load the association matrix
    association_df = pandas.read_csv('./PN_components/association_matrix.tsv',
                                     sep='\t', index_col=0)
    association_matrix = association_df.values
    association_matrix[association_matrix==0] = -1e-08
    # Row and column clustering object
    row_clustering_data={'data':association_matrix, 'method':'ward', 
                     'metric':'euclidean', 'ratio':0.0}
    col_clustering_data={'data':association_matrix.transpose(), 'method':'ward', 
                         'metric':'euclidean', 'ratio':0.1}
    
    # x and y labels for the heatmap
    y_labels = [i for i in list(association_df.index)]
    y_labels_for_plot = [i.capitalize() for i in y_labels]
    y_labels_dict = dict(zip(y_labels, y_labels_for_plot))
    
    x_labels = [i for i in list(association_df.columns.values)]
    x_labels_for_plot = [species_names_mapping[i] for i in x_labels]
    x_labels_dict = dict(zip(x_labels, x_labels_for_plot))
    
    colors = [taxonomy_colors[species_taxonomy_mapping[i]] for i in x_labels]    
    colors_data = {'colors':colors, 'ratio':0.02}
    colors_legend_data = {}
    colors_legend_data.update({'patches': [[], []]})
    for key, value in taxonomy_colors.items():
        colors_legend_data['patches'][0].append(key.capitalize())
        p = Line2D([0], [0], marker='o', color=value, markerfacecolor=value,
                   markersize=14, label=key.capitalize())
        colors_legend_data['patches'][1].append(p)
    
    colors_legend_data.update({'title':'Taxonomy'})
    colors_legend_data.update({'bbox_anchor':(0.205,0.26), 'fontsize':12, 
                               'handlelength':1, 'handletextpad':0.25,
                               'handleheight':1.5, 'title_size':14})
    

    legend_ticks = numpy.arange(0, 4.5, 0.5)    
    legend_labels = ["{:.1f}".format(f) for f in legend_ticks]

    legend_title = 'Association\nScore'
    legend_data = {'x':0.09, 'y':0.33, 'w':0.03, 'h':0.25, 'labels':legend_labels, 
                   'labels_size':12, 'cbar_kws':{'ticks':legend_ticks}, 
                   'title':legend_title, 'title_size':14}
    
    x_specific_labels_format = {}
    for label in x_labels:
        tmp_dict = {'color': taxonomy_colors[species_taxonomy_mapping[label]],
                    'weight':600}
        x_specific_labels_format.update({x_labels_dict[label]:tmp_dict})
    x_axis_data = {'labels': x_labels_for_plot,
                   'specific_labels_format':x_specific_labels_format, 
                   'fontdict':{'size':1, 'rotation':90}}
    
    y_specific_labels_format = {}
    for i, label in enumerate(y_labels):
            tmp_dict = {'color': '#000000', 'weight':300}
            y_specific_labels_format.update({y_labels_dict[label]:tmp_dict})
    y_axis_data = {'labels': y_labels_for_plot,
                   'specific_labels_format':y_specific_labels_format,
                   'fontdict':{'size':8}}


    heatmap_data={'data':association_matrix, 'type':'features', 'x_ratio':1.5, 'y_ratio':1.2}
    cmap = LinearSegmentedColormap.from_list("my_colormap", ('#cdcdcd', '#000000'), N=80)
    cmap.set_under('#ffffff')
    c = heatmap.Clustergram(heatmap_data, figsize=(8,8), cmap=cmap, 
                            y_axis_data=y_axis_data,
                            x_axis_data=x_axis_data,
                            row_clustering_data=row_clustering_data,
                            col_clustering_data=col_clustering_data,
                            row_colors_data = None,
                            col_colors_data = colors_data,
                            colors_legend_data = colors_legend_data, vmin=-1e-08,
                            norm=Normalize(0,4), legend_data=legend_data, 
                            linecolor='#ffffff', linewidth=0.0)
    c.construction()
    c.set_coords()
    c.set_labels()
    c.clustergrid.ax_heatmap.set_frame_on(False)
    c.clustergrid.ax_heatmap.tick_params(pad=1, width=0.05, length=0.3, direction='out')
    c.clustergrid.cax.set_title(label='Association\nScore', pad=12, size=14)
    
    c.clustergrid.savefig('./PN_components/association_matrix.png', 
                          dpi=600, format='png', pad_inches=0.02)
    #c.clustergrid.savefig('./PN_components/association_matrix.tiff', 
                          #dpi=600, format='tiff', pad_inches=0.02)

