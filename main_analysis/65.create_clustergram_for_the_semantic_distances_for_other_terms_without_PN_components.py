#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import pandas
import numpy
import os
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import matplotlib
matplotlib.rcParams['font.sans-serif'] = ['FreeSans', ]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['axes.titlepad'] = 4
sys.path.append('../')
import heatmap
from core_functions import taxonomy_colors


if __name__ == '__main__':
    
    '''
    DESCRIPTION:

    This script is used to illustrate a semantic distance matrix as a heatmap 
    and calculate the clustergram of species, based on the Ward's minimum variance 
    method. The illustrated distances matrices belong to the semantic profiles of 
    the 20 conserved mechanisms, without the inclusion of PN components. 
    The 'heatmap' module is used for this process (called from the upper 
    directory). This is a custom module which is based on 'seaborn.clustergram'
    and configures appropriately the plot for these specific data.
    '''
    
    main_dir = './other_terms_analysis/'
    other_terms = os.listdir(main_dir+'pathway_analysis_inputs/')
    output_dir = main_dir + 'semantic_analysis_without_PN_components/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # run the process for each term in the main_dir
    for other_term in other_terms:

        df = pandas.read_csv(output_dir+other_term+'/'+'final_distance_matrix.tsv', 
                             sep='\t', index_col=0)
        matrix = df.values
        
        # Define some paramters for the 'heatmap' module
        row_clustering_data={'data':matrix, 'method':'ward', 'metric':None, 'ratio':0.10}
        col_clustering_data={'data':matrix, 'method':'ward', 'metric':None, 'ratio':0.10}
        labels = [i for i in list(df.index)]
        labels_for_plot = [i.capitalize() for i in labels]
        labels_dict = dict(zip(labels, labels_for_plot))
        
        # Color parameters
        df = pandas.read_csv('./files/species.tsv', sep='\t')
        species_taxonomy_mapping = dict(zip(df.species_name.tolist(), df.taxonomy.tolist()))
        colors = [taxonomy_colors[species_taxonomy_mapping[i]] for i in labels]    
        colors_data = {'colors':colors, 'ratio':0.02}
        colors_legend_data = {}
        colors_legend_data.update({'patches': [[], []]})
        for key, value in taxonomy_colors.items():
            colors_legend_data['patches'][0].append(key.capitalize())
            p = Line2D([0], [0], marker='o', color=value, markerfacecolor=value,
                       markersize=14, label=key.capitalize())
            colors_legend_data['patches'][1].append(p)
        
        
        colors_legend_data.update({'title':'Taxonomy\n'})
        colors_legend_data.update({'bbox_anchor':(0.205,0.275), 'fontsize':12, 
                                   'handlelength':1, 'handletextpad':0.25,
                                   'handleheight':1.5, 'title_size':14})
        
        legend_labels = numpy.arange(0,1.01,0.2).round(2)
        legend_title = 'Normalized\nDistance\n'
        legend_data = {'x':0.1, 'y':0.35, 'w':0.03, 'h':0.25, 'labels':legend_labels, 
                       'labels_size':12, 'cbar_kws':{'ticks':legend_labels}, 
                       'title':legend_title, 'title_size':14}
        
        # Specify x and y labels format
        specific_labels_format = {}
        for label in labels:
            tmp_dict = {'color': taxonomy_colors[species_taxonomy_mapping[label]],
                        'weight':900}
            specific_labels_format.update({labels_dict[label]:tmp_dict})
        
        x_axis_data = {'labels': labels_for_plot, 'specific_labels_format':specific_labels_format, 
                       'fontdict':{'size':0.4, 'rotation':90}}
        y_axis_data = {'labels': labels_for_plot, 'specific_labels_format':specific_labels_format,
                       'fontdict':{'size':0.4}}
                       
        # Define some necessary paramters for the 'heatmap' module
        heatmap_data={'data':matrix, 'type':'distances', 'x_ratio':0.8, 'y_ratio':0.8}  
        cmap = LinearSegmentedColormap.from_list("my_colormap", ('#eaeaea', '#000000'), N=100)
        
        # Call the function to initiate the clustergram construction
        c = heatmap.Clustergram(heatmap_data, figsize=(8,8), cmap=cmap, 
                                y_axis_data=y_axis_data,
                                x_axis_data=x_axis_data,
                                row_clustering_data=row_clustering_data,
                                col_clustering_data=col_clustering_data,
                                row_colors_data = colors_data,
                                col_colors_data = colors_data,
                                colors_legend_data = colors_legend_data,
                                norm=Normalize(0,1), legend_data=legend_data, 
                                linecolor='#e0e0e0', linewidth=0.00)
        
        # Plot construction
        c.construction()
        c.set_coords()
        c.set_labels()
        c.clustergrid.ax_heatmap.xaxis.set_tick_params(width=0.3, length=1.5, pad=1)
        c.clustergrid.ax_heatmap.yaxis.set_tick_params(width=0.3, length=1.5, pad=1)
        c.clustergrid.cax.set_title(label='Normalized\nDistance', pad=12, size=14)
        # Save it as png or tiff
        #c.clustergrid.savefig(main_dir+hsp+'_heatmap.tiff', dpi=800, format='tiff',
        #                  pad_inches=0.02)
        c.clustergrid.savefig(output_dir+other_term+'/'+'clustergram.png', dpi=900, format='png',
                              pad_inches=0.02)
