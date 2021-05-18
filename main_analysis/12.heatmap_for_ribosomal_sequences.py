# -*- coding: utf-8 -*-

import sys
import pandas
import numpy
sys.path.append('../')
import heatmap
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
import matplotlib
matplotlib.rcParams['font.sans-serif'] = ['FreeSans', ]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['axes.titlepad'] = 4


if __name__ == '__main__':
    
    df = pandas.read_csv('./rRNAs_analysis/final_distance_matrix.tsv', sep='\t', 
                         index_col=0)
    matrix = df.values
    
    row_clustering_data={'data':matrix, 'method':'ward', 'metric':None, 'ratio':0.15}
    col_clustering_data={'data':matrix, 'method':'ward', 'metric':None, 'ratio':0.15}
    labels = [i for i in list(df.index)]
    labels_for_plot = [i.capitalize() for i in labels]
    labels_dict = dict(zip(labels, labels_for_plot))
    
    colors_df = pandas.read_csv('./files/colors.csv', sep=',')
    colors_taxonomy_mapping = dict(zip(colors_df.taxonomy.tolist(), colors_df.color.tolist()))
    
    df = pandas.read_csv('./files/species.tsv', sep='\t')
    species = df['abbreviation'].tolist()
    species_taxonomy_mapping = dict(zip(df.species_name.tolist(), df.taxonomy.tolist()))
    colors = [colors_taxonomy_mapping[species_taxonomy_mapping[i]] for i in labels]
    
    colors_data = {'colors':colors, 'ratio':0.02}
    colors_legend_data = {}
    colors_legend_data.update({'patches': [[], []]})
    for key, value in colors_taxonomy_mapping.items():
        colors_legend_data['patches'][0].append(key.capitalize())
        p = Line2D([0], [0], marker='o', color=value, markerfacecolor=value,
                   markersize=14, label=key.capitalize())
        colors_legend_data['patches'][1].append(p)
    
    
    colors_legend_data.update({'title':'Taxonomy\n'})
    colors_legend_data.update({'bbox_anchor':(0.165,0.45), 'fontsize':10, 
                               'handlelength':1.2, 'handletextpad':0.5,
                               'handleheight':2, 'title_size':12})
    
    legend_labels = numpy.arange(0,1.01,0.2).round(2)
    legend_title = 'Normalized\nDistance\n'
    legend_data = {'x':0.09, 'y':0.47, 'w':0.03, 'h':0.25, 'labels':legend_labels, 
                   'labels_size':10, 'cbar_kws':{'ticks':legend_labels}, 
                   'title':legend_title, 'title_size':12}
                   
    
    specific_labels_format = {}
    for label in labels:
        tmp_dict = {'color': colors_taxonomy_mapping[species_taxonomy_mapping[label]],
                    'weight':900}
        specific_labels_format.update({labels_dict[label]:tmp_dict})
    
    x_axis_data = {'labels': labels_for_plot, 'specific_labels_format':specific_labels_format, 
                   'fontdict':{'size':0.4, 'rotation':90}}
    y_axis_data = {'labels': labels_for_plot, 'specific_labels_format':specific_labels_format,
                   'fontdict':{'size':0.4}}
                   
    heatmap_data={'data':matrix, 'type':'distances', 'x_ratio':0.8, 'y_ratio':0.8}
    cmap = LinearSegmentedColormap.from_list("my_colormap", ('#eaeaea', '#000000'), N=10)
    
    c = heatmap.Clustergram(heatmap_data, figsize=(8,8), cmap=cmap, 
                            y_axis_data=y_axis_data,
                            x_axis_data=x_axis_data,
                            row_clustering_data=row_clustering_data,
                            col_clustering_data=col_clustering_data,
                            row_colors_data = colors_data,
                            col_colors_data = colors_data,
                            colors_legend_data = colors_legend_data,
                            vmin=0.0, vmax=1, legend_data=legend_data, 
                            linecolor='#e0e0e0', linewidth=0.005)
    c.construction()
    c.set_coords()
    c.set_labels()
    
    c.clustergrid.ax_heatmap.xaxis.set_tick_params(width=0.3, length=1.5, pad=1)
    c.clustergrid.ax_heatmap.yaxis.set_tick_params(width=0.3, length=1.5, pad=1)
    #c.clustergrid.savefig('./rRNAs_analysis/rRNAs_heatmap.tiff', dpi=800, format='tiff',
    #                      pad_inches=0.02)
    c.clustergrid.savefig('./rRNAs_analysis/rRNAs_heatmap.png', dpi=800, format='png',
                          pad_inches=0.02)