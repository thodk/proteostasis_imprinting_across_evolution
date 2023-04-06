#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import pandas
import numpy
from scipy.spatial import distance
from scipy.cluster import hierarchy
from sklearn.metrics import homogeneity_score
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
SH.y_align = 0.3



if __name__ == '__main__':
    
    '''
    DESCRIPTION:

    This script creates the 4rth plot of tha manuscript, in order to
    compare the silhouette and homogeneity scores among the 20 conserved mechanisms,
    PN and ribosomal sequences, when they are used as criteria for the separation
    of the 3 taxonomic kingdoms.
    '''

    categories_dict = {
        'eukaryotes':{
                'column':'phylum',
                'key':'Phylum',
                },
        'bacteria':{
                'column':'phylum',
                'key':'Phylum',
                },
        'archaea':{
                'column':'class',
                'key':'Class'
                }
    }
    
    scenarios_dict = {
        'rRNA':{
            'dir':'./rRNAs_analysis/',
            'color':'#56862e',
            'style':'dotted',
            'label': 'Ribosomal Sequences'
        },
        'proteostasis':{
            'dir':'./PN_analysis/',
            'color':'#026690',
            'style':'solid',
            'label': 'Proteostasis Semantic Network'
        },
        'hsp40':{
            'dir':'./HSPs_analysis/hsp40/',
            'color':'#a62c25',
            'style':'dashed',
            'label': 'HSP40 Sequences'
        },
        'hsp70':{
            'dir':'./HSPs_analysis/hsp70/',
            'color': '#aa5325',
            'style':'dashdot',
            'label': 'HSP70 Sequences'
        }
    }
    
    
    subcategories = ['archaea', 'bacteria', 'eukaryotes']
    scenarios = ['rRNA', 'proteostasis', 'hsp70', 'hsp40']
    
    results = {}
    for subcategory in subcategories:
        species_df = pandas.read_csv('./files/species_taxonomies.tsv', sep='\t')
        species_df = species_df.loc[species_df.taxonomy == subcategory, :]
        
        column = categories_dict[subcategory]['column']
        species_taxonomy = dict(zip(species_df.species_name, species_df.loc[:,column]))
    
        ref_labels = []
        for species in species_df.species_name:
            tax = species_taxonomy[species]
            ref_labels.append(tax)
    
        length = len(list(set(ref_labels)))
        N = int(numpy.math.ceil(1.2*length))
        categories_dict.setdefault(subcategory, {}).update({'classes':length})    
    
        for scenario in scenarios:
            main_dir = scenarios_dict[scenario]['dir']
            tmp_data = []
            df = pandas.read_csv(main_dir+'final_distance_matrix.tsv', sep='\t', index_col=0)
            df = df.loc[species_df.species_name, species_df.species_name]
            matrix = df.values
            data = []
            for n in range(N, 1, -1):
                condensed_matrix = distance.squareform(matrix, force='to_vector')
                Z = hierarchy.linkage(condensed_matrix, method='ward')
                labels = list(hierarchy.fcluster(Z, criterion='maxclust',t=n))
                hs = round(homogeneity_score(labels_true=ref_labels, labels_pred=labels),2)
                tmp_data.append([n,hs])
            tmp_df = pandas.DataFrame(tmp_data, columns=['n_clusters', 'aim'])
            results.setdefault(subcategory, {}).update({scenario:tmp_df})                 
    
    


    fig = plt.figure(1, figsize=(6,7), frameon=False)
    grids = GridSpec(11, 1, hspace=0.4, wspace=0)
    step = 3
    
    for i, subcategory in enumerate(subcategories):
        ax = fig.add_subplot(grids[i*step+i:(i+1)*step+i,0])
        for scenario in scenarios:
            tmp = results[subcategory][scenario]
            ax.plot(tmp.n_clusters, tmp.aim, label=scenario, linewidth=1.5, 
                    linestyle='solid', color=scenarios_dict[scenario]['color'])
        ax.axvline(x=categories_dict[subcategory]['classes'], linewidth=1, color='black')
        ax.grid(True, linestyle='--', alpha=0.5, linewidth=0.7)
        ax.set_ylim([0,1.1])
        if tmp.n_clusters.max() > 16:
            ax.set_xticks(range(tmp.n_clusters.min(), tmp.n_clusters.max()+2,2))
        else:
            ax.set_xticks(tmp.n_clusters)
        ax.tick_params(labelsize=8)
        #ax.set_yticklabels(ax.get_yticklabels(), fontdict={'size':8})
        ax.set_title(subcategory.capitalize() + ' ('+categories_dict[subcategory]['key']+'-level)', fontsize=9)
        ax.set_xlabel('Number of Clusters', fontsize=8)
        ax.set_ylabel('Homogeneity Score', fontsize=8)
        
    
    patches = []
    patches_labels = []
    
    p = Line2D([0], [0], marker='o', color=scenarios_dict['rRNA']['color'], 
               markerfacecolor=scenarios_dict['rRNA']['color'],
               markersize=9, label=scenarios_dict['rRNA']['label'])
    patches.append(p)
    patches_labels.append(scenarios_dict['rRNA']['label'])
    
    p = Line2D([0], [0], marker='o', color=scenarios_dict['proteostasis']['color'], 
               markerfacecolor=scenarios_dict['proteostasis']['color'],
               markersize=9, label=scenarios_dict['proteostasis']['label'])
    patches.append(p)
    patches_labels.append(scenarios_dict['proteostasis']['label'])
    
    p = Line2D([0], [0], marker='o', color=scenarios_dict['hsp40']['color'], 
               markerfacecolor=scenarios_dict['hsp40']['color'],
               markersize=9, label=scenarios_dict['hsp40']['label'])
    patches.append(p)
    patches_labels.append(scenarios_dict['hsp40']['label'])
    
    p = Line2D([0], [0], marker='o', color=scenarios_dict['hsp70']['color'], 
               markerfacecolor=scenarios_dict['hsp70']['color'],
               markersize=9, label=scenarios_dict['hsp70']['label'])
    patches.append(p)
    patches_labels.append(scenarios_dict['hsp70']['label'])
    
    
    
    plt.legend(handles=patches, 
               labels=patches_labels, 
               bbox_to_anchor=(0.23, -0.85, 0.6, 0.5),
               fontsize=7, ncol=2, fancybox=False, framealpha=0, 
               handlelength=1, handletextpad=1,
               handleheight=1,
               labelspacing=1,
               handler_map={matplotlib.lines.Line2D: SH})
    
    main_dir = './subtaxonomies_figure/'
    if not os.path.exists(main_dir):
        os.makedirs(main_dir)
    #results_df.to_csv(main_dir+'scores.tsv', sep='\t')
    #plt.savefig(main_dir+'hs_plots.tiff', dpi=800, format='tiff', bbox_inches='tight', pad_inches=0.05)
    plt.savefig(main_dir+'hs_plots.png', dpi=800, format='png', bbox_inches='tight', pad_inches=0.05)
    #plt.savefig(main_dir+'hs_plots.pdf', dpi=600, format='pdf', bbox_inches='tight', pad_inches=0.05)
        