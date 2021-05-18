# -*- coding: utf-8 -*-

import pandas
import os
import sys
import numpy
import json
import operator
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
from sklearn.manifold import MDS
sys.path.append('../')
from core_functions import remove_unannotated
from core_functions import construct_graph_from_mongo
import core_classes
import heatmap
import matplotlib
matplotlib.rcParams['font.sans-serif'] = ['FreeSans', ]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['axes.titlepad'] = 4
matplotlib.rcParams['xtick.major.pad']='0'
matplotlib.rcParams['ytick.major.pad']='0'


def worker_for_components(clusters_to_keep):
    
    to_keep_terms = []
    for cluster in clusters_to_keep:
        cluster_obj = GCP.entries[cluster]
        members = cluster_obj.get_all_members()
        to_keep_terms.extend(members)

    enriched_terms = [i for j in outputs_mapping.values() for i in j]
    enriched_terms = sorted(list(set(enriched_terms)))

    to_remove_terms = []
    for term in enriched_terms:
        new_terms = substitutions_dict[term]
        if len(set(new_terms).intersection(to_keep_terms)) == 0:
            to_remove_terms.append(term)
        else:
            pass    

    final_terms = list(set(enriched_terms).difference(to_remove_terms))    
    output_df = terms_sem_sim_matrix.loc[final_terms,final_terms]
    return output_df


def worker_for_species(tmp_terms_sim_matrix):

    final_terms = tmp_terms_sim_matrix.index.tolist()
    output_matrix = numpy.ones((len(species), len(species)))
    
    for i in range(len(species)):
        terms1 = list(set(outputs_mapping[species[i]]).intersection(final_terms))
        for j in range(i+1,len(species)):
            terms2 = list(set(outputs_mapping[species[j]]).intersection(final_terms))
            print i, j
            if len(terms1) == 0 and len(terms2) == 0:
                pair_sim = 0
            else:
                pair_matrix_1 = tmp_terms_sim_matrix.loc[terms1, terms2]
                pair_matrix_2 = tmp_terms_sim_matrix.loc[terms2, terms1]                
                pair_sim = semantics.get_average_best_matches([pair_matrix_1, pair_matrix_2])
            output_matrix[i,j] = output_matrix[j,i] = round(pair_sim, 3)
    
    tmp = [species_names[o] for o in species]
    
    zero_line_indices = list(numpy.where(numpy.all(output_matrix == 0, axis=1))[0])
    output_matrix = numpy.delete(arr=output_matrix, axis=0, obj=zero_line_indices)
    output_matrix = numpy.delete(arr=output_matrix, axis=1, obj=zero_line_indices)        
    tmp = numpy.delete(numpy.array(tmp), obj=zero_line_indices)
    
    output_matrix = 1 - output_matrix
    numpy.fill_diagonal(output_matrix, 0)
    output_matrix = pandas.DataFrame(output_matrix, index=tmp, columns=tmp)
    return output_matrix


def worker_for_single_scatter_plot(matrix_df, output):        
    mds_obj = MDS(dissimilarity='precomputed', n_components=2, n_jobs=4, 
                  random_state=100, eps=1e-6, n_init=10)
    colors = [colors_taxonomy_mapping[species_names_tax_mapping[i]] for i in matrix_df.index]
    coords = mds_obj.fit_transform(matrix_df.values)
    return coords, colors


def worker_for_all_scatter_plot(ax, coords, colors, title):
    ax.scatter(coords[:,0], coords[:,1], c=colors, s=5, linewidths=.0, edgecolor='black',  alpha=0.85)
    ax.set_xlim([-lim,lim])
    ax.set_ylim([-lim,lim])
    plt.xticks(ticks, size=0)
    plt.yticks(ticks, size=0)
    plt.title(title, size=8)
    plt.xticks(ticks, size=5)
    plt.xlabel('Eigen-component 1', size=8, labelpad=1)
    plt.yticks(ticks, size=5)
    plt.ylabel('Eigen-component 2', size=8, labelpad=1)
    plt.tick_params(labelsize=6, pad=0.75, width=0.75, length=2, direction='out')
    #plt.setp(ax.get_xticklabels(), visible=False)
    #plt.setp(ax.get_yticklabels(), visible=False)
    ax.grid(True, linestyle='--', alpha=0.5, linewidth=0.3)
    
    





if __name__ == '__main__':

    main_dir = './PN_components/'
    clustering_dir = main_dir+'terms_clustering/' 
    if not os.path.exists(clustering_dir):
        os.makedirs(clustering_dir)

    G = construct_graph_from_mongo('GO_P', mongo_database='background')
    
    with open('./PN_analysis/standardized_graph/GO_P_terms_substitutions.json', 'r') as f:
        substitutions_dict = json.load(f)

    all_terms = list(set([i for i in substitutions_dict.keys()]))
    to_remove_terms = list(set(G.entries.keys()).difference(all_terms))
    G = remove_unannotated(G, to_remove_terms)
    semantics = core_classes.Semantics(G)


    print 'LOADED GRAPH', len(G.entries)

    '''
    STEP 1: Construction of association matrix and heatmap between the main
    PN components and species
    '''

    # load pathway analysis results
    pathway_analysis_outputs = os.listdir('./PN_analysis/pathway_analysis_outputs/')
    outputs_mapping = {}
    for f in pathway_analysis_outputs:
        df = pandas.read_csv('./PN_analysis/pathway_analysis_outputs/'+f, sep=',')
        tmp_species = f.replace('_GO_P.csv', '')
        terms = df.term_id.tolist()
        ics = [ (term, semantics.get_information_content(term, criterion='graph_corpus')) for term in terms]
        ics = sorted(ics, key=operator.itemgetter(1), reverse=True)
        filter_terms = [k[0] for k in filter(lambda x: x[1] >= 0.25, ics)]
        df = df.loc[df.term_id.isin(filter_terms), :]
        length = df.shape[0]
        if length < 100:
            pass
        elif df.loc[df.corrected_pvalue <= 0.05].shape[0] >= 100:
            df = df.loc[df.corrected_pvalue <= 0.05]
        else:
            df = df.iloc[0:100,:]
        all_terms = df.term_id
        outputs_mapping.update({tmp_species: all_terms})

    print 'LOADED PA RESULTS'
    # substitute terms with their ancestors from the standardized GO-BP
    new_outputs_mapping = {}
    for tmp_species, terms in outputs_mapping.items():
        new_terms = []
        for term in terms:
            new_terms.extend(substitutions_dict[term])
        new_terms = list(set(new_terms))
        new_outputs_mapping.update({tmp_species:new_terms})
    final_terms = list(set([i for l in new_outputs_mapping.values() for i in l]))
    print 'CLUSTERING INIT', len(final_terms)
    # Resnik clustering
    GCP = core_classes.GraphClusteringProcessMICA(entries=final_terms, graph_instance=G, 
                                                  metric='resnik', criterion='graph_corpus',
                                                  replace_method='mica')
    GCP.clusteringT(threshold=0.175)
    clusters = sorted(GCP.entries.keys())
    species = sorted(new_outputs_mapping.keys())
    print 'CLUSTERING FINISH', len(clusters)
    clustering_df = pandas.DataFrame({'cluster_id':[], 'term_id':[] })
    for cluster, obj in GCP.entries.items():
        members = obj.get_all_members()
        for m in members:
            clustering_df = clustering_df.append({'cluster_id':cluster,
                                                  'term_id':m},
                                                  ignore_index=True)
    clustering_df.to_csv(clustering_dir+'clustering.csv', sep=',')

    # Construct the association matrix
    M = numpy.zeros((len(clusters),len(species)))
    print 'association matrix', M.shape
    for i, cluster in enumerate(clusters):
        obj = GCP.entries[cluster]
        members = obj.get_all_members()
        for j, tmp_species in enumerate(species):
            if len(set(new_outputs_mapping[tmp_species]).intersection(members)) > 0:
                M[i,j] = 1
            else:
                M[i,j] = 0

    association_df = pandas.DataFrame(M, columns=species, index=clusters)            
    common_components = list(association_df.index[association_df.apply(lambda x: sum(x)/float(len(x)) >= 0.9, axis=1)])
    different_componenets = list(set(clusters).difference(common_components))
    defs = [G.get_entry_obj(i).definition for i in clusters]
    association_df = pandas.DataFrame(M, columns=species, index=defs)
    association_df.to_csv(main_dir+'PN_components_matrix.tsv', sep='\t')
    association_matrix = association_df.values


    # Construct the heatmap of the association matrix
    colors_df = pandas.read_csv('./files/colors.csv', sep=',')
    colors_taxonomy_mapping = dict(zip(colors_df.taxonomy, colors_df.color))
    
    species_df = pandas.read_csv('./files/species.tsv', sep='\t')
    species_names = dict(zip(species_df.abbreviation, species_df.species_name))
    species_abbr_tax_mapping = dict(zip(species_df.abbreviation, species_df.taxonomy))
    species_names_tax_mapping = dict(zip(species_df.species_name, species_df.taxonomy))
    
    row_clustering_data={'data':association_matrix, 'method':'ward', 'metric':'hamming', 'ratio':0.0}
    col_clustering_data={'data':association_matrix.transpose(), 'method':'ward', 'metric':'hamming', 
                         'ratio':0.1}
    y_labels = [i for i in list(association_df.index)]
    y_labels_for_plot = [i.capitalize() for i in y_labels]
    y_labels_dict = dict(zip(y_labels, y_labels_for_plot))
    
    x_labels = [i for i in list(association_df.columns.values)]
    x_labels_for_plot = [species_names[i] for i in x_labels]
    x_labels_dict = dict(zip(x_labels, x_labels_for_plot))
    
    colors = [colors_taxonomy_mapping[species_abbr_tax_mapping[i]] for i in x_labels]
    
    colors_data = {'colors':colors, 'ratio':0.04}
    colors_legend_data = {}
    colors_legend_data.update({'patches': [[], []]})
    for key, value in colors_taxonomy_mapping.items():
        colors_legend_data['patches'][0].append(key.capitalize())
        p = Line2D([0], [0], marker='o', color=value, markerfacecolor=value,
                   markersize=12, label=key.capitalize())
        colors_legend_data['patches'][1].append(p)
    
    colors_legend_data.update({'title':'Taxonomy'})
    colors_legend_data.update({'bbox_anchor':(0.165,0.5), 'fontsize':8, 
                               'handlelength':0.7, 'handletextpad':1,
                               'handleheight':2, 'title_size':10, 'y_align':0.4})
    
    legend_labels = ['NO', 'YES']
    legend_title = 'Association'
    legend_data = {'x':0.1, 'y':0.5, 'w':0.02, 'h':0.15, 'labels':legend_labels, 
                   'labels_size':8, 'cbar_kws':{'ticks':[0.25, 0.75]}, 
                   'title':legend_title, 'title_size':10}
    
    x_specific_labels_format = {}
    for label in x_labels:
        tmp_dict = {'color': colors_taxonomy_mapping[species_abbr_tax_mapping[label]],
                    'weight':600}
        x_specific_labels_format.update({x_labels_dict[label]:tmp_dict})
    
    x_axis_data = {'labels': x_labels_for_plot,
                   'specific_labels_format':x_specific_labels_format, 
                   'fontdict':{'size':1, 'rotation':90}}
    
    y_specific_labels_format = {}
    for i, label in enumerate(y_labels):
        if clusters[i] in common_components:
            tmp_dict = {'color': '#B30000', 'weight':600}
            y_specific_labels_format.update({y_labels_dict[label]:tmp_dict})
    
    y_axis_data = {'labels': y_labels_for_plot,
                   'specific_labels_format':y_specific_labels_format,
                   'fontdict':{'size':5}}
    
    heatmap_data={'data':association_matrix, 'type':'features', 'x_ratio':1.5, 'y_ratio':1.2}
    cmap = LinearSegmentedColormap.from_list("my_colormap", ('#eaeaea', '#000000'), N=2)
    
    c = heatmap.Clustergram(heatmap_data, figsize=(7,5.5), cmap=cmap, 
                            y_axis_data=y_axis_data,
                            x_axis_data=x_axis_data,
                            row_clustering_data=row_clustering_data,
                            col_clustering_data=col_clustering_data,
                            row_colors_data = None,
                            col_colors_data = colors_data,
                            colors_legend_data = colors_legend_data,
                            vmin=0.0, vmax=1, legend_data=legend_data, 
                            linecolor='#e0e0e0', linewidth=0.0)
    c.construction()
    c.set_coords()
    c.set_labels()
    c.clustergrid.ax_heatmap.set_frame_on(False)
    c.clustergrid.ax_heatmap.tick_params(pad=1, width=0.05, length=0.3, direction='out')

    output = main_dir+'figure3_A'
    c.clustergrid.savefig(output+'.png', dpi=800, format='png')
    #c.clustergrid.savefig(main_dir+output+'.tiff', dpi=800, format='tiff')
    c.clustergrid.savefig(output+'.pdf', dpi=800, format='pdf')
    





    '''
    STEP 2: Constrution of scatter plots for the different partitions of the 
    PN components, based on the MDS
    '''

    '''
    scaffold_similarity_matrix = pandas.read_csv('./PN_analysis/standardized_graph/mean_sem_sim_matrix.csv', 
                                                 sep=',', index_col=0)
    enriched_terms = [i for j in outputs_mapping.values() for i in j]
    enriched_terms = sorted(list(set(enriched_terms)))
    terms_sem_sim_matrix = numpy.ones((len(enriched_terms), len(enriched_terms)))
    for t1 in range(len(enriched_terms)):
        term1 = enriched_terms[t1]
        new_term1 = substitutions_dict[term1]
        for t2 in range(t1+1, len(enriched_terms)):
            term2 = enriched_terms[t2]
            new_term2 = substitutions_dict[term2]
            sim = scaffold_similarity_matrix.loc[new_term1, new_term2].values.mean()
            terms_sem_sim_matrix[t1,t2] = terms_sem_sim_matrix[t2,t1] = round(sim,3)
    
    terms_sem_sim_matrix = pandas.DataFrame(terms_sem_sim_matrix, 
                                            columns=enriched_terms, 
                                            index=enriched_terms)
    '''
    terms_sem_sim_matrix = pandas.read_csv('./PN_components/enriched_terms_semantic_similarities.csv',
                                           sep=',', index_col=0)


print 'SCATTER PLOTS'
all_componenets = clusters
all_components_sim_matrix = worker_for_components(all_componenets)
species_sim_matrix = worker_for_species(all_components_sim_matrix)
species_sim_matrix.to_csv(main_dir+'all_components_sim_matrix.csv', sep=',')
#species_sim_matrix = pandas.read_csv(main_dir+'all_components_sim_matrix.csv', sep=',', index_col=0)
coords1, colors1 = worker_for_single_scatter_plot(species_sim_matrix, main_dir+'all_componenents')
df1 = pandas.DataFrame(coords1, index=species)
df1.loc[:,'colors'] = colors1
df1.to_csv(main_dir+'all_componenets_coords.csv', sep=',')

common_components_sim_matrix = worker_for_components(common_components)
species_sim_matrix2 = worker_for_species(common_components_sim_matrix)
species_sim_matrix2.to_csv(main_dir+'common_components_sim_matrix.csv', sep=',')
#species_sim_matrix2 = pandas.read_csv(main_dir+'common_components_sim_matrix.csv', sep=',', index_col=0)
coords2, colors2 = worker_for_single_scatter_plot(species_sim_matrix2, main_dir+'common')
df2 = pandas.DataFrame(coords2, index=species)
df2.loc[:,'colors'] = colors2
df2.to_csv(main_dir+'common_components_coords.csv', sep=',')


different_components_sim_matrix = worker_for_components(different_componenets)
species_sim_matrix3 = worker_for_species(different_components_sim_matrix)
species_sim_matrix3.to_csv(main_dir+'differential_components_sim_matrix.csv', sep=',')
#species_sim_matrix3 = pandas.read_csv(main_dir+'differential_components_sim_matrix.csv', sep=',', index_col=0)
coords3, colors3 = worker_for_single_scatter_plot(species_sim_matrix3, main_dir+'common')
df3 = pandas.DataFrame(coords3, index=species)
df3.loc[:,'colors'] = colors3
df3.to_csv(main_dir+'differential_components_coords.csv', sep=',')




plt.close()
lim = 0.7
ticks = numpy.arange(-0.6, 0.6+0.1, 0.2)

fig = plt.figure(1, figsize=(7,2), frameon=False)
grids = GridSpec(1,3, hspace=0, wspace=0.3)


ax = fig.add_subplot(grids[0,1])
worker_for_all_scatter_plot(ax, coords2, colors2, 'Common Components of\nProteostasis Network')

ax = fig.add_subplot(grids[0,2])
worker_for_all_scatter_plot(ax, coords3, colors3, 'Differential Components of\nProteostasis Network')

ax = fig.add_subplot(grids[0,0])
worker_for_all_scatter_plot(ax, coords1, colors1, 'Entire Proteostasis\nNetwork')


output = main_dir+'figure3_B'
plt.savefig(output+'.png', dpi=800, format='png', bbox_inches='tight', pad_inches=0.05)
#plt.savefig(output+'.tiff', dpi=800, format='tiff', bbox_inches='tight', pad_inches=0.05)
plt.savefig(output+'.pdf', dpi=600, format='pdf', bbox_inches='tight', pad_inches=0.05)

