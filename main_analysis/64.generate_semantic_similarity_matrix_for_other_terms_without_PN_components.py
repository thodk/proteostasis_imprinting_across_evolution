#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import os
import pandas
import numpy
import operator
import collections
import json
sys.path.append('../')
from core_functions import remove_unannotated
from core_functions import construct_graph_from_mongo
import core_classes



if __name__ == '__main__':

    '''
    DESCRIPTION: Calculate the semantic similarity distances based on the 
    standardized GO graph.
    '''    

    main_dir = './other_terms_analysis/'
    other_terms = os.listdir(main_dir+'pathway_analysis_inputs/')
    output_dir = main_dir + 'semantic_analysis_without_PN_components/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    G = construct_graph_from_mongo('GO_P')
    with open('./PN_analysis/standardized_graph/GO_P_terms_substitutions.json', 'r') as f:
        substitutions_dict = json.load(f)    
    all_terms = list(set([i for i in substitutions_dict.keys()]))
    to_remove_terms = list(set(G.entries.keys()).difference(all_terms))
    G = remove_unannotated(G, to_remove_terms)
    semantics = core_classes.Semantics(G)
    
    
    
    for other_term in other_terms[14:]:
        
        mapping = {}
        bim_outputs = os.listdir(main_dir+'pathway_analysis_outputs/'+other_term+'/')
        for f in bim_outputs:
            df = pandas.read_csv(main_dir+'pathway_analysis_outputs/'+other_term+'/'+f, sep='\t')
            species = f.replace('_GO_P.tsv', '')
            # filter ics
            terms = df.term_id.tolist()
    
            ics = [ (term, semantics.get_information_content(term, criterion='graph_corpus')) for term in terms]
            ics = sorted(ics, key=operator.itemgetter(1), reverse=True)
            filter_terms = [k[0] for k in filter(lambda x: x[1] >= 0.1, ics)]
            df = df.loc[df.term_id.isin(filter_terms), :]
            if df.shape[0] < 100:
                pass
            elif df.loc[df.corrected_pvalue <= 0.05].shape[0] >= 100:
                df = df.loc[df.corrected_pvalue <= 0.05]
            else:
                df = df.iloc[0:100,:]
            all_terms = df.term_id
            mapping.update({species: all_terms})
    
        names_df = []
        for key in mapping.keys():
            names_df.append([key, '_'.join(key.split('_')[0:2])])
        names_df = pandas.DataFrame(names_df, columns=['sub','species'])
    
        new_mapping = {}
        for species, tmp_df in names_df.groupby('species'):
            N = tmp_df.shape[0]
            terms_pool = []
            for sub in tmp_df.loc[:,'sub'].tolist():
                sub_terms = mapping[sub]
                terms_pool.extend(sub_terms)
            if N > 1:
                terms_pool = list(filter(lambda x: x[1] > numpy.math.ceil(0.2*N), collections.Counter(terms_pool).items()))
                terms_pool = [t[0] for t in terms_pool]
            else:
                pass
            new_mapping.update({species:terms_pool})
    
    
        PN_df = pandas.read_csv('./PN_components/PN_terms.tsv', sep='\t',
                                index_col=0)
        PN_terms = PN_df.term_id.tolist()
        for tmp_species, tmp_terms in new_mapping.items():
            new_tmp_terms = []
            for tmp_term in tmp_terms:        
                new_term = substitutions_dict[tmp_term]
                new_term = list(set(new_term).difference(PN_terms))
                if len(new_term) == 0:
                    pass
                else:
                    new_tmp_terms.append(tmp_term)
            new_mapping.update({tmp_species:new_tmp_terms})
    
        enriched_terms = [t for l in new_mapping.values() for t in l]
        enriched_terms = sorted(list(set(enriched_terms)))
        enriched_terms = pandas.DataFrame(enriched_terms, columns=['term_id'])
    
        index_mapping = {}
        for tmp_species, tmp_terms in new_mapping.items():
            tmp_indices = enriched_terms.loc[enriched_terms.term_id.isin(tmp_terms), :].index.tolist()
            if len(tmp_indices) > 0:
                index_mapping.update({tmp_species: tmp_indices})
            else:
                index_mapping.update({tmp_species: []})
                
        sem_sim_matrix = pandas.read_csv('PN_analysis/standardized_graph/mean_sem_sim_matrix.csv',
                                         sep=',', index_col=0)
        enriched_terms_sim_matrix = numpy.zeros((len(enriched_terms), len(enriched_terms)))
        for i1, row1 in enriched_terms.iterrows():
            new_term1 = list(set(substitutions_dict[row1.term_id]).difference(PN_terms))
            for i2, row2 in enriched_terms.iloc[i1+1:,].iterrows():
                new_term2 = list(set(substitutions_dict[row2.term_id]).difference(PN_terms))
                sim = sem_sim_matrix.loc[new_term1, new_term2].values.mean()
                enriched_terms_sim_matrix[i1,i2] = enriched_terms_sim_matrix[i2,i1] = round(sim,3)
            print(i1, enriched_terms.shape[0])
        #enriched_terms_sim_matrix_df = pandas.DataFrame(enriched_terms_sim_matrix, 
        #                                                index=enriched_terms, 
        #                                                columns=enriched_terms)
        #enriched_terms_sim_matrix_df.to_csv(main_dir+'enriched_terms_sem_similarities.tsv', 
        #                                    sep='\t')
    
        # Similarities for species
        species = sorted(new_mapping.keys())
        species_sim_matrix = numpy.zeros((len(species), len(species)))
        for i in range(len(species)):
            terms1_indices = index_mapping[species[i]]
            for j in range(i+1,len(species)):
                print(other_term, i, species[i], species[j], len(species))
                terms2_indices = index_mapping[species[j]]
                pair_matrix_1 = enriched_terms_sim_matrix[terms1_indices, ]
                pair_matrix_1 = pair_matrix_1[:, terms2_indices]
                pair_matrix_2 = enriched_terms_sim_matrix[terms2_indices, ]
                pair_matrix_2 = pair_matrix_2[:, terms1_indices]
                if len(pair_matrix_1) == 0 or len(pair_matrix_2) == 0:
                    pair_sim = 0 
                else:
                    pair_sim = semantics.get_average_best_matches([pair_matrix_1, pair_matrix_2])
                species_sim_matrix[i,j] = species_sim_matrix[j,i] = round(pair_sim, 3)
    
        # distance matrix
        species_dist_matrix = 1 - species_sim_matrix
        numpy.fill_diagonal(species_dist_matrix, 0)
    
        species_df = pandas.read_csv('./files/species.tsv', sep='\t', index_col=0)
        names_dict = dict(zip(species_df.abbreviation, species_df.species_name))
        columns = [names_dict[i] for i in species]
    
        if not os.path.exists(output_dir+other_term+'/'):
            os.makedirs(output_dir+other_term+'/')
        final_df = pandas.DataFrame(species_dist_matrix, index=columns, 
                                    columns=columns)
        final_df.to_csv(output_dir+other_term+'/final_distance_matrix.tsv', sep='\t')
        
        
    

    
