#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import os
import pandas
import numpy
import operator
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

    main_dir = './PN_analysis/'
    
    G = construct_graph_from_mongo('GO_P')
    with open(main_dir+'standardized_graph/GO_P_terms_substitutions.json', 'r') as f:
        substitutions_dict = json.load(f)
    all_terms = list(set([i for i in substitutions_dict.keys()]))
    to_remove_terms = list(set(G.entries.keys()).difference(all_terms))
    G = remove_unannotated(G, to_remove_terms)
    semantics = core_classes.Semantics(G)
    
    bim_outputs = os.listdir(main_dir+'pathway_analysis_outputs/')
    mapping = {}
    for f in bim_outputs:

        df = pandas.read_csv(main_dir+'pathway_analysis_outputs/'+f, sep='\t')
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
    
    enriched_terms = [i for j in mapping.values() for i in j]
    enriched_terms = sorted(list(set(enriched_terms)))

    sem_sim_matrix = pandas.read_csv(main_dir+'standardized_graph/mean_sem_sim_matrix.csv',
                                     sep=',', index_col=0)
    enriched_terms_sim_matrix = numpy.ones((len(enriched_terms), len(enriched_terms)))
    for t1 in range(len(enriched_terms)):
        print(t1, len(enriched_terms))
        term1 = enriched_terms[t1]
        new_term1 = substitutions_dict[term1]
        for t2 in range(t1+1, len(enriched_terms)):
            term2 = enriched_terms[t2]
            new_term2 = substitutions_dict[term2]
            sim = sem_sim_matrix.loc[new_term1, new_term2].values.mean()
            enriched_terms_sim_matrix[t1,t2] = enriched_terms_sim_matrix[t2,t1] = round(sim,3)
            
    enriched_terms_sim_matrix_df = pandas.DataFrame(enriched_terms_sim_matrix, 
                                                    index=enriched_terms, 
                                                    columns=enriched_terms)
    enriched_terms_sim_matrix_df.to_csv(main_dir+'enriched_terms_sem_similarities.tsv', 
                                        sep='\t')


    species = sorted(mapping.keys())
    matrix = numpy.zeros((len(species), len(species)))
    for i in range(len(species)):
        print(i, len(species))
        terms1 = mapping[species[i]]
        for j in range(i+1,len(species)):
            terms2 = mapping[species[j]]
            pair_matrix_1 = enriched_terms_sim_matrix_df.loc[terms1, terms2]
            pair_matrix_2 = enriched_terms_sim_matrix_df.loc[terms2, terms1]
            pair_sim = semantics.get_average_best_matches([pair_matrix_1, pair_matrix_2])
            matrix[i,j] = matrix[j,i] = round(pair_sim, 3)
    species_df = pandas.read_csv('./files/species.tsv', sep='\t', index_col=0)
    names = dict(zip(species_df.abbreviation, species_df.species_name))
    c = [names[i] for i in species]
    f = 1 - matrix
    numpy.fill_diagonal(f, 0)
    pandas.DataFrame(f, index=c, columns=c).to_csv(main_dir+'final_distance_matrix.tsv', sep='\t')
    

    
