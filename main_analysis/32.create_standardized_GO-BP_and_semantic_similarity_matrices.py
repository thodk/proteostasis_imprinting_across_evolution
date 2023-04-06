#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import pandas
import numpy
import json
import os
sys.path.append('../')
from core_functions import remove_unannotated
from core_functions import construct_graph_from_mongo
from core_functions import get_mapping_from_mongo
from core_semantics import Semantics


if __name__ == '__main__':
    
    main_dir = './PN_analysis/standardized_graph/'
    if not os.path.exists(main_dir):
        os.makedirs(main_dir)
    
    # Load the initial GO-BP graph
    G = construct_graph_from_mongo('GO_P')
    
    # Find GO BP terms annotated for at least one species
    species_df = pandas.read_csv('./files/species.tsv', sep='\t', index_col=0)
    union = []
    for species in species_df.abbreviation:
        mapping = get_mapping_from_mongo('GO_P', species, corrected_mapping=True)
        union = list(set(union).union(mapping.keys()))

    # Terms to remove
    to_remove_terms = list(set(G.entries.keys()).difference(union))
    G = remove_unannotated(G, to_remove_terms)
    semantics = Semantics(G)
    terms_details = []
    for term in G.entries.keys():
        sm = semantics.get_semantic_value(term, 'graph_corpus')
        ic = semantics.get_information_content(term, 'graph_corpus')
        terms_details.append([term, sm, ic, G.get_entry_obj(term).definition])
    semantics_df = pandas.DataFrame(terms_details, columns=['term_id', 
                                    'semantic_value', 'information_content',
                                    'definition'])
    semantics_df.to_csv(main_dir+'semantics.tsv', sep='\t')

    high_ic = round(numpy.percentile(semantics_df.information_content, 20),3)
    low_ic = 0.3
    high_sem_value = round(numpy.percentile(semantics_df.semantic_value,20),3)

    low_sem_value = 0
    substitutions_dict = {}
    for term in G.entries.keys():
        new_terms = semantics.get_ancestors_from_bounded_graph(term, low_ic=low_ic,
                                                               high_ic=high_ic,
                                                               high_sem_value=high_sem_value,
                                                               low_sem_value=low_sem_value)
        substitutions_dict.update({term:new_terms})
    
    final_terms = list(set([i for j in substitutions_dict.values() for i in j]))
    for term, sub_terms in substitutions_dict.items():
        if sub_terms[0] == term:
            continue
        else:
            pass
        if len(set(G.get_entry_obj(term).descendants).intersection(final_terms)) > 0:
            substitutions_dict[term] = [term]
    final_terms = list(set([i for j in substitutions_dict.values() for i in j]))


    with open(main_dir+'GO_P_terms_substitutions.json', 'w') as f:
        json.dump(substitutions_dict, f)

    stand_G = remove_unannotated(G, list(set(G.entries.keys()).difference(final_terms)))
    semantics = Semantics(stand_G)
    
    final_terms_data = []
    for term in final_terms:
        final_terms_data.append([term, stand_G.get_entry_obj(term).definition,
                                 str(semantics.get_information_content(term, 'graph_corpus')),
                                 str(semantics.get_semantic_value(term, 'graph_corpus'))
                                 ])

    tmp_df = pandas.DataFrame(final_terms_data, columns=['term_id', 'definition',
                                                         'ic', 'semantic_value'])
    tmp_df = tmp_df.sort_values(by='ic', ascending=True)
    tmp_df.to_csv(main_dir+'final_terms.csv')
    terms = list(sorted(final_terms))
    
    resnik_mica_matrix = semantics.get_pairwise_similarity(terms, 'resnik', 
                                                           ancestors_set='mica')
    resnik_xgrasm_matrix = semantics.get_pairwise_similarity(terms, 'resnik', 
                                                           ancestors_set='xgrasm')                                                       
    agg_ic_matrix = semantics.get_pairwise_similarity(terms, 'aggregate_ic', 
                                                      ancestors_set='mica')
    mean_matrix = (resnik_mica_matrix + resnik_xgrasm_matrix + agg_ic_matrix)/3.
    
    pandas.DataFrame(resnik_mica_matrix, columns=terms, index=terms).to_csv(main_dir+'resnik_mica_matrix.csv', sep=',')
    pandas.DataFrame(resnik_xgrasm_matrix, columns=terms, index=terms).to_csv(main_dir+'resnik_xgrasm_matrix.csv', sep=',')
    pandas.DataFrame(agg_ic_matrix, columns=terms, index=terms).to_csv(main_dir+'agg_ic_matrix.csv', sep=',')
    pandas.DataFrame(mean_matrix, columns=terms, index=terms).to_csv(main_dir+'mean_sem_sim_matrix.csv', sep=',')

