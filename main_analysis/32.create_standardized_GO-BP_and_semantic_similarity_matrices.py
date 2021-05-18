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
import core_classes


if __name__ == '__main__':
    
    main_dir = './PN_analysis/standardized_graph/'
    if not os.path.exists(main_dir):
        os.makedirs(main_dir)
    
    '''
    STEP 1
    Create the standardized version of GO-BP, by filtering out unannotated
    terms and those with extremely high IC and SV values    
    '''
    G = construct_graph_from_mongo('GO_P') # initial GO-BP graph
    
    # find terms annotated for at least one species
    species_df = pandas.read_csv('./files/species.tsv', sep='\t', index_col=0)
    union = []
    for species in species_df.abbreviation:
        mapping = get_mapping_from_mongo('GO_P', species, corrected_mapping=True)
        union = list(set(union).union(mapping.keys()))

    # terms to remove
    to_remove_terms = list(set(G.entries.keys()).difference(union))
    G = remove_unannotated(G, to_remove_terms)
    semantics = core_classes.Semantics(G)
    terms_details = []
    for term in G.entries.keys():
        sm = semantics.get_semantic_value(term, 'graph_corpus')
        ic = semantics.get_information_content(term, 'graph_corpus')
        terms_details.append([term, sm, ic])
    semantics_df = pandas.DataFrame(terms_details, columns=['term_id', 
                                    'semantic_value', 'information_content'])

    high_ic = round(numpy.percentile(semantics_df.information_content, 20),3)
    low_ic = 0.25
    high_sem_value = round(numpy.percentile(semantics_df.semantic_value,20),3)
    low_sem_value = 0

    substitutions_dict = {}
    final_terms = []
    for term in G.entries.keys():
        new_terms = semantics.get_ancestors_from_bounded_graph(term, low_ic=0.25,
                                                               high_ic=high_ic,
                                                               high_sem_value=high_sem_value,
                                                               low_sem_value=0)
        substitutions_dict.update({term:new_terms})
        final_terms.extend(new_terms)
        
    final_terms = list(set(final_terms))
    with open(main_dir+'GO_P_terms_substitutions.json', 'w') as f:
        json.dump(substitutions_dict, f)



    '''
    STEP 2
    Construct semantic similarity matrices for the terms of standardized GO-BP  
    '''
    
    all_terms = list(set([i for i in substitutions_dict.keys()]))
    final_terms = list(set([i for j in substitutions_dict.values() for i in j]))
    final_terms_data = []
    for term in final_terms:
        final_terms_data.append([term, G.get_entry_obj(term).definition,
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
    
    
    





