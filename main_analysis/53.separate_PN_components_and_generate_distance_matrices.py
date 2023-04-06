#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas
import os
import sys
import numpy
import json
import operator
sys.path.append('../')
sys.path.append('/home/thodoris/Projects/proteostasis/proteostasis_imprinting_across_evolution/')
import core_classes
from core_functions import remove_unannotated
from core_functions import construct_graph_from_mongo


def worker(list_of_terms, key):
    
    tmp_sem_sim_matrix = sem_sim_matrix.loc[list_of_terms, list_of_terms]
    species = sorted(terms_df.species.unique())
    output_matrix = numpy.ones((len(species), len(species)))
    for i in range(len(species)):
        terms1 = species_terms[species[i]][key]
        for j in range(i+1,len(species)):
            terms2 = species_terms[species[j]][key]
            if len(terms1) == 0 and len(terms2) == 0:
                pair_sim = 0
            else:
                pair_matrix_1 = tmp_sem_sim_matrix.loc[terms1, terms2]
                pair_matrix_2 = tmp_sem_sim_matrix.loc[terms2, terms1]                
                pair_sim = semantics.get_average_best_matches([pair_matrix_1, pair_matrix_2])
            output_matrix[i,j] = output_matrix[j,i] = round(pair_sim, 3)

    species_names = [species_names_mapping[s] for s in species]
    output_matrix = 1 - output_matrix
    numpy.fill_diagonal(output_matrix, 0)
    output_matrix = pandas.DataFrame(output_matrix, index=species_names,
                                     columns=species_names)
    return output_matrix


if __name__ == '__main__':

    #os.chdir('/home/thodoris/Projects/proteostasis/proteostasis_imprinting_across_evolution/main_analysis/')
    # Load the GO BP graph from MongoDB
    G = construct_graph_from_mongo('GO_P', mongo_database='background')
    
    # Load the GO_P terms substitution dictionary
    filename = './PN_analysis/standardized_graph/GO_P_terms_substitutions.json'
    with open(filename, 'r') as f:
        substitutions_dict = json.load(f)
    
    # Find all the annotated terms (which are included as keys in the above 
    # substitutions dictionary) and create a new graph, keeping only these 
    # terms
    all_terms = list(set([i for i in substitutions_dict.keys()]))
    to_remove_terms = list(set(G.entries.keys()).difference(all_terms))
    G = remove_unannotated(G, to_remove_terms)
    semantics = core_classes.Semantics(G)


    # STEP 1
    # Load the results of pathway analysis and store them in a pandas data frame
    pathway_analysis_outputs = os.listdir('./PN_analysis/pathway_analysis_outputs/')
    tmp_dfs = []
    for f in pathway_analysis_outputs:
        tmp_df = pandas.read_csv('./PN_analysis/pathway_analysis_outputs/'+f, sep='\t')
        tmp_species = f.replace('_GO_P.tsv', '')
        terms = tmp_df.term_id.tolist()
        ics = [(term, semantics.get_information_content(term, criterion='graph_corpus')) for term in terms]
        ics = sorted(ics, key=operator.itemgetter(1), reverse=True)
        # 1. Terms with overall IC < 0.1 are filtered out.
        filter_terms = [k[0] for k in filter(lambda x: x[1] >= 0.1, ics)]
        tmp_df = tmp_df.loc[tmp_df.term_id.isin(filter_terms), :]
        # 2. An extra variable is created ('score') which is equal to 
        # score = -log10(corrected_pvalue) + 1
        tmp_df.loc[:, 'score'] = -numpy.log10(tmp_df.corrected_pvalue) + 1
        # 3. Filter terms dataframe, based on term set size and statistics
        if tmp_df.shape[0] < 100:
            pass
        elif tmp_df.loc[tmp_df.corrected_pvalue <= 0.05].shape[0] >= 100:
            tmp_df = tmp_df.loc[tmp_df.corrected_pvalue <= 0.05]
        else:
            tmp_df = tmp_df.iloc[0:100,:]
        tmp_df.loc[:, 'species'] = [tmp_species]*len(tmp_df)
        tmp_df = tmp_df.loc[:, ['species', 'term_id']]
        tmp_dfs.append(tmp_df)
    terms_df = pandas.concat(tmp_dfs, axis=0)
    terms_df.reset_index(drop=True, inplace=True)
    
    
    # STEP 2
    # Load the NMF components and split terms in the conserved and differential 
    # parts. Terms with non-zero coefficients in all components belong to the
    # conserved part and all the other in the differential one. 
    nmf_df = pandas.read_csv('./PN_components/nmf/components_matrix.tsv', 
                             sep='\t', index_col=0)
    nmf_common = nmf_df[(nmf_df > 0).all(axis=1)]
    nmf_differential = nmf_df[~(nmf_df > 0).all(axis=1)]
    common_terms = []
    differential_terms = []
    for term, obj in G.entries.items():
        if obj.definition in nmf_common.index:
            common_terms.append(term)
        elif obj.definition in nmf_differential.index:
            differential_terms.append(term)


    # STEP 3
    # Create two dictionaries, for common and differential terms for each species.
    # If an enriched term (from pathway analysis results) is or has an ancestor
    # in the list of common_terms, then it belongs to that category, otherwise
    # it will be included in the differential part.
    species_terms = {}
    for tmp_species, tmp_df in terms_df.groupby(by='species'):
        tmp_common = []
        tmp_differential = []
        for term in tmp_df.term_id.tolist():
            ancestors = G.get_entry_obj(term).ancestors + [term]
            if len(set(ancestors).intersection(common_terms)) > 0:
                tmp_common.append(term)
            else:
                tmp_differential.append(term)
        tmp_common = list(set(tmp_common))
        tmp_differential = list(set(tmp_differential))
        species_terms.update({tmp_species: {'common' : tmp_common, 
                                            'differential': tmp_differential}})


    # STEP 4
    # Create a distance matrix for the species, based on the sematnic similarities
    # of common and differential components
    
    df = pandas.read_csv('./files/species.tsv', sep='\t', index_col=0)
    species_names_mapping = dict(zip(df.abbreviation, df.species_name))
    
    sem_sim_matrix = pandas.read_csv('./PN_analysis/enriched_terms_sem_similarities.tsv',
                                     sep='\t', index_col=0)
    
    all_common_terms = [t for s,v in species_terms.items() for t in v['common']]
    all_common_terms = list(set(all_common_terms))
    all_differential_terms = [t for s,v in species_terms.items() for t in v['differential']]
    all_differential_terms = list(set(all_differential_terms))
    
    len(all_common_terms), len(all_differential_terms)

    if not os.path.exists('./PN_components/common_and_differential/'):
        os.makedirs('./PN_components/common_and_differential/')
    
    dis_matrix = worker(all_common_terms, 'common')
    filename = './PN_components/common_and_differential/common_components_distance_matrix.tsv'
    dis_matrix.to_csv(filename, sep='\t')
    filename = './PN_components/common_and_differential/differential_components_distance_matrix.tsv'
    dis_matrix = worker(all_differential_terms, 'differential')
    dis_matrix.to_csv(filename, sep='\t')
    


    