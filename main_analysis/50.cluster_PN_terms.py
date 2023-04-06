#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas
import os
import sys
import numpy
import json
import operator
import collections
sys.path.append('../')
sys.path.append('/home/thodoris/Projects/proteostasis/proteostasis_imprinting_across_evolution/')
import core_classes
from core_functions import remove_unannotated
from core_functions import construct_graph_from_mongo

os.chdir('/home/thodoris/Projects/proteostasis/proteostasis_imprinting_across_evolution/main_analysis')
if __name__ == '__main__':

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
        # 3. An extra variable is created to show the significance of a term, based on 
        # pvalues distribution and the amount of enriched terms. The 'significant' terms 
        # have been used in the previous step for the construction of PN-clustergram. 
        # All the other terms are kept as 'non-significant'.
        if tmp_df.shape[0] < 100:
            tmp_df.loc[:, 'significance'] = [1]*len(tmp_df)
        elif tmp_df.loc[tmp_df.corrected_pvalue <= 0.05].shape[0] >= 100:
            #tmp_df = tmp_df.loc[tmp_df.corrected_pvalue <= 0.05]
            tmp_df.loc[:, 'significance'] = tmp_df.corrected_pvalue.apply(lambda x: 1 if x <= 0.05 else 0)
        else:
            tmp_df.loc[:, 'significance'] = tmp_df.index.map(lambda x: 1 if x < 100 else 0)
        tmp_df.loc[:, 'species'] = [tmp_species]*len(tmp_df)
        tmp_df = tmp_df.loc[:, ['species', 'term_id', 'score', 'significance']]
        tmp_dfs.append(tmp_df)
    terms_df = pandas.concat(tmp_dfs, axis=0)
    terms_df.reset_index(drop=True, inplace=True)


    # STEP 2
    # Preliminary terms correction based on the substitutions dictionary
    # 1. Initially, get the corrected terms and create a new column
    tmp_dfs = []
    for tmp_term, tmp_term_df in terms_df.groupby(by='term_id'):
        corrected_terms = substitutions_dict[tmp_term]
        for corrected_term in corrected_terms:
            tmp_term_df.loc[:,'corrected_term_id'] = [corrected_term]*len(tmp_term_df)
            tmp_dfs.append(tmp_term_df.copy(deep=True))
    draft_corr_terms_df = pandas.concat(tmp_dfs, axis=0)
    draft_corr_terms_df.reset_index(drop=True, inplace=True)
    # 2. Then, for each species, get terms' scores and keep the max ones
    tmp_dfs = []
    for tmp_species, tmp_species_df in draft_corr_terms_df.groupby(by='species'):
        scores_df = tmp_species_df.groupby(by='corrected_term_id').score.mean().to_frame()
        significances_df = tmp_species_df.groupby(by='corrected_term_id').significance.max().to_frame()
        changed_terms_df = tmp_species_df.groupby(by='corrected_term_id').term_id.apply(list).to_frame()
        changed_terms_df.columns = ['terms_from_pathway_analysis']
        tmp_species_df = pandas.concat([scores_df, significances_df, changed_terms_df], ignore_index=False, axis=1)
        tmp_species_df.reset_index(drop=False, inplace=True)
        tmp_species_df.loc[:,'species'] = [tmp_species]*len(tmp_species_df)
        tmp_dfs.append(tmp_species_df)
    corr_terms_df = pandas.concat(tmp_dfs, axis=0)
    corr_terms_df.reset_index(drop=True, inplace=True)


    
    # STEP 3
    # Split the corr_terms_df (which contains the results of pathway analysis 
    # with the corrected terms) regarding the 3 taxonomic kingdoms
    kingdoms = ['archaea', 'bacteria', 'eukaryotes']
    kingdoms_terms_dict = {}
    for kingdom in kingdoms:
        # 1. Get the dataframe for the specific kingdom
        if kingdom in ['archaea', 'bacteria']:
            kingdom_df = corr_terms_df.loc[corr_terms_df.species.apply(lambda s: s.startswith(kingdom))]
        else:
            kingdom_df = corr_terms_df.loc[corr_terms_df.species.apply(lambda s: not (s.startswith('bac') or s.startswith('arc')))]
        # 2. Find terms with significance = 1 and calculate the size of kingdom
        kingdom_enriched_terms = kingdom_df.loc[kingdom_df.significance == 1, 'corrected_term_id'].unique().tolist()
        #kingdom_enriched_terms = kingdom_df.loc[kingdom_df.significance.isin([1,0]), 'corrected_term_id'].unique().tolist()

        kingdom_size = len(kingdom_df.species.unique().tolist())
        # 3. Examine the terms of kingdom and separate them based on the following
        # criterion: Terms significantly enriched in more or less that the 90% 
        # of the species
        above_90per_terms = []
        below_90per_terms = []
        for tmp_term in kingdom_enriched_terms:
            bool1 = kingdom_df.corrected_term_id == tmp_term
            bool2 = kingdom_df.significance == 1
            #bool2 = kingdom_df.significance.isin([1,0])
            amount = len(kingdom_df.loc[(bool1) & (bool2), 'species'].unique())
            if amount/kingdom_size >= 0.9:
                above_90per_terms.append(tmp_term)
            else:
                below_90per_terms.append(tmp_term)
        kingdoms_terms_dict.update({kingdom: {'above_90per_terms':above_90per_terms,
                                              'below_90per_terms':below_90per_terms,
                                              'kingdom_df':kingdom_df,
                                              'kingdom_size':kingdom_size}})


    # STEP 4
    # Filtering of the list of significantly enriched terms per kingdom. Each  
    # term, whose descendants are also included in the list, is substituted by 
    # its descendants. The final goal is to create a list with uniquely 
    # semantically defined terms
    above90_correction = []
    for kingdom, kingdom_dict in kingdoms_terms_dict.items():
        kingdom_df = kingdom_dict['kingdom_df']
        above_90per_terms = kingdom_dict['above_90per_terms']
        above_90per_kingdom_df = kingdom_df.loc[kingdom_df.corrected_term_id.isin(above_90per_terms), :]
        tmp_dfs = []    
        for term_id in above_90per_terms:
            term_df = above_90per_kingdom_df.loc[above_90per_kingdom_df.corrected_term_id == term_id,:]
            descendants = list(set(G.get_entry_obj(term_id).descendants).intersection(above_90per_terms))
            if len(descendants) > 0:
                for desc in descendants:
                    term_df.loc[:,'corrected_term_id'] = [desc]*len(term_df)
                    term_df.loc[:,'term_from_st_graph'] = [term_id]*len(term_df)
                    tmp_dfs.append(term_df.copy(deep=True))
            else:
                term_df.loc[:,'corrected_term_id'] = [term_id]*len(term_df)
                term_df.loc[:,'term_from_st_graph'] = [term_id]*len(term_df)
                tmp_dfs.append(term_df.copy(deep=True))      
        above_90per_kingdom_df = pandas.concat(tmp_dfs, axis=0)
        above_90per_kingdom_df.reset_index(drop=True, inplace=True)
        above90_correction.append(above_90per_kingdom_df)
    above90_correction_df = pandas.concat(above90_correction, axis=0)
    above90_correction_df.reset_index(drop=True, inplace=True)
    above90_correction_terms = list(above90_correction_df.corrected_term_id.unique())


    # STEP 5
    # Filtering of the non-significantly enriched terms per kingdom.
    
    # STEP 5A
    # Find the generic GO BP terms with distance 2 from the root and no
    # descendants in the list (uniquely semantically defined)
    generic_terms = []
    for term in G.get_entry_obj('GO:0008150').children:
        generic_terms.extend(G.get_entry_obj(term).children)
    generic_terms = list(set(generic_terms))
    to_remove = []
    for term in generic_terms:
        if len(set(generic_terms).intersection(G.get_entry_obj(term).descendants)) > 0:
            to_remove.append(term)
    generic_terms = list(set(generic_terms).difference(to_remove))

    # STEP 5B
    # Examine the below_90per terms for each species and keep only those with
    # significance = 1. Find their generic ancestors (defined above). In the next
    # step they will be substituted by these generic terms
    below90_correction_terms = []
    for kingdom, kingdom_dict in kingdoms_terms_dict.items():
        kingdom_terms = []
        kingdom_df = kingdom_dict['kingdom_df']
        kingdom_size = kingdom_dict['kingdom_size']
        below_90per_terms = kingdom_dict['below_90per_terms']
        below_90per_kingdom_df = kingdom_df.loc[kingdom_df.corrected_term_id.isin(below_90per_terms), :]
        # 1. Get the highly enriched below90 terms
        below_90per_kingdom_df = below_90per_kingdom_df.loc[below_90per_kingdom_df.significance == 1, :]
        #below_90per_kingdom_df = below_90per_kingdom_df.loc[below_90per_kingdom_df.significance.isin([0,1]), :]
        # 2. Examine each specis separately
        for tmp_species, tmp_species_df in below_90per_kingdom_df.groupby(by='species'):
            tmp_species_terms = []
            for tmp_term in tmp_species_df.corrected_term_id.unique():
                # 3. Substitute the term with its generic ancestors
                ancestors = set(generic_terms).intersection(G.get_entry_obj(tmp_term).ancestors)
                if len(ancestors) > 0:
                    for ancestor in ancestors:
                        tmp_species_terms.append(ancestor)
            tmp_species_terms = list(set(tmp_species_terms))
            kingdom_terms = kingdom_terms + tmp_species_terms
        kingdom_terms = list(collections.Counter(kingdom_terms).items())
        # 3. Keep only the generic terms with percentage grater than 50% in the kingdom
        kingdom_terms = list(filter(lambda x: x[1]>0.5*kingdom_size, kingdom_terms))
        kingdom_terms = [x[0] for x in kingdom_terms]
        for k in kingdom_terms:
            print(kingdom, G.get_entry_obj(k).definition)
        below90_correction_terms.extend(kingdom_terms)
    
    # The final clusters
    clusters = list(set(below90_correction_terms + above90_correction_terms))
    print('number of clusters:', len(clusters))
    #for term in clusters:
        #print(G.get_entry_obj(term).definition)

    # STEP 5C
    # Substitute the below90 terms with their ancestors in the clusters list.
    # If a below90 term has not any ancestor, then it is filtered out, as it is
    # neither significantly enriched nor associated with any other term that is
    # significant for the examined kingdom
    below90_correction_df = []
    for kingdom, kingdom_dict in kingdoms_terms_dict.items():
        kingdom_df = kingdom_dict['kingdom_df']
        below_90per_terms = kingdom_dict['below_90per_terms']
        below_90per_kingdom_df = kingdom_df.loc[kingdom_df.corrected_term_id.isin(below_90per_terms), :]
        tmp_dfs = []
        for term_id in below_90per_terms:
            term_df = below_90per_kingdom_df.loc[below_90per_kingdom_df.corrected_term_id == term_id,:]
            ancestors = set(clusters).intersection(G.get_entry_obj(term_id).ancestors)
            if len(ancestors) > 0:
                for ancestor in ancestors:
                    term_df.loc[:,'corrected_term_id'] = [ancestor]*len(term_df)
                    term_df.loc[:,'term_from_st_graph'] = [term_id]*len(term_df)
                    tmp_dfs.append(term_df.copy(deep=True))
            else:
                pass
        below_90per_kingdom_df = pandas.concat(tmp_dfs, axis=0)
        below_90per_kingdom_df.reset_index(drop=True, inplace=True)
        below90_correction_df.append(below_90per_kingdom_df)
    
    below90_correction_df = pandas.concat(below90_correction_df, axis=0)
    below90_correction_df.reset_index(drop=True, inplace=True)


    # STEP 6
    # Create the final data frame for the heatmap
    
    # 1. Concatenate the transformed data frames for above90 and below90 terms
    tmp_df = pandas.concat([above90_correction_df, below90_correction_df], axis=0)
    tmp_df.reset_index(drop=True, inplace=True)
    
    # 2. Correct the score and significance values - keep the maximum
    tmp_dfs = []
    for tmp_species, tmp_species_df in tmp_df.groupby(by='species'):
        scores_df = tmp_species_df.groupby(by='corrected_term_id').score.mean().to_frame()
        significances_df = tmp_species_df.groupby(by='corrected_term_id').significance.max().to_frame()
        changed_terms_df = tmp_species_df.groupby(by='corrected_term_id').term_from_st_graph.apply(list).to_frame()
        changed_terms_df.columns = ['terms_from_clustering']
        tmp_species_df = pandas.concat([scores_df, significances_df, changed_terms_df], ignore_index=False, axis=1)
        tmp_species_df.reset_index(drop=False, inplace=True)
        tmp_species_df.loc[:,'species'] = [tmp_species]*len(tmp_species_df)
        tmp_dfs.append(tmp_species_df)
    final_df = pandas.concat(tmp_dfs, axis=0)
    final_df.reset_index(drop=True, inplace=True)
    
    if not os.path.exists('./PN_components/'):
        os.makedirs('./PN_components/')
    tmp_df.to_csv('./PN_components/corrected_pathway_analysis_results.tsv',
                    sep='\t')
    final_df.to_csv('./PN_components/final_pathway_analysis_results.tsv',
                    sep='\t')
    

    # STEP 7
    # The final association matrix
    terms = sorted(final_df.corrected_term_id.unique().tolist())
    species = sorted(final_df.species.unique().tolist())
    M = numpy.zeros((len(terms),len(species)))
    '''
    for i, tmp_species in enumerate(species):
        tmp_species_df = final_df.loc[final_df.species == tmp_species,:]
        indexes = []
        scores = []
        for j, row in tmp_species_df.iterrows():
            try:
                indexes.append(terms.index(row.corrected_term_id))
                scores.append(row.score)
            except ValueError:
                pass
        M[indexes, i] = scores
    '''
    for i, tmp_species in enumerate(species):
        for j, tmp_term in enumerate(terms):
            # Search for the score in final_df
            tmp_species_df1 = final_df.loc[final_df.species == tmp_species,:]
            score = tmp_species_df1.loc[tmp_species_df1.corrected_term_id == tmp_term, 'score']
            if len(score) == 1:
                M[j,i] = score
                continue
            else:
                pass
            # Search for the score in corr_terms_df
            tmp_species_df2 = corr_terms_df.loc[corr_terms_df.species == tmp_species,:]
            score = tmp_species_df2.loc[tmp_species_df2.corrected_term_id == tmp_term, 'score']                
            if len(score) == 1:
                M[j,i] = score
                continue
            else:
                pass
            descendants = G.get_entry_obj(tmp_term).descendants[:]
            tmp_desc = list(set(tmp_species_df1.corrected_term_id).intersection(descendants))
            if len(tmp_desc) > 0:
                score = tmp_species_df1.loc[tmp_species_df1.corrected_term_id.isin(tmp_desc), 'score'].mean()
                M[j,i] = score
                continue
            else:
                pass
            tmp_desc = list(set(tmp_species_df2.corrected_term_id).intersection(descendants))
            if len(tmp_desc) > 0:
                score = tmp_species_df2.loc[tmp_species_df2.corrected_term_id.isin(tmp_desc), 'score'].mean()
                M[j,i] = score
                continue
            else:
                pass            

    defs = [G.get_entry_obj(i).definition for i in terms]
    association_df = pandas.DataFrame(M, columns=species, index=defs)
    association_df.to_csv('./PN_components/association_matrix.tsv', sep='\t')



    terms_df = terms_df.loc[terms_df.significance > 0]
    PN_terms = []
    for cluster in clusters:
        descendants = set(G.get_entry_obj(cluster).descendants + [cluster])
        tmp_terms = list(descendants.intersection(terms_df.term_id.tolist()))
        tmp_terms_2 = []
        for t in tmp_terms:
            tmp_terms_2.extend([t])
            tmp_terms_2.extend(substitutions_dict[t])
        PN_terms.extend(tmp_terms_2)
    PN_terms = list(set(PN_terms))
    
    F = open('./PN_components/PN_terms.tsv', 'w')
    F.write('\tterm_id\tdefinition\n')
    for i, term in enumerate(PN_terms):
        F.write(str(i) + '\t' + term + '\t' + G.get_entry_obj(term).definition + '\n')
    F.close()