#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import json
import os
import pandas
import numpy
import pymongo
import sys
sys.path.append('../')
from core_functions import mongodb
from core_functions import remove_unannotated
from core_functions import construct_graph_from_mongo
import core_classes

if __name__ == '__main__':
    

    terms = ['GO:0006457', 'GO:0051641', 'GO:0008033', 'GO:0006974', 
             'GO:0009260', 'GO:0008104', 'GO:0051252', 'GO:0032259', 
             'GO:0006629', 'GO:0006355', 'GO:0022607', 'GO:0033554', 
             'GO:0015031', 'GO:0006096', 'GO:0060255', 'GO:0016310', 
             'GO:0006605', 'GO:0006310', 'GO:0006281', 'GO:0046034']
    
    
    # Load the GO BP graph from MongoDB
    G = construct_graph_from_mongo('GO_P', mongo_database='background')
    
    # Load the GO_P terms substitution dictionary
    filename = './PN_analysis/standardized_graph/GO_P_terms_substitutions.json'
    with open(filename, 'r') as f:
        substitutions_dict = json.load(f)
    
    
    len(G.entries)
    
    # Find all the annotated terms (which are included as keys in the above 
    # substitutions dictionary) and create a new graph, keeping only these 
    # terms
    all_terms = list(set([i for i in substitutions_dict.keys()]))
    to_remove_terms = list(set(G.entries.keys()).difference(all_terms))
    G = remove_unannotated(G, to_remove_terms)
    semantics = core_classes.Semantics(G)
    
    len(G.entries)
    
    selected_terms_df = []
    for term in terms:
        ic = semantics.get_information_content(term, criterion='graph_corpus')
        selected_terms_df.append([term, G.get_entry_obj(term).definition, ic])
    selected_terms_df = pandas.DataFrame(selected_terms_df,
                                         columns=['term_id', 'definition', 'ic'])


    # STEP 1
    # Load PN analysis gene lists and their sizes
    PN_data = {}
    for f in os.listdir('./PN_analysis/pathway_analysis_inputs/'):
        F = open('./PN_analysis/pathway_analysis_inputs/'+f, 'r')
        size = len(F.readlines())
        species = f.replace('.txt', '')
        PN_data.update({species: size})

    # STEP 2
    # Load other terms annotation data
    client = pymongo.MongoClient()
    db = client[mongodb]
    numpy.random.seed(1234)
    other_species_data = {}
    for tmp_species in PN_data.keys():
        print(tmp_species)
        collection = db[tmp_species+'_mapping_GO_P_corrected']
        results = collection.aggregate([
            {'$match': {'term_accession': {'$in': selected_terms_df.term_id.tolist()}}},
            {'$group': {'_id': "$term_accession", 
                        'list': {'$addToSet': '$gene_symbol'}}}
            ])
        size = PN_data[tmp_species]
        for entry in results:
            gene_set = entry['list']
            if len(gene_set) <= size:
                pass
            else:
                gene_set = list(numpy.random.choice(gene_set, size=size, replace=False))
            other_species_data.setdefault(entry['_id'], {}).update({tmp_species:gene_set})
    client.close()
    
    other_species_data.keys()
    # STEP 3:
    # Write the gene sets as inputs for pathway analysis
    main_dir = './other_terms_analysis/'
    for i, row in selected_terms_df.iterrrows():
        definition = row.definition.replace(' ', '_')
        tmp_dir = main_dir+'pathway_analysis_inputs/'+definition+'/'
        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)



        #F = open(tmp_dir+tmp_species+'.txt', 'w')
        #F.write('\n'.join(gene_set))
        #F.close()
        