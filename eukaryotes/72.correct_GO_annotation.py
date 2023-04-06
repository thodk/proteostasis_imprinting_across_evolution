#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
from constant_variables import mongodb
import pandas
import os
import pymongo
import copy
sys.path.append('../')
import core_functions
import core_classes




if __name__ == '__main__':
    
    '''
    DESCRIPTION:
    
    This script is used to correct the GO annotation based on the True Path Rule.
    New (corrected) MongoDB collections are created and the old ones are
    removed.
    '''
    
    if not os.path.exists('./graph_instances'):
        os.makedirs('./graph_instances')
        
    df = pandas.read_csv('./files/final_species_set_for_analysis.tsv', 
                         sep='\t', index_col=0)
    
    client = pymongo.MongoClient()
    db = client[mongodb]
    ontology = 'GO_P'
    
    for k, species in enumerate(df.abbreviation.tolist()):

        collection = db["_".join([species, "mapping", ontology, "corrected"])]

        # Get the GO annotation
        G = core_functions.construct_graph_from_mongo(ontology, connections="briefly")
        initial_mapping = core_functions.get_mapping_from_mongo(ontology, species, 
                                                                corrected_mapping=False)
        G.set_reference_pool(initial_mapping, mapping_correction=False)
        to_remove_list = core_functions.find_unannotated(G)
        G = core_functions.remove_unannotated(G, to_remove_list)
        G_entries = copy.deepcopy(G.entries)
        
        # Correct the GO annotation based on the True Path Rule        
        G_corrected = core_classes.Graph(ontology, species)
        G_corrected.entries = copy.deepcopy(G_entries)
        G_corrected.set_reference_pool(initial_mapping, mapping_correction=True)
        corrected_mapping = []
        for term, genes_list in G_corrected.reference_mapping.items():
            for gene in genes_list:
                entry = {"term_accession":term, "gene_symbol":gene}
                corrected_mapping.append(entry)

        # Remove the initial collection of GO annotation and create a new one
        # the for corrected annotation
        collection = db["_".join([species, "mapping", ontology, "corrected"])]
        collection.remove({})
        collection.insert(corrected_mapping)
        collection = db["_".join([species, "mapping", ontology])]
        collection.remove({})
        collection.drop()
        
        print(species, ontology)
        print("initial mapping :", len(initial_mapping))
        print("initial amount of terms: ", len(G.entries))
        print("corrected mapping: ", len(corrected_mapping))
        print("corrected amount of terms: ", len(G_corrected.entries))

    client.close()
