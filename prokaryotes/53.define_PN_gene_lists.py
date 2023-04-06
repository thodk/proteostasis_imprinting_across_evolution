#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas
import pymongo
import os
import sys
sys.path.append('../')
from constant_variables import define_main_dir
from core_functions import mongodb


def worker_for_prokaryotes(prokaryotic_class):
    
    # Get the species from this specific prokaryotic 'category'
    main_dir = main_dir = define_main_dir(prokaryotic_class)+'preprocessing/'
    results_dir = './data/'+prokaryotic_class+'/data/PN_gene_lists/'
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    
    species_df = pandas.read_csv(main_dir+'03_species_with_rRNA_and_HSPs_annotation.tsv',
                                 sep='\t', index_col=0)
    species_dict = dict(zip(species_df.species_name, species_df.abbreviation))
    
    # Connect to MongoDB to find the genes that have been characterized as
    # proteostasis-related prokaryotic genes and their are included in the
    # annotation of these species.
    client = pymongo.MongoClient()
    db = client[mongodb]
    for species, abbreviation in species_dict.items():
        collection = db[abbreviation + '_mapping_GO_P_corrected']
        results = collection.find({'gene_symbol':{'$in':prokaryotic_genes}})
        tmp_genes = []
        for i in results:
            tmp_genes.append(i['gene_symbol'])
        tmp_genes = list(set(tmp_genes))
        F = open(results_dir+abbreviation+'.txt', 'w')
        F.write('\n'.join(tmp_genes))
        F.close()
    

    

if __name__ == '__main__':

    '''
    DESCRIPTION:
        
    This script creates the input files for pathway analysis.
    '''
    
    prokaryotic_genes = pandas.read_csv('./files/prokaryotic_genes_for_proteostasis.tsv',
                                        sep='\t', index_col=0).gene_symbol.tolist()

    worker_for_prokaryotes('bacteria')
    worker_for_prokaryotes('archaea')

