#!/usr/bin/python
# -*- coding: utf-8 -*-
from constant_variables import mongodb
import pandas
import pymongo
import re
import sys
sys.path.append('../')
import core_functions




if __name__ == '__main__':
    
    '''
    DESCRIPTION:
    
    This script read the GO annotation of each species in the 'tsv' folder
    and creates a MongoDB collection. Finally only the collections for GO BP
    will be stored, as the analysis will be performed using the its genomic
    annotation. GO CC and GO MF collections are removed.
    '''
    
    df = pandas.read_csv('./files/final_species_set_for_analysis.tsv', 
                         sep='\t', index_col=0)
    
    for species in df.abbreviation.tolist():
        BM = core_functions.BiomartMapping(species, ontology='GO', 
                                           mongo_database=mongodb,
                                           tsv_files_dir='./tsv/')
        BM.read_mapping_obo()
        BM.import_to_mongo_obo()
    
    client = pymongo.MongoClient()
    db = client['background']
    for collection in db.collection_names():
        if re.search('GO_C|GO_F', collection):
            tmp_col = db[collection]
            tmp_col.drop()