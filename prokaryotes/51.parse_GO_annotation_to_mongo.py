#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import re
import pymongo
sys.path.append('../')
import core_functions
import pandas




if __name__ == '__main__':
    
    '''
    DESCRIPTION:
    
    This script read the GO annotation of each species in the 'annotation_data' folder
    and creates a MongoDB collection. Finally only the collections for GO BP
    will be stored, as the analysis will be performed using the its genomic
    annotation. GO CC and GO MF collections are removed.
    '''
    
    for prokaryotes_class in ['archaea', 'bacteria']:
    
        filename = './data/'+prokaryotes_class+'/preprocessing/03_species_with_rRNA_and_HSPs_annotation.tsv'
        annotation_folder = './tsv/'
        df = pandas.read_csv(filename, sep='\t', index_col=0)
        for species in df.abbreviation.tolist():
            BM = core_functions.BiomartMapping(species, ontology='GO', 
                                               mongo_database=core_functions.mongodb,
                                               tsv_files_dir=annotation_folder)
            BM.read_mapping_obo()
            BM.import_to_mongo_obo()

    client = pymongo.MongoClient()
    db = client[core_functions.mongodb]
    
    for collection in db.collection_names():
        if re.search('GO_C|GO_F', collection):
            tmp_col = db[collection]
            tmp_col.drop()