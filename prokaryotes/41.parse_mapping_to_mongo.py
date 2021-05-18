# -*- coding: utf-8 -*-
"""
Created on Wed May 30 15:43:14 2018

@author: thodoris
"""

import sys
import re
import pymongo
sys.path.append('../../../BioInfoMiner2/bim2/python_scripts/core/')
import core_classes
import core_functions
sys.path.append('../../../BioInfoMiner2/bim2/python_scripts/core/databases/')
import databases_core
import pandas




def worker(category):
    
    df = pandas.read_csv('./'+prokaryotes_class+'/preprocessing/STEP3_ranked_organisms_with_names.tsv', 
                         sep='\t', index_col=0)

    abbr_names = df.abbreviation.tolist()

    for organism in abbr_names:
        print organism
        BM = databases_core.BiomartMapping(organism, 'GO', mongo_database='BIM_background',
                                           tsv_files_dir='./'+prokaryotes_class+'/data/annotation_data/')
        BM.read_mapping_obo()
        BM.import_to_mongo_obo()



if len(sys.argv) == 1:
    for prokaryotes_class in ['archaea', 'bacteria']:
        worker(prokaryotes_class)        
else:
    prokaryotes_class = sys.argv[1]
    worker(prokaryotes_class)       





client = pymongo.MongoClient()
db = client['BIM_background']

for collection in db.collection_names():
    if re.search('GO_C|GO_F', collection):
        tmp_col = db[collection]
        tmp_col.drop()