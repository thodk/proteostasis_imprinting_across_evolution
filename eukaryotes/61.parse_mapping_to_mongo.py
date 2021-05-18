# -*- coding: utf-8 -*-
"""
Created on Tue May 22 16:17:08 2018

@author: thodoris
"""

import sys
import pandas
import pymongo
sys.path.append('../../../BioInfoMiner2/bim2/python_scripts/core/')
import core_functions
import core_classes
sys.path.append('../../../BioInfoMiner2/bim2/python_scripts/core/databases/')
import databases_core
import os
import re


df = pandas.read_csv('./files/organisms_for_analysis.tsv', sep='\t', index_col=0)
organisms = df.abbreviation.tolist()
print len(organisms)

def biomart_mapping(organisms, ontology):
    for organism in organisms:
        BM = databases_core.BiomartMapping(organism, ontology, mongo_database='BIM_background',
                                           tsv_files_dir = './tsv/')
        BM.read_mapping_obo()
        BM.import_to_mongo_obo()

biomart_mapping(organisms, 'GO')




client = pymongo.MongoClient()
db = client['BIM_background']

for collection in db.collection_names():
    if re.search('GO_C|GO_F', collection):
        tmp_col = db[collection]
        tmp_col.drop()