# -*- coding: utf-8 -*-
"""
Created on Wed May 30 18:05:09 2018

@author: thodoris
"""

import sys
import pandas
import os
sys.path.append('../../../BioInfoMiner2/bim2/python_scripts/core/')
import core_functions
import core_classes
import bioinfominer_interface




output_dir = './bim_outputs/'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

files = os.listdir('./bim_inputs/')
ontology = 'GO_P'


for f in files:
    
    label = os.path.basename(f).split('.')[0]

    if os.path.isfile(output_dir + label + '_GO_P.csv'):
        continue
   
    F = open('./bim_inputs/'+f, 'r')
    input_list = [i.replace('\n', '').replace('\r', '') for i in F.readlines()]

    print  label, len(input_list)
    
    input_list = [{'gene_symbol':i} for i in input_list]

    ea_results = bioinfominer_interface.steptwo(input_list, ontology=ontology, 
                                                organism=label,
                                                corrected_mapping=True,
                                                hyper_pvalue_threshold=0.1,
                                                corrected_pvalue_threshold=1)
    
    cols = ['rank', 'term_id', 'term_definition', 'enrichment', 
            'hyper_pvalue', 'corrected_pvalue', 'genes']
    


    df = pandas.DataFrame(ea_results['GO_P'], columns=cols)
    df.to_csv(output_dir + label + '_GO_P.csv', index=None)


