#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import pandas
import os
sys.path.append('../../../BioInfoMiner/core/')
import bioinfominer_interface



if __name__ == '__main__':
    
    mongodb = 'background'
    main_dir = './other_terms_analysis/'
    terms = os.listdir(main_dir+'pathway_analysis_inputs/')
    
    for term in terms[15:]:

        output_dir = main_dir+'pathway_analysis_outputs/'+term+'/'
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        files = os.listdir(main_dir+'pathway_analysis_inputs/'+term+'/')
        ontology = 'GO_P'
    
        for en, f in enumerate(files):
            print(term, en, len(files))
            label = os.path.basename(f).split('.')[0]
            tmp_species = '_'.join(label.split('_')[0:2])
            #if os.path.isfile(output_dir + label + '_GO_P.csv'):
            #    continue
           
            F = open(main_dir+'pathway_analysis_inputs/'+term+'/'+f, 'r')
            input_list = [i.replace('\n', '').replace('\r', '') for i in F.readlines()]
            
            input_list = [{'gene_symbol':i} for i in input_list]
            ea_results = bioinfominer_interface.steptwo(input_list, ontology=ontology, 
                                                        organism=tmp_species,
                                                        corrected_mapping=True,
                                                        hyper_pvalue_threshold=0.1,
                                                        corrected_pvalue_threshold=1,
                                                        mongo_database=mongodb)
            
            cols = ['rank', 'term_id', 'term_definition', 'enrichment', 
                    'hyper_pvalue', 'corrected_pvalue', 'genes']    
            df = pandas.DataFrame(ea_results['GO_P'], columns=cols)
            df.to_csv(output_dir + label + '_GO_P.tsv', index=None, sep='\t')


