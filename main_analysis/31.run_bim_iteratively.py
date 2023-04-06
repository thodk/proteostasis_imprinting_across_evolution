# -*- coding: utf-8 -*-
import sys
import pandas
import os
sys.path.append('../../../BioInfoMiner/core/')
import bioinfominer_interface



if __name__ == '__main__':
    
    mongodb = 'background'
    
    output_dir = './PN_analysis/pathway_analysis_outputs/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    files = os.listdir('./PN_analysis/pathway_analysis_inputs/')
    ontology = 'GO_P'
    
    for f in files:
        
        label = os.path.basename(f).split('.')[0]
        #if os.path.isfile(output_dir + label + '_GO_P.csv'):
        #    continue
       
        F = open('./PN_analysis/pathway_analysis_inputs/'+f, 'r')
        input_list = [i.replace('\n', '').replace('\r', '') for i in F.readlines()]
        print(label, len(input_list))
        
        input_list = [{'gene_symbol':i} for i in input_list]
        ea_results = bioinfominer_interface.steptwo(input_list, ontology=ontology, 
                                                    organism=label,
                                                    corrected_mapping=True,
                                                    hyper_pvalue_threshold=0.1,
                                                    corrected_pvalue_threshold=1,
                                                    mongo_database=mongodb)
        
        cols = ['rank', 'term_id', 'term_definition', 'enrichment', 
                'hyper_pvalue', 'corrected_pvalue', 'genes']    
        df = pandas.DataFrame(ea_results['GO_P'], columns=cols)
        df.to_csv(output_dir + label + '_GO_P.tsv', index=None, sep='\t')


