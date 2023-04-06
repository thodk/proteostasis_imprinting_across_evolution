#!/usr/bin/python
# -*- coding: utf-8 -*-
import pandas
import os
import re
from constant_variables import prokaryotic_classes
from constant_variables import define_main_dir


def worker(prokaryotic_class):

    main_dir = define_main_dir(prokaryotic_class)+'preprocessing/'
    files = [main_dir+'scanning/'+i for i in os.listdir(main_dir+'scanning/')]
    files = filter(lambda x: re.search('details', x), files)
    details_list = []
    for f in files:
        df = pandas.read_csv(f, sep='\t')
        details_list.append(df)

    details_df = pandas.concat(details_list, ignore_index=True)
    pos_bool1 = details_df.hsp70_bool == True
    pos_bool2 = details_df.hsp40_bool == True
    details_df = details_df.loc[((pos_bool1) & (pos_bool2)), :]
    details_df = details_df.sort_values(by=['terms_amount', 'genes_amount'], ascending=False)   
    details_df.to_csv(main_dir+'01_ranked_species_based_on_genomic_annotation.tsv', sep='\t')




if __name__ == "__main__":
   
    
    '''
    DESCRIPTION:
    
    Given the execution of step 0, here, proteomes without annotated HSP40 and 
    HSP70 are filtered out and the final list is ranked according to the amount 
    of annotated terms and the number of genes.
    
    Output: 01_ranked_species_based_on_genomic_annotation.tsv
    '''
    
    # Execute the 'worker' function for each prokaryotic class
    [worker(c) for c in prokaryotic_classes]
