#!/usr/bin/python
# -*- coding: utf-8 -*-
import pandas
from constant_variables import define_main_dir
from constant_variables import ensembl_classes


def worker(ensembl_class):
     
    main_dir = define_main_dir(ensembl_class)
    files_dir = main_dir+'files/'

    # '02_organisms_with_rRNA_annotation.tsv' file should have been created 
    # manually!
    df = pandas.read_csv(main_dir+'02_species_with_rRNA_annotation.tsv',
                         sep='\t', index_col=0)
    datasets_df = pandas.read_csv(files_dir+'/ensembl_datasets.tsv', sep='\t')
    datasets_df = datasets_df.loc[:, ['abbreviation', 'dataset']]
    merged_df = pandas.merge(df, datasets_df, on='abbreviation', how='inner')
    merged_df = merged_df.rename({'abbreviation': 'species_abbreviation',
                                  'species': 'species_name'}, axis=1)
    
    # Create a new column, 'abbreviation', which will work as an id for each 
    # species:
    # examples: ensembl_0, ensembl_1, ensembl_2, plant_10, plant_11 etc. 
    merged_df.loc[:,'abbreviation'] = [ensembl_class+'_'+str(i) for i in range(len(merged_df))]
    merged_df.to_csv(main_dir+'03_species_with_rRNA_and_HSPs_annotation.tsv', sep='\t')




if __name__ == "__main__":
    
    '''
    DESCRIPTION:

    Before the execution of this step, it is mandatory to create manually the 
    '02_organisms_with_rRNA_annotation.tsv' file. This file contains information
    about the species with at least one annotated ribosomal RNA sequence.

    The rRNA sequences should have been retrieved manually from the ENA 
    repository. The whole process is quite time consuming. 

    Finally the 'worker' function examines which species have both HSPs and rRNA
    annotation and write that information in the 
    '03_species_with_rRNA_and_HSPs_annotation.tsv' file. These species will be 
    further examined in the following steps, based on their genomic GO annotation 
    to define the set that will be used for the analysis.
    '''

    # Execute the 'worker' function for each Ensembl class
    [worker(c) for c in ensembl_classes]