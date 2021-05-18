# -*- coding: utf-8 -*-

import pandas
from constant_variables import define_main_dir
from constant_variables import parameters_df

def worker(ensembl_class):
    
    main_dir = define_main_dir(ensembl_class)
    f1 = main_dir+'files/species_with_lack_of_annotation.tsv'
    df1 = pandas.read_csv(f1, sep='\t', index_col=0)    
    df1 = df1.loc[df1.loc[:, 'pass'] == False, :]
    to_remove = df1.species.unique()

    f2 = main_dir+'03_species_with_rRNA_and_HSPs_annotation.tsv'
    df2 = pandas.read_csv(f2, sep='\t', index_col=0) 
    all_species = df2.species_abbreviation.tolist()
    
    final_species = list(set(all_species).difference(to_remove))
    final_df = df2.loc[df2.species_abbreviation.isin(final_species), :]
    final_df.to_csv(main_dir+'04_species_with_genomic_annotation.tsv', sep='\t')

    return final_df


'''
DESCRIPTION:

1. Reading of lack_of_annotation files
2. Definition of final set of species that will be used for the analysis
'''

if __name__ == "__main__":

    global_df = pandas.DataFrame(data={'abbreviation':[], 'organism':[], 
                                       'eukaryotes_class':[]})
    
    for ensembl_class in parameters_df.loc[:,'class'].tolist():    
        df = worker(ensembl_class)    
        df.loc[:,'ensembl_class'] = [ensembl_class]*len(df)
        global_df = pandas.concat([global_df, df], axis=0, ignore_index=True, sort=True)
        global_df = global_df.reset_index(drop=True)
        
    global_df.to_csv('./files/final_species_set_for_analysis.tsv', sep='\t')
