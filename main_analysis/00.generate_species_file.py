#!/usr/bin/python
# -*- coding: utf-8 -*-
import pandas
import os




if __name__ == '__main__':
    
    '''
    DESCRIPTION:

    This script gets the 'final_species_set_for_analysis.tsv' files
    from eukaryotes and prokaryotes directories and create a new one,
    which record all the necessary information about the species that will
    be used for the analysis.
    
    Output: 'species.tsv' file in 'files' directory.
    '''
    
    # Retrieve the data for eukaryotes
    filename = '../eukaryotes/files/final_species_set_for_analysis.tsv'
    eukaryotes_df = pandas.read_csv(filename, sep='\t', index_col=0)
    eukaryotes_df = eukaryotes_df.loc[:, ['species_name', 'species_abbreviation',
                                          'abbreviation', 'ensembl_class', 'rRNA']]
    eukaryotes_df.loc[:, 'taxonomy'] = ['eukaryotes']*len(eukaryotes_df)
    
    # Retrieve the data for prokaryotes
    filename = '../prokaryotes/files/final_species_set_for_analysis.tsv'    
    prokaryotes_df = pandas.read_csv(filename, sep='\t', index_col=0)
    prokaryotes_df = prokaryotes_df.loc[:, ['species_name', 'abbreviation', 
                                            'prokaryotes_class', 'rRNA']]    
    prokaryotes_df.rename(columns={'prokaryotes_class':'taxonomy'}, inplace=True)
    
    # Create the concatenated data frame and store it
    df = pandas.concat([eukaryotes_df, prokaryotes_df], axis=0, ignore_index=True,
                       sort=True)
    df = df.loc[:, ['species_name','species_abbreviation', 'abbreviation',
                    'taxonomy', 'ensembl_class', 'rRNA']]
    if not os.path.exists('./files'):
        os.makedirs('./files')
    df.to_csv('./files/species.tsv', sep='\t')


