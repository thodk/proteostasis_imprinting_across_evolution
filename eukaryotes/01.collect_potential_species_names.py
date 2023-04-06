#!/usr/bin/python
# -*- coding: utf-8 -*-
import pandas
from constant_variables import ensembl_classes
from constant_variables import model_species_classes
from constant_variables import define_main_dir


def worker(ensembl_class):

    main_dir = define_main_dir(ensembl_class)
    files_dir = main_dir+'files/'
    
    # Load the file with species abreviations, created in the previous task
    dataset_names_df = pandas.read_csv(files_dir+'ensembl_datasets.tsv',
                                       sep='\t', index_col=0)

    # Create a new tsv file with all the potential names for each abbreviation
    F = open(files_dir+'potential_species_names.tsv', 'w')
    F.write('\tabbreviation\tpotential_names\n')
    for i, abbr in enumerate(dataset_names_df.abbreviation):
        try:        
            if model_species_classes[abbr] == ensembl_class:
                pass
            else:
                continue
        except KeyError:
            pass
        try:
            F.write(str(i)+'\t'+abbr+'\t'+'; '.join(background_names_dict[abbr])+'\n')
        except KeyError:
            F.write(str(i)+'\t'+abbr+'\t'+ 'missing\n') # no name for that abbreviation
    F.close()




if __name__ == "__main__":
    
    '''
    DESCRIPTION:
    
    The previous task returned a tsv file for each ensembl class, containing 
    the name and the version of each dataset, while abbreviation names for the 
    respective species were generated from the strings of dataset names. This 
    abbreviation is important for the following steps, because it will be used 
    as the exclusive id for each species.
    
    Except from the abbreviation name, the main latin name of species should 
    be recorded. However, sometimes an abbreviation corresponds to two or more 
    latin names. This script uses an extra tsv file 
    ('background_mapping_of_species_names_and_abbreviations.tsv'), loaded in 
    'species_names_df' variable, which contains thousands of species names 
    and their abbreviations. Using this mapping, a 'potential names' set will 
    be created for each abbreviation and then, the user need to decide manually 
    the correct latin name for each species, by filtering the output 
    'potential_species_names.tsv' file.
    '''
    
    # Read the background mapping of species names and abbreviations
    tmp = './files/background_mapping_of_species_names_and_abbreviations.tsv'
    species_names_df = pandas.read_csv(tmp, sep='\t', index_col=0)
    
    # Construct a dictionary: {abbreviation: [matched species names]}
    background_names_dict = {}
    for abbr, name in zip(species_names_df.abbreviation, species_names_df.species):
        background_names_dict.setdefault(abbr, []).append(name)
    
    # Execute the 'worker' function for each Ensembl class
    [worker(c) for c in ensembl_classes]
