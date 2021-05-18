# -*- coding: utf-8 -*-

import pandas
from constant_variables import ensembl_classes
from constant_variables import model_species_classes
from constant_variables import define_main_dir


def worker(ensembl_class):

    main_dir = define_main_dir(ensembl_class)
    files_dir = main_dir+'files/'
    # load the file with species abreviations, created in the previous task
    dataset_names_df = pandas.read_csv(files_dir+'ensembl_datasets.tsv',
                                       sep='\t', index_col=0)

    # create a new tsv file with all the potential names for each abbreviation
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


'''
DESCRIPTION:

The previous task returned a tsv file for each ensembl class, containing the 
name and the version of each dataset, while abbreviation names for the respective 
species were generated from the strings of dataset names. That abbreviation
is important for the following steps, because it will be used as the exclusive 
id for the species.

Except from the abbreviation, the latin name of species should be recorded.
However, sometimes an abbreviation corresponds to two or more latin names. 
That current script uses an extra tsv file, loaded in 'species_names_df'
variable, which contains thousands of species names and their abbreviations.
Using that mapping, a 'potential names' set will be created for each abbreviation
and then, the user should decide the correct latin name, filtering manually 
the 'potential_species_names.tsv' file.

'''

if __name__ == "__main__":

    # 1. read the background mapping of species names and abbreviations
    # 2. construct a dict: {abbreviation: [matched species names]}
    tmp = './files/background_mapping_of_species_names_and_abbreviations.tsv'
    species_names_df = pandas.read_csv(tmp, sep='\t', index_col=0)
    background_names_dict = {}
    for abbr, name in zip(species_names_df.abbreviation, species_names_df.species):
        background_names_dict.setdefault(abbr, []).append(name)

    # execute the worker for each ensembl class to create the 
    # 'species_potential_names.tsv' file
    [worker(c) for c in ensembl_classes]
