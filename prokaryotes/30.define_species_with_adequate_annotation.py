#!/usr/bin/python
# -*- coding: utf-8 -*-
import pandas
from constant_variables import define_main_dir
from constant_variables import prokaryotic_classes


def worker(prokaryotic_class):

    main_dir = define_main_dir(prokaryotic_class)+'preprocessing/'
    species_df = pandas.read_csv(main_dir+'02_species_with_rRNA_annotation.tsv', 
                                 sep='\t', index_col=0)
    
    # Keep only entries with rRNA data 
    species_df = species_df.loc[~pandas.isnull(species_df.rRNA), :] 
    
    # Define species latin abbreviations
    species_dict = {}
    for species in species_df.species.tolist():
        fields = species.split(' (')[0].split(' ')
        if len(fields) == 2:
            species_abbr = fields[0][0].capitalize()+''.join(fields[0][1:]) + ' ' + fields[1]   
        else:
            species_abbr = fields[0][0].capitalize()+''.join(fields[0][1:]) + ' ' + ' '.join(fields[1:])  
        species_dict.setdefault(species_abbr, []).append(species)
    
    # Correct species names taking into account the strains
    species_names_dict = {}
    for species_abbr, species_list in species_dict.items():
        for species in species_list:
            species_name = species.split('(')[0]
            try:
                strain = species.split('(')[1]
                strain = strain.replace('strain ', '').split('/')[0].replace(')', '')
            except (ValueError, IndexError):
                strain = None

            species_name = species_abbr
            if strain != None:
                species_name = species_name + ' ' + strain
            else:
                pass
        
            species_name = species_name.replace('.', '')
            species_names_dict.update({species:species_name})
    
    # Create the new file for the annotated species
    species_df.loc[:, 'species_name'] = [species_names_dict[s] for s in species_df.species.tolist()]
    species_df = species_df.reset_index(drop=True)
    species_df.loc[:, 'abbreviation'] = [prokaryotic_class+'_'+str(i) for i in range(species_df.shape[0])]
    species_df.to_csv(main_dir+'03_species_with_rRNA_and_HSPs_annotation.tsv', sep='\t')




if __name__ == "__main__":
    
    '''
    DESCRIPTION:

    This script filters out species without HSPs and rRNA annotation from the 
    '02_species_with_rRNA_annotation.tsv' file and creates the
    '03_species_with_rRNA_and_HSPs_annotation' for each prokaryotic class.
    
    Then, a unified data frame for bacteria and archaea is created 
    ('final_species_set_for_analysis.tsv') in the 'files' folder.

    The following steps will use this set of species to create the GO annotation 
    files and MongoDB collections.
    '''
    
    # Execute the 'worker' function for each prokaryotic class
    [worker(c) for c in prokaryotic_classes]
    
    

    global_df = pandas.DataFrame(data={'abbreviation':[], 'prokaryotes_class':[]})

    for c in ['bacteria', 'archaea']:
        main_dir = define_main_dir(c)+'preprocessing/'
        filename = main_dir+'03_species_with_rRNA_and_HSPs_annotation.tsv'
        df = pandas.read_csv(filename, sep='\t', index_col=0)
        df.loc[:,'prokaryotes_class'] = [c]*len(df)
        global_df = pandas.concat([global_df, df], axis=0, ignore_index=True, sort=True)
        global_df = global_df.reset_index(drop=True)
        global_df.to_csv('./files/final_species_set_for_analysis.tsv', sep='\t')
    
    
