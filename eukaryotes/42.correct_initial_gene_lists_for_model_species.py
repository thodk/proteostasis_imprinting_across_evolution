#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
import pandas
import numpy
from constant_variables import define_main_dir
from constant_variables import model_species_df, model_species


def get_inverse_synonyms(input_file):
    output = {}
    df = pandas.read_csv(input_file, sep='\t')
    primaries = list(df.iloc[:,1])
    synonyms = list(df.iloc[:,2])
    for i in range(len(synonyms)):
        if not pandas.isnull(synonyms[i]):
            tmp_synonyms = synonyms[i].split(' ')
            for entry in tmp_synonyms:
                if not pandas.isnull(entry) and not entry=='' and not entry==';':
                    output.update({entry:primaries[i]})
    return output


def get_translation(input_file):
    translation = {}
    description = {}
    translation_file = pandas.read_csv(input_file, sep='\t')
    gene_symbols = list(translation_file.loc[:,'external_gene_name'])
    ensembl_ids = list(translation_file.loc[:,'ensembl_gene_id'])
    description_strings = list(translation_file.iloc[:,2])
    for i in range(len(gene_symbols)):
        if not pandas.isnull(gene_symbols[i]):
            translation.setdefault('gene_symbols', {}).setdefault(gene_symbols[i], []).append(ensembl_ids[i])
            description.update({gene_symbols[i]:description_strings[i]})
        if not pandas.isnull(ensembl_ids[i]):
            translation.setdefault('ensembl_ids', {}).setdefault(ensembl_ids[i], []).append(gene_symbols[i])
            description.update({ensembl_ids[i]:description_strings[i]})
    return translation, description


def worker(row):

    species = row.species.values[0]
    ensembl_class = row.loc[:,'class'].values[0]
    main_dir = define_main_dir(ensembl_class) +'/annotation/model_species/'
    
    # Load the manually curated gene lists
    initial_data_dir = main_dir+'initial_gene_lists/'
    
    # Load the 'translation' and 'synonyms' files for the gene symbols of 
    # each model species
    additional_data_dir = main_dir+'additional_data/'
    
    # Output directory, where the corrected gene list will be stored    
    corrected_data_dir = main_dir+'corrected_gene_lists/'
    if not os.path.exists(corrected_data_dir):
        os.makedirs(corrected_data_dir)

    # The main part - gene names correction
    inverse_synonyms = get_inverse_synonyms(additional_data_dir+species+'_synonyms.tsv')
    translation, description = get_translation(additional_data_dir+species+'_translation.tsv')
    F = open(initial_data_dir+species+'.txt', 'r')
    bim_genes = [i.replace('\n', '') for i in F.readlines()]

    ver1 = [] # ver1 contains all the possible primary gene symbols
    for gene in bim_genes:
        try:
            translation['gene_symbols'][gene]
            ver1.append(gene)
        except KeyError:
            try:
                ver1.append(inverse_synonyms[gene]) 
            except KeyError:
                try:
                    tmp_gs = translation['ensembl_ids'][gene][0] 
                    if tmp_gs != '' and not pandas.isnull(tmp_gs):
                        ver1.append(tmp_gs)
                    else:
                        ver1.append(numpy.nan)
                except KeyError:
                    ver1.append(numpy.nan)

    ver2 = []
    for gene in ver1:
        try:
           ens = translation['gene_symbols'][gene]
           ens = list(set(ens))
           ver2.append(' '.join(ens))
        except KeyError:
            ver2.append(numpy.nan)

    data={'initial_symbols' : bim_genes,
          'gene_symbols_(primaries)' : ver1,
          'ensembl_ids' : ver2}
    df = pandas.DataFrame(data, columns=['initial_symbols', 
                                         'gene_symbols_(primaries)',
                                         'ensembl_ids' ])
    ds = []
    for i, row in df.iterrows():
        try:
            d = description[row['gene_symbols_(primaries)']]
        except KeyError:
            try:
                d = description[row['ensembl_ids'].split(' ')[0]]
            except (KeyError, AttributeError):
                d = numpy.nan
        ds.append(d)
    df.loc[:,'description'] = ds
    bool1 = pandas.isnull(df.loc[:,'gene_symbols_(primaries)'])
    bool2 = pandas.isnull(df.loc[:,'ensembl_ids'])
    df = df.loc[~(bool1 | bool2)]
    df.to_csv(corrected_data_dir+species+'.tsv', sep='\t', index=False)




if __name__ == "__main__":
    
    '''
    DESCRIPTION:

    Corretion of the initial gene lists, using the retrieved mappings from the 
    previous steps (41A and 41B). The final proteostasis related gene 
    sets will be stored in the 
    '/data/class/preprocessing/annotation/model_species/corrected_gene_lists/'
    folder of each eukaryotic class.
    '''
    
    # Execute the 'worker' function for each model species
    for species in model_species:
        row = model_species_df.loc[model_species_df.species == species,:]
        worker(row)






