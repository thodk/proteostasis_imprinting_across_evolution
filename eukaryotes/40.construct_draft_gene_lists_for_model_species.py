#!/usr/bin/python
# -*- coding: utf-8 -*-
import pandas
import os



if __name__ == "__main__":
    
    '''
    DESCRIPTION:

    This script reads the initial_table_evol_proteostasis_manual_curation.csv' 
    file (stored in 'files') and creates a gene list for each model species. 
    The input file was created manually and it does not have correct tabular 
    format. The script retrieves the necessary data for each model species and 
    creates a file in the directory
    '/data/class/preprocessing/annotation/model_species/initial_gene_lists/ for
    each ensembl class.

    These initial (draft) lists will be corrected in the following steps, in 
    order to include only primary gene names (gene symbols). The final corrected 
    lists will be used as reference to retrieve gene homologies for all the 
    other species of the same taxonomic class and create the proteostasis 
    related gene sets for them.
    '''
    
    model_species_df = pandas.read_csv('./files/model_species.tsv', 
                                       index_col=0, sep='\t')
    indexes = [i+1 for i in model_species_df.index.tolist()]
    species_dict = dict(zip(indexes, model_species_df.species.tolist()))
    
    # 1. Read the Table_Evol_Proteostasis_manual_curation file which contains
    # manually curated lists of genes for model species and prokaryotes
    # 2. Save the data in a dictionary, with model species names as keys    
    data = {}
    proteostasis_file = open('./files/initial_table_evol_proteostasis_manual_curation.csv', 'r')
    next(proteostasis_file)
    for line in proteostasis_file:
        fields = line.split(',')
        for i in range(len(fields)):
            try:
                species_dict[i]
                if fields[i] != '' and not pandas.isnull(fields[i]):
                    gene = fields[i].replace(' ', '')
                    data.setdefault(species_dict[i], []).append(gene)
            except KeyError:
                pass
    proteostasis_file.close()
    
    # 1. For each model species, create a file with proteostasis-relates genes
    # 2. Save it in '/data/class/preprocessing/annotation/model_organisms/initial_datasets/'
    for model_species, genes in data.items():
        row = model_species_df.loc[model_species_df.species == model_species,:]
        mart = row.mart.values[0]
        ensembl_class = row.loc[:,'class'].values[0]
        
        p = '/preprocessing/annotation/model_species/initial_gene_lists/'
        main_dir = './data/'+ensembl_class+p
        
        if not os.path.exists(main_dir):
            os.makedirs(main_dir)
    
        F = open(main_dir+model_species+'.txt', 'w')
        F.write('\n'.join(list(set(genes))))
        F.close()