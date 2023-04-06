#!/usr/bin/python
# -*- coding: utf-8 -*-
import pandas
import os
import sys
sys.path.append('../')
from core_functions import uniprot_query
from core_functions import correct_uniprot_dataframe
from constant_variables import define_main_dir
from constant_variables import model_species_df, model_species


def worker(row):
    species = row.species.values[0]
    ensembl_class = row.loc[:,'class'].values[0]
    
    main_dir = define_main_dir(ensembl_class)
    output_dir = main_dir+'/annotation/model_species/additional_data/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    species_names_df = pandas.read_csv(main_dir+'03_species_with_rRNA_and_HSPs_annotation.tsv', 
                                   sep='\t', index_col=0)
    tmp = species_names_df[species_names_df.species_abbreviation == species]
    species_name = tmp.species_name.values[0]
    print(species_name)
    
    # Send the query to uniprot and get the gene symbols and alternative names
    query_columns = ['genes(PREFERRED)', 'genes(ALTERNATIVE)']
    df_columns = ['primary_gene_names', 'alternative_gene_names']
    query = 'organism:'+species_name
    success, data, column_names = uniprot_query(query, query_columns)
    if len(data) == 0 or data == ['']:
        uniprot_df = pandas.DataFrame(dict( (i,[]) for i in df_columns ))
    else:
        entries = []
        for s in data[1:]:
            fields = s.split('\t')
            if len(fields) < len(column_names):
                continue
            entries.append(fields)
        uniprot_df = pandas.DataFrame(entries)
        uniprot_df.columns = df_columns

    uniprot_final_df = correct_uniprot_dataframe(uniprot_df, df_columns,
                                                 'primary_gene_names')
    uniprot_final_df.to_csv(output_dir+species+'_synonyms.tsv', sep='\t')




if __name__ == "__main__":

    '''
    DESCRIPTION:
    
    This script calls a request function to send a query to the UniProt database 
    in order to retrieve alternative gene names for each gene symbol of model
    species. This task is necessary, as the manually curated gene lists of 
    proteostasis contain various types of gene names (primary and secondary 
    symbols, ensembl ids, entrez ids etc). All of them should be translated 
    into primary gene symbols.
    '''

    # Execute the 'worker' function for each model species
    for species in model_species:
        row = model_species_df.loc[model_species_df.species == species,:]
        worker(row)


