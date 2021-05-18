# -*- coding: utf-8 -*-

import pandas
import sys
import os
sys.path.append('../')
from core_functions import uniprot_query
from core_functions import correct_uniprot_dataframe
from core_functions import str_modification
from core_functions import filter_hsp_data
from constant_variables import ensembl_classes
from constant_variables import define_main_dir


def worker(ensembl_class):
    
    main_dir = define_main_dir(ensembl_class)
    files_dir = main_dir+'files/'
    hsp40_dir = main_dir+'uniprot_data/hsp40_raw_files/'
    hsp70_dir = main_dir+'uniprot_data/hsp70_raw_files/'
    
    # create hsp40_dir & hsp70_dir to save the retrieved HSPs annotation
    if not os.path.exists(hsp40_dir):
        os.makedirs(hsp40_dir)
        
    if not os.path.exists(hsp70_dir):
        os.makedirs(hsp70_dir)

    # load the tsv file with organism names and abbreviations
    species_names_df = pandas.read_csv(files_dir+'final_species_names.tsv', sep='\t', index_col=0)

    # 1. defined the features that will be retrieved from the uniprot database 
    # for each entry
    # 2. alterantive names for uniprot features, to use them as column names 
    # in the final output tsv files
    query_columns = ['id', 'reviewed', 'organism', 'genes(PREFERRED)', 'protein names']
    df_columns = ['id', 'reviewed', 'species_name', 'primary_gene_names', 'protein_name']

    for i, entry in species_names_df.iterrows():

        # retrieve data using the uniprot API, criteria: query and columns
        query = 'organism:'+entry.species_name 
        success, data, column_names = uniprot_query(query, query_columns)
        
        if len(data) == 0 or data == ['']: # empty response
            continue
        else: # non-empty response
            entries = []
            for e in data[1:]:
                entries.append(e.split('\t')) # split the entries with \t
            uniprot_tmp_df = pandas.DataFrame(entries)
            uniprot_tmp_df.columns = df_columns
            # use the 'correct_uniprot_dataframe' function to duplicate an 
            # entry if it contains more that one primary gene names in the 
            # respective field
            uniprot_final_df = correct_uniprot_dataframe(uniprot_tmp_df, 
                                                         ['primary_gene_names'],
                                                         'primary_gene_names')
            # 1. filter entries to keep only those whose primary gene name
            # is related to HSP40 or HSP70
            # 2. filter entries to keep only those for the specific organism
            # (sometimes uniprot returns proteins which belong to other species)
            # 3. store the data in tsv format in /preprocessing/uniprot_hsp40(or 70)_raw_files/
            tmp = str_modification(entry.species_name)
            
            hsp70_df = filter_hsp_data(uniprot_final_df, 'hsp70')
            hsp70_df = hsp70_df.loc[hsp70_df.species_name.str.contains(tmp, regex=True), :]
            hsp70_df.to_csv(hsp70_dir+entry.abbreviation+'.tsv', sep='\t')

            hsp40_df = filter_hsp_data(uniprot_final_df, 'hsp40')
            hsp40_df = hsp40_df.loc[hsp40_df.species_name.str.contains(tmp, regex=True), :]
            hsp40_df.to_csv(hsp40_dir+entry.abbreviation+'.tsv', sep='\t')
        

'''
DESCRIPTION:

That script retrieves annotation for HSPs from the UniProt database. The worker 
function is executed for the species of each ensembl class separately. The final 
goal is to detect which species have adequate annotation for HSPs and filter out 
all the others. The output lists of annotated HSPs will be combined with similar 
data retrieved from the Ensembl database (see 10B.get_HSPs_annotation_from_ensembl.R)

The request process is imported from the core_functions module. Also, some
additional functions of core_functions are applied to filter the retireved
data and keep only trustworthy annotations.

'''

if __name__ == "__main__":
    # execute the worker for each ensembl class to create the 
    # 'species_potential_names.tsv' file
    [worker(c) for c in ensembl_classes]









  
  
  
    
