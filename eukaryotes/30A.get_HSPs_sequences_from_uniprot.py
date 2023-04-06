#!/usr/bin/python
# -*- coding: utf-8 -*-
import pandas
import os
import sys
sys.path.append('../')
import core_functions
from constant_variables import define_main_dir
from constant_variables import ensembl_classes


def worker(ensembl_class, hsp):
    
    main_dir = define_main_dir(ensembl_class)
    output_dir = main_dir+'/uniprot_data/'+hsp+'_sequences/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Load the file which contains the species which have passed the previous
    # filtering steps (HSPs and rRNA)
    species_names_df = pandas.read_csv(main_dir+'03_species_with_rRNA_and_HSPs_annotation.tsv', 
                                   sep='\t', index_col=0)
    
    # The query attributes
    columns = ['id', 'genes(PREFERRED)', 'sequence']
    # Alternative names for the query arrtibutes
    df_columns = ['id', 'primary_gene_names', 'sequence']

    for i, entry in species_names_df.iterrows():
        f = main_dir+'filtered_HSPs_annotation/'+entry.species_abbreviation+'.tsv'
        df = pandas.read_csv(f, sep='\t', index_col=0)
        df = df.loc[((df.hsp==hsp) & (df.resource=='uniprot')),:]
        if len(df) == 0:
            uniprot_df = pandas.DataFrame(dict((i,[]) for i in df_columns))
        else:
            query = '+OR+'.join(['id:'+i for i in df.gene_symbol.tolist()])
            success, data, column_names = core_functions.uniprot_query(query, columns)
            # If no data are returned create an empty data frame
            if len(data) == 0 or data == ['']: 
                uniprot_df = pandas.DataFrame(dict( (i,[]) for i in df_columns ))
            # If data are returned, check their validity and create a data frame
            else:  
                entries = []
                for s in data[1:]:
                    fields = s.split('\t')
                    if len(fields) < len(column_names):
                        continue
                    entries.append(fields)
                uniprot_df = pandas.DataFrame(entries)
                uniprot_df.columns = df_columns
        
        # 1. Save the data as tsv and not fasta files
        # 2. Use as file name the id of each species, e.g. plant_0, fungi_8
        uniprot_df.to_csv(output_dir+entry.abbreviation+'.tsv', sep='\t')
      



if __name__ == "__main__":
    
    '''
    DESCRIPTION:

    This script retrieves the HSP sequences from the UniProt database. The final 
    files are in tsv format, containing the protein id, gene symbol and sequence 
    for each entry.

    All the files are stored in the '/uniprot_data/hsp40_sequences/' and 
    '/uniprot_data/hsp70_sequences/'subfolder of each ensembl class directory.

    Species are named with a unique id, based on their 'ensembl_class', in order 
    to avoid to use their abbreviations during the analysis. This nomenclature 
    is simpler and provides a safer way to create the necessary Mongo collections 
    in the following steps, insted of using species names:

    example: ensembl_0, ensembl_10, ensembl_20, plant_2, plant_13 etc 
    '''
    
    # Execute the 'worker' function for each Ensembl class and HSP40
    [worker(c, 'hsp40') for c in ensembl_classes]
    # Execute the 'worker' function for each Ensembl class and HSP70
    [worker(c, 'hsp70') for c in ensembl_classes]    
 
