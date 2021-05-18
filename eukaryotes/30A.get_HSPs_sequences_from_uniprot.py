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

    # load the file which contains the species which have passed the previous
    # filtering steps (HSPs and rRNA)
    species_names_df = pandas.read_csv(main_dir+'03_species_with_rRNA_and_HSPs_annotation.tsv', 
                                   sep='\t', index_col=0)
    
    # query for protein, gene ids and amino acid sequence
    columns = ['id', 'genes(PREFERRED)', 'sequence']
    df_columns = ['id', 'primary_gene_names', 'sequence']

    for i, entry in species_names_df.iterrows():
        f = main_dir+'filtered_hsp_annotation/'+entry.species_abbreviation+'.tsv'
        df = pandas.read_csv(f, sep='\t', index_col=0)
        df = df.loc[((df.hsp==hsp) & (df.resource=='uniprot')),:]
        if len(df) == 0:
            uniprot_df = pandas.DataFrame(dict((i,[]) for i in df_columns))
        else:
            query = '+OR+'.join(['id:'+i for i in df.gene_symbol.tolist()])
            success, data, column_names = generic_functions.uniprot_query(query, columns)
            # if no data was returned create an empty data frame
            if len(data) == 0 or data == ['']: 
                uniprot_df = pandas.DataFrame(dict( (i,[]) for i in df_columns ))
            # if data was returned, check their validity and create a data frame
            else:  
                entries = []
                for s in data[1:]:
                    fields = s.split('\t')
                    if len(fields) < len(column_names):
                        continue
                    entries.append(fields)
                uniprot_df = pandas.DataFrame(entries)
                uniprot_df.columns = df_columns
        # save the data as tsv and not fasta files
        # use as file name the specific abbreviation code of each organism,
        # e.g. plant_0, plant_10, fungi_8 etc
        uniprot_df.to_csv(output_dir+entry.abbreviation+'.tsv', sep='\t')
      

'''
DESCRIPTION:

That script retrieves the HSP sequences from the UniProt database. The final files
are in tsv format, containing the protein id, gene symbol and sequence for each
entry.

All the files are stores in '/uniprot_data/'+hsp+'_sequences/'

Species are named with a unique id, based on their ensembl_class, in order to
avoid to use their abbreviations, during the analysis. That nomenclature is simpler
and provides a more efficient way to create Mongo collections using species names:

example: ensembl_0, ensembl_10, ensembl_20, plant_2, plant_13 etc 
'''

if __name__ == "__main__":
    [worker(c, 'hsp40') for c in ensembl_classes]
    [worker(c, 'hsp70') for c in ensembl_classes]    
 
