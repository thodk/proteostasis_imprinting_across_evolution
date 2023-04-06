#!/usr/bin/python
# -*- coding: utf-8 -*-
import pandas
import os
import sys
from constant_variables import define_main_dir
from constant_variables import prokaryotic_classes
sys.path.append('../')
from core_functions import uniprot_query


def worker(prokaryotic_class, hsp):

    main_dir = define_main_dir(prokaryotic_class)
    hsp_dir = main_dir+'data/'+hsp+'_fasta/'
    if not os.path.exists(hsp_dir):
        os.makedirs(hsp_dir)
    
    species_df = pandas.read_csv(main_dir+'preprocessing/03_species_with_rRNA_and_HSPs_annotation.tsv', 
                                   sep='\t', index_col=0)
    columns = ['id', 'sequence']
    
    # For each species (i.e. abbreviation name) there is a string for the 
    # annotated HSPs (hsp_string). The functions inside the loop split that
    # script and send a query to the UniProt database in order to get the
    # respective sequences.
    for abbreviation, hsp_string in zip(species_df.abbreviation, species_df.loc[:, hsp]):
        
        hsp_ids = hsp_string.split(', ')      
        query = '+OR+'.join(['id:'+i for i in hsp_ids])        
        success, data, column_names = uniprot_query(query, columns)
        
        if len(data) == 0:
            continue
        else:
            entries = []
            for s in data[1:]:
                fields = s.split('\t')
                if len(fields) < len(column_names):
                    continue
                entries.append(fields)
            
            # Fasta file creation
            tmp_df = pandas.DataFrame(entries, columns=columns)
            F = open(hsp_dir+abbreviation+'.fasta', 'w')
            for i, row in tmp_df.iterrows():
                if i == 0:
                    F.write('>' + row.id + '\n' + row.sequence)
                else:
                    F.write('\n>' + row.id + '\n' + row.sequence)            
            F.close()




if __name__ == "__main__":

    '''
    DESCRIPTION:
    
    This script collects the HSP sequences from the UniProt database. It reads
    the '03_species_with_rRNA_and_HSPs_annotation.tsv' file to get the protein
    id for all the HSPs and then sends the query to retrieve the sequence data. 
    The final files are in fasta format, containing the protein id and sequence 
    string for each entry.
    
    All the files are stored in '/data/class/data/'hsp40_fasta/' and
    '/data/class/data/'hsp70_fasta/'.
    
    Species are named with a unique id, like in the case of eukaryotes, in order 
    to avoid to use their abbreviations during the analysis. This nomenclature 
    is safer and provides a more efficient way to create the respective MongoDB 
    collections:
    
    example: bacteria_0, bacteria_10, archaea_20 etc
    '''

    # Execute the 'worker' function for each prokaryotic class and each HSPs
    # family
    for c in prokaryotic_classes:
        worker(c, 'hsp40')
        worker(c, 'hsp70')
     
        

