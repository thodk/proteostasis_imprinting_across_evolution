#!/usr/bin/python
# -*- coding: utf-8 -*-
import pandas
import os
from constant_variables import define_main_dir
from constant_variables import ensembl_classes


def collect_data(main_dir, organism, hsp):
    
    # Uniprot
    data = []    
    try:
        df = pandas.read_csv(main_dir+'uniprot_data/'+hsp+'_sequences/'+organism+'.tsv', 
                             sep='\t', index_col=0)
        zip_tup = zip(df.id.tolist(), df.primary_gene_names.tolist(), df.sequence.tolist())
        for trip in zip_tup:
            if pandas.isnull(trip[2]):
                continue
            data.append([str(trip[0])+'_'+str(trip[1]) , trip[2].replace('*', '')])
    except (AttributeError, IOError):
        pass
    
    # Ensembl
    try:
        df = pandas.read_csv(main_dir+'ensembl_data/'+hsp+'_sequences/'+organism+'.tsv', 
                             sep='\t', index_col=0)
        zip_tup = zip(df.index.tolist(), df.gene_symbol.tolist(), df.peptide.tolist())
        for trip in zip_tup:
            if pandas.isnull(trip[2]) or pandas.isnull(trip[1]):
                continue
            data.append([str(trip[0])+'_'+str(trip[1]) , trip[2].replace('*', '')])
    except (AttributeError, IOError):
        pass
    
    return data


def worker(ensembl_class):

    main_dir = define_main_dir(ensembl_class)
    hsp40_dir = './data/'+ensembl_class+'/data/hsp40_fasta/'
    hsp70_dir = './data/'+ensembl_class+'/data/hsp70_fasta/'
    species_df = pandas.read_csv(main_dir+'03_species_with_rRNA_and_HSPs_annotation.tsv',
                                   sep='\t')
    
    # Create the folders to save the final fasta files
    if not os.path.exists(hsp40_dir):
        os.makedirs(hsp40_dir)
    if not os.path.exists(hsp70_dir):
        os.makedirs(hsp70_dir)    
    
    # For each species, collect the valid sequences of HSP40 and HSP70 from
    # UniProt and Ensembl subfolders and create a unified fasta file
    for organism in species_df.abbreviation.tolist():
        data = collect_data(main_dir, organism, 'hsp70')
        if len(data) != 0: 
            F = open(hsp70_dir+organism+'.fasta', 'w')
            for entry in data:
                if entry[1] == 'Sequence unavailable':
                    continue
                F.write('>'+entry[0]+'\n')
                F.write(entry[1]+'\n')
            F.close()
        
        data = collect_data(main_dir, organism, 'hsp40')
        if len(data) != 0: 
            F = open(hsp40_dir+organism+'.fasta', 'w')
            for entry in data:
                if entry[1] == 'Sequence unavailable':
                    continue
                F.write('>'+entry[0]+'\n')
                F.write(entry[1]+'\n')
            F.close()




if __name__ == "__main__":
    
    '''
    DESCRIPTION:

    This script collects the sequences which have been retrieved from UniProt 
    and Ensembl and creates the final fasta files. The files are stored in the 
    'data/class/data/hsp40/' and /'data/class/data/hsp70' directories of each 
    ensembl class. They constitute the final dataset that will be used for the
    analysis of HSP sequences.
    '''
    
    # Execute the 'worker' function for each ensembl class
    [worker(c) for c in ensembl_classes]



  
    
    
