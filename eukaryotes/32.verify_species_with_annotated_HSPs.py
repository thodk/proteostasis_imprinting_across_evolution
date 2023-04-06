#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import pandas
from constant_variables import define_main_dir
from constant_variables import ensembl_classes


def worker(ensembl_class):

    main_dir = define_main_dir(ensembl_class)
    hsp40_dir = './data/'+ensembl_class+'/data/hsp40_fasta/'
    hsp70_dir = './data/'+ensembl_class+'/data/hsp70_fasta/'
    hsp40_fasta = os.listdir(hsp40_dir)
    hsp70_fasta = os.listdir(hsp70_dir)
    
    species_df = pandas.read_csv(main_dir+'03_species_with_rRNA_and_HSPs_annotation.tsv',
                                 index_col=0, sep='\t')
    to_remove = []
    for abbreviation in species_df.abbreviation.to_list():
        key = abbreviation+'.fasta'
        bool1 = key in hsp40_fasta
        bool2 = key in hsp70_fasta
        if not (bool1 and bool2):
            to_remove.append(abbreviation)
            
    if len(to_remove) == 0:
        return
    else:
        tmp = species_df.loc[species_df.loc[:, 'abbreviation'].isin(to_remove)]
        species_df.drop(tmp.index, inplace=True)
        species_df.to_csv(main_dir+'03_species_with_rRNA_and_HSPs_annotation.tsv',
                          sep='\t')
       
        
        

if __name__ == '__main__':
    
    # Execute the 'worker' function for each ensembl class
    [worker(c) for c in ensembl_classes]
