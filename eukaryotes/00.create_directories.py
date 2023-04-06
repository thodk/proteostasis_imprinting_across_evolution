#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas
import os



if __name__ == '__main__':
    
    df = pandas.read_csv('./files/ensembl_parameters.tsv', sep='\t')
    eukaryotic_classes = df.loc[:, 'class']
    for c in eukaryotic_classes:
        os.makedirs('./data/'+c+'/data/rRNA_fasta/', exist_ok=True)
        os.makedirs('./data/'+c+'/data/hsp40_fasta/', exist_ok=True)
        os.makedirs('./data/'+c+'/data/hsp70_fasta/', exist_ok=True)
        os.makedirs('./data/'+c+'/preprocessing/files/', exist_ok=True)