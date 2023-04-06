#!/usr/bin/python
# -*- coding: utf-8 -*-
import pandas
import os
import shutil


def collect_data_for_eukaryotes(hsp):

    eukaryotes_df = species_df.loc[species_df.taxonomy == 'eukaryotes', :]
    classes = eukaryotes_df.loc[:,'ensembl_class'].unique()

    fasta_files = []
    for c in classes:    
        ref_dir = '../eukaryotes/data/'+c+'/data/'+hsp+'_fasta/'
        fasta_files = fasta_files + [ref_dir+i for i in os.listdir(ref_dir)]
    
    fasta_files = [i for i in fasta_files if os.path.basename(i).split('.')[0] in 
                   eukaryotes_df.abbreviation.tolist()]

    for f in fasta_files:
        name = os.path.basename(f)
        shutil.copy(f, './HSPs_analysis/'+hsp+'/sequences/'+name)


def collect_data_for_prokaryotes(category, hsp):
    main_dir = '../prokaryotes/data/'+category+'/data/'+hsp+'_fasta/'
    fasta_files = [main_dir+i for i in os.listdir(main_dir)]
    for f in fasta_files:
        name = os.path.basename(f)
        shutil.copy(f, './HSPs_analysis/'+hsp+'/sequences/'+name)



        
if __name__ == '__main__':
    
    '''
    DESCRIPTION:
        
    Copy and paste the fasta files that contain the HSP sequences from the
    'data' repositories of eukaryotes & prokaryotes.
    
    Output: The HSPs_analysis folder which contains two main subfolders, one 
    for HSP40s and one for HSP70s.
    '''
    
    if not os.path.exists('./HSPs_analysis/hsp40/sequences/'):
        os.makedirs('./HSPs_analysis/hsp40/sequences/')
        
    if not os.path.exists('./HSPs_analysis/hsp70/sequences/'):
        os.makedirs('./HSPs_analysis/hsp70/sequences/')
    
    species_df = pandas.read_csv('./files/species.tsv', sep='\t')
    
    collect_data_for_eukaryotes('hsp40')
    collect_data_for_eukaryotes('hsp70')
    collect_data_for_prokaryotes('bacteria', 'hsp40')
    collect_data_for_prokaryotes('bacteria', 'hsp70')
    collect_data_for_prokaryotes('archaea', 'hsp40')
    collect_data_for_prokaryotes('archaea', 'hsp70')
