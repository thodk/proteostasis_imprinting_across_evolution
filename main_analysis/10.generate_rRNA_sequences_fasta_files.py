#!/usr/bin/python
# -*- coding: utf-8 -*-
import pandas
import os
import shutil




if __name__ == '__main__':
    
    '''
    DESCRIPTION:
        
    Copy and paste the fasta files that contain the rRNA sequences from the
    'data' repositories of eukaryotes & prokaryotes
    
    Output: The rRNA_analysis folder which contains another one ('sequences')
    where 409 (the size of species set) fasta files are stored.
    '''
    
    if not os.path.exists('./rRNAs_analysis/sequences/'):
        os.makedirs('./rRNAs_analysis/sequences/')

    species_df = pandas.read_csv('./files/species.tsv', sep='\t')

    # Copy and paste the rRNA fasta files for eukaryotes    
    eukaryotes_df = species_df.loc[species_df.taxonomy == 'eukaryotes', :]
    classes = eukaryotes_df.loc[:,'ensembl_class'].unique()
    fasta_files = []
    for c in classes:
        ref_dir = '../eukaryotes/data/'+c+'/data/rRNA_fasta/'
        fasta_files = fasta_files + [ref_dir+i for i in os.listdir(ref_dir)]
    
    fasta_files = [i for i in fasta_files if os.path.basename(i).split('.')[0] in 
                   eukaryotes_df.rRNA.tolist()]
    for f in fasta_files:
        name = os.path.basename(f)
        shutil.copy(f, './rRNAs_analysis/sequences/'+name)

    # Copy and paste the rRNA fasta files for bacteria
    bacteria_df = species_df.loc[species_df.taxonomy == 'bacteria', :]
    ref_dir = '../prokaryotes/data/bacteria/data/rRNA_fasta/'
    fasta_files = fasta_files + [ref_dir+i for i in os.listdir(ref_dir)]
    fasta_files = [i for i in fasta_files if os.path.basename(i).split('.')[0] in 
                   bacteria_df.rRNA.tolist()]
    for f in fasta_files:
        name = os.path.basename(f)
        shutil.copy(f, './rRNAs_analysis/sequences/'+name)

    # # Copy and paste the rRNA fasta files for archaea
    bacteria_df = species_df.loc[species_df.taxonomy == 'archaea', :]
    ref_dir = '../prokaryotes/data/archaea/data/rRNA_fasta/'
    fasta_files = fasta_files + [ref_dir+i for i in os.listdir(ref_dir)]
    fasta_files = [i for i in fasta_files if os.path.basename(i).split('.')[0] in 
                   bacteria_df.rRNA.tolist()]
    for f in fasta_files:
        name = os.path.basename(f)
        shutil.copy(f, './rRNAs_analysis/sequences/'+name)
