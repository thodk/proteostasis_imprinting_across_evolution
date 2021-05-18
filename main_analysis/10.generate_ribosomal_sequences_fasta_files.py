# -*- coding: utf-8 -*-

import pandas
import os
import shutil



if __name__ == '__main__':
    
    if not os.path.exists('./rRNAs_analysis/sequences/'):
        os.makedirs('./rRNAs_analysis/sequences/')

    # eukaryotes
    species_df = pandas.read_csv('./files/species.tsv', sep='\t')
    eukaryotes_df = species_df.loc[species_df.taxonomy == 'eukaryotes', :]
    classes = eukaryotes_df.loc[:,'ensembl_mart'].unique()
    fasta_files = []
    for c in classes:
        ref_dir = '../eukaryotes/data/'+c+'/data/rRNA_fasta/'
        fasta_files = fasta_files + [ref_dir+i for i in os.listdir(ref_dir)]
    
    fasta_files = [i for i in fasta_files if os.path.basename(i).split('.')[0] in 
                   eukaryotes_df.rRNA.tolist()]
    for f in fasta_files:
        name = os.path.basename(f)
        shutil.copy(f, './rRNAs_analysis/sequences/'+name)


    # bacteria
    bacteria_df = species_df.loc[species_df.taxonomy == 'bacteria', :]
    ref_dir = '../prokaryotes/bacteria/data/rRNA_fasta/'
    fasta_files = fasta_files + [ref_dir+i for i in os.listdir(ref_dir)]
    fasta_files = [i for i in fasta_files if os.path.basename(i).split('.')[0] in 
                   bacteria_df.rRNA.tolist()]
    for f in fasta_files:
        name = os.path.basename(f)
        shutil.copy(f, './rRNAs_analysis/sequences/'+name)

    # archaea
    bacteria_df = species_df.loc[species_df.taxonomy == 'archaea', :]
    ref_dir = '../prokaryotes/archaea/data/rRNA_fasta/'
    fasta_files = fasta_files + [ref_dir+i for i in os.listdir(ref_dir)]
    fasta_files = [i for i in fasta_files if os.path.basename(i).split('.')[0] in 
                   bacteria_df.rRNA.tolist()]
    for f in fasta_files:
        name = os.path.basename(f)
        shutil.copy(f, './rRNAs_analysis/sequences/'+name)
