#!/usr/bin/python
# -*- coding: utf-8 -*-
import pandas
import os
import shutil


def worker_for_prokaryotes(category):
    
    main_dir = '../prokaryotes/data/'+category+'/data/PN_gene_lists/'
    files = [main_dir+f for f in os.listdir(main_dir)]
    for f in files:
        new_f = f.split('/')[-1]
        shutil.copy(f, './PN_analysis/pathway_analysis_inputs/'+new_f)

def worker_for_eukaryotes(category):
    
    main_dir = '../eukaryotes/data/'+category+'/data/PN_gene_lists/'
    files = [main_dir+f for f in os.listdir(main_dir)]
    for f in files:
        new_f = f.split('/')[-1]
        shutil.copy(f, './PN_analysis/pathway_analysis_inputs/'+new_f)   



if __name__ == '__main__':

    '''
    DESCRIPTION:
        
    This script creates the input files for pathway analysis.
    
    '''
    
    species_df = pandas.read_csv('./files/species.tsv', sep='\t')
    if not os.path.exists('./PN_analysis/pathway_analysis_inputs'):
        os.makedirs('./PN_analysis/pathway_analysis_inputs')
    
    ensembl_parameters_df = pandas.read_csv('../eukaryotes/files/ensembl_parameters.tsv',
                                            sep='\t', index_col=0)
    for c in ensembl_parameters_df.loc[:,'class'].tolist():
        worker_for_eukaryotes(c)
    worker_for_prokaryotes('bacteria')
    worker_for_prokaryotes('archaea')

