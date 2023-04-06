#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas
import pymongo
import os
import sys
sys.path.append('../')
from constant_variables import define_main_dir
from constant_variables import model_species_df, model_species
from core_functions import mongodb


def worker(row):
    
    ensembl_class = row.loc[:,'class'].values[0]
    preprocessing_dir = define_main_dir(ensembl_class)
    results_dir = os.path.abspath('./data/'+ensembl_class+'/data/PN_gene_lists')+'/'
    
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    df = pandas.read_csv(preprocessing_dir+'04_species_with_genomic_annotation.tsv', sep='\t', index_col=0) 
    species_dict = dict(zip(df.species_abbreviation.tolist(), df.abbreviation.tolist()))

    ref_dir = preprocessing_dir+'annotation/model_species/corrected_gene_lists/'
    reference_files = os.listdir(ref_dir)
    # Create the reference lists for the proteostasis-related genes of
    # model species
    reference_mapping = {}
    for reference_file in reference_files:
        reference_species = reference_file.split('.')[0]        
        df = pandas.read_csv(ref_dir+reference_file, sep='\t')
        if reference_species == 'xtropicalis': # if xtropicalis get the ensembl ids
            genes = df.loc[:,'ensembl_ids'].tolist()
            genes = [j for j in genes if not pandas.isnull(j)]
            genes = [i for j in genes for i in str(j).split(' ')]
        else: # else get the gene symbols
            genes = df.loc[:,'gene_symbols_(primaries)'].tolist()
        genes = list(filter(lambda x: not pandas.isnull(x), genes)) # remove NANs
        reference_mapping.update({reference_species:genes})
    
    # Get the homologies mappings and create the gene sets for all the other species
    # given their reference model species.
    homo_dir = preprocessing_dir+'/annotation/homologies/'
    for species, abbreviation in species_dict.items():
        if species in reference_eukaryotes:
            continue
        species_genes = []
        for ref_species, ref_genes in reference_mapping.items():
            tmp_df = pandas.read_csv(homo_dir+species+'_'+ref_species+'_homologs.tsv',
                                     sep='\t', index_col=0)
            tmp_df = tmp_df.loc[tmp_df.external_gene_name.isin(ref_genes), :]
            if ensembl_class == 'ensembl':
                homologies = tmp_df.loc[:,species+'_homolog_ensembl_gene'].tolist()
            else:
                try:
                    homologies = tmp_df.loc[:,species+'_eg_homolog_ensembl_gene'].tolist()
                except KeyError:
                    homologies = tmp_df.loc[:,species+'_eg_homolog_associated_gene_name'].tolist()
            homologies = list(filter(lambda x: not pandas.isnull(x), homologies))
            species_genes = species_genes + homologies
        species_genes = list(set(species_genes))
        F = open(results_dir+abbreviation+'.txt', 'w')
        F.write('\n'.join(species_genes))
        F.close()

    # Create the files for model species
    for ref_species, ref_genes in reference_mapping.items():
        F = open(results_dir+species_dict[ref_species]+'.txt', 'w')
        F.write('\n'.join(ref_genes))
        F.close()
    



if __name__ == '__main__':

    '''
    DESCRIPTION:
        
    This script creates the input files for pathway analysis.
    
    '''
    
    ensembl_parameters_df = pandas.read_csv('./files/ensembl_parameters.tsv',
                                            sep='\t', index_col=0)
    reference_eukaryotes = []
    for c in ensembl_parameters_df.loc[:,'class'].tolist():
        ref_dir = '../eukaryotes/data/'+c+'/preprocessing/annotation/model_species/corrected_gene_lists/'
        ref_files = os.listdir(ref_dir)
        reference_eukaryotes = reference_eukaryotes + [i.split('.')[0] for i in ref_files]
    
    # Execute the 'worker' function for each model species
    for species in model_species:
        row = model_species_df.loc[model_species_df.species == species,:]
        worker(row)