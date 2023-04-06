#!/usr/bin/python
# -*- coding: utf-8 -*-
import pandas
import os
from constant_variables import define_main_dir
from constant_variables import ensembl_classes


def collect_data(main_dir, species_abbr, hsp):
    try:
        f = main_dir+'ensembl_data/'+hsp+'_raw_files/'+species_abbr+'.tsv'
        df = pandas.read_csv(f, sep='\t', index_col=0)
        ensembl_ids = df.ensembl_gene_id.tolist()
    except (IOError, AttributeError):
        ensembl_ids = []
    try:
        f = main_dir+'uniprot_data/'+hsp+'_raw_files/'+species_abbr+'.tsv'
        df = pandas.read_csv(f, sep='\t', index_col=0)
        uniprot_ids = df.id.tolist()
    except (IOError, AttributeError):
        uniprot_ids = []
    return {'ensembl':ensembl_ids, 'uniprot':uniprot_ids}
        

def worker(ensembl_class):
    
    main_dir = define_main_dir(ensembl_class)
    output_dir = main_dir+'filtered_HSPs_annotation/'
    files_dir = main_dir+'files/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Load the file of species names
    df = pandas.read_csv(files_dir+'final_species_names.tsv', sep='\t', index_col=0)
    
    # 1. For each abbreviation search if there is available HSPs annotation
    # 2. Usage of the 'collect_data' function to search in the 'raw files' folders 
    annotation_results = []
    data = {}
    for species_abbr, species in zip(df.abbreviation.tolist(), df.species_name.tolist()):
        hsp70_annotation = collect_data(main_dir, species_abbr, 'hsp70')
        hsp40_annotation = collect_data(main_dir, species_abbr, 'hsp40')
        data.update({species_abbr : {'hsp70':hsp70_annotation, 
                                     'hsp40':hsp40_annotation}})
        
        # True if count of annotated hsp70 > 0
        s70 = False if sum([len(i) for i in hsp70_annotation.values()]) == 0 else True
        # True if count of annotated hsp40 > 0
        s40 = False if sum([len(i) for i in hsp40_annotation.values()]) == 0 else True
        annotation_results.append([species_abbr, species, s70, s40])

    columns = ['abbreviation', 'species', 'hsp70_bool', 'hsp40_bool']
    df = pandas.DataFrame(annotation_results, columns=columns)

    pos_bool1 = df.hsp70_bool == True 
    pos_bool2 = df.hsp40_bool == True
    filtered_df = df.loc[((pos_bool1) & (pos_bool2)), :]
    filtered_df.to_csv(main_dir+'01_species_with_HSPs_annotation.tsv', sep='\t')
    
    # Save the HSPs annotation for each species in the 'filtered_hsp_annotation'
    # subfolder. Each data frame should contain three columns: protein id, 
    # hsp family name and database name (Ensembl or UniProt)
    for species_abbr in filtered_df.abbreviation.tolist():
        tmp_annotation = data[species_abbr]
        tmp_list = []
        for hsp, tmp_dict in tmp_annotation.items():
            # key: ensembl or uniprot
            # hsp: hsp40 or hsp70
            # ids: protein ids
            for key, ids in tmp_dict.items():
                tmp_list = tmp_list + [[i,hsp,key] for i in ids]
        tmp_df = pandas.DataFrame(tmp_list, columns=['gene_symbol', 'hsp', 'resource'])
        tmp_df.to_csv(output_dir+species_abbr+'.tsv', sep='\t')
   



if __name__ == "__main__":

    '''
    DESCRIPTION:
    
    This script reads the HSP annotations in the 'raw_files' subfolders, which 
    have been retrieved with the scripts 10A and 10B. The annotated entries
    are gathered together to examine which species have adequate HSP annotation 
    and which should be removed.
    
    Two outputs:
    
    1. '01_organisms_with_HSPs_annotation.tsv': This file contains only the 
    species with adequate annotation. They have passed the first filtering. 
    In order to included them in the final set for the analysis, they should 
    also pass the filtering for rRNA data and genomic GO annotation (in the 
    following steps).
    
    2. 'filtered_HSPs_annotation': This subfolder in 'preprocessing' contains 
    an HSP annotation file for each species:
    
    example:
    
    gene_symbol	hsp	     resource
    Q03751          hsp40	     uniprot
    Q9VPQ2	     hsp40	     uniprot
    P82910	     hsp70	     uniprot
    Q9VG58	     hsp70	     uniprot
    
    The gene ids will be used to retrieve the respective fasta files in the 
    following steps.
    '''

    # Execute the 'worker' function for each ensembl class
    [worker(c) for c in ensembl_classes] 
