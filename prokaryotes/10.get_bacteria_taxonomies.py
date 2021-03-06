# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 15:16:34 2020

@author: thodoris
"""

import pandas
import numpy
import os
import sys
import re
import multiprocessing
sys.path.append('../')
from core_functions import uniprot_query_with_limit
from constant_variables import define_main_dir



def worker(_id, ref_proteomes_list):

    query_columns = ['organism', 'lineage(PHYLUM)', 'lineage(CLASS)']
    df_columns = ['proteome_id', 'species', 'phylum', 'class']
                  
    if not os.path.isfile(main_dir+'taxonomy_details_'+_id+'.txt'):
        F_details = open(main_dir+'taxonomy_details_'+_id+'.txt', 'a')
        details = df_columns
        F_details.write('\t'.join(details)+'\n')
        F_details.close()
        proteomes = []
    else:
        df = pandas.read_csv(main_dir+'taxonomy_details_'+_id+'.txt', sep='\t')
        proteomes = df.proteome_id.tolist()
    
    for i, ref_proteome in enumerate(ref_proteomes_list):
        if ref_proteome in proteomes:
            continue
        query = 'proteome:'+ref_proteome
        success, data, column_names = uniprot_query_with_limit(query, query_columns, 5)
        if success == False:
            continue
        else:
            for entry in data[1:]:
                f = entry.split('\t')
                details = [ref_proteome, f[0], f[1], f[2]]
                F_details = open(main_dir+'taxonomy_details_'+_id+'.txt', 'a')
                F_details.write('\t'.join(details)+'\n')
                F_details.close()
                break

    queue.put({i:'ok'})



'''
DESCRIPTION:

That process is executed only for bacteria, in order to get their taxonomic
classification and select at least one species from all the Phyla. The script 
reads the 01_ranked_species_based_on_genomic_annotation.tsv file.

The whole process is parallelized, using the multiprocessing package (default
number of processes is 32). So 32 different output files are created in the
'scanning' directory. All these will be concatenated in the following steps
to clarify the species that will be used in the analysis. 



Output: the directory 'scanning_taxonomies' and bacteria_taxonomies.tsv
'''

if __name__ == "__main__":
    prokaryotic_class='bacteria'  
    
    main_dir = define_main_dir(prokaryotic_class)+'preprocessing/'
    if not os.path.exists(main_dir+'scanning_taxonomies/'):
        os.makedirs(main_dir+'scanning_taxonomies/')
    
    f = main_dir+'01_ranked_species_based_on_genomic_annotation.tsv'
    ref_proteomes_df = pandas.read_csv(f, sep='\t')
    ref_proteomes_list = ref_proteomes_df.proteome_id.tolist()
    
    queue = multiprocessing.Queue()
    procs = []
    nprocs = 32
    batch = int(numpy.math.floor(len(ref_proteomes_list)/float(nprocs)))
    for i in range(nprocs):
        if i == 31:
            tmp = ref_proteomes_list[i*batch:]
        else:
            tmp = ref_proteomes_list[i*batch:(i+1)*batch]
        p = multiprocessing.Process(target=worker, args=(str(i).zfill(2), tmp,))
        p.start()
        procs.append(p)
        
    for p in procs:
        p.join()


    files = [main_dir+'scanning_taxonomies/'+i for i in os.listdir(main_dir+'scanning_taxonomies/')]
    files = filter(lambda x: re.search('taxonomy_details', x), files)
    details_list = []
    for f in files:
        df = pandas.read_csv(f, sep='\t')
        details_list.append(df)

    details_df = pandas.concat(details_list, ignore_index=True)
    details_df = details_df.reset_index(drop=True)
    details_df.to_csv(main_dir+'bacteria_taxonomies.tsv', sep='\t')












