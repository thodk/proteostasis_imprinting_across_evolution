# -*- coding: utf-8 -*-
"""
Created on Wed May 30 15:11:21 2018

@author: thodoris
"""

import pandas
import sys
import os
import multiprocessing
import numpy
sys.path.append('../')
from core_functions import construct_graph_from_mongo
from core_functions import uniprot_query
from core_functions import correct_uniprot_dataframe
from constant_variables import prokaryotic_classes
from constant_variables import define_main_dir


def worker(_id, input_list):

    query_columns = ['id','organism', 'genes(PREFERRED)', 
                     'genes(ALTERNATIVE)', 'go-id']

    df_columns = ['entry', 'organism', 'primary_gene_names', 
                  'synonym_gene_names', 'GO_terms']
                  
    annotation_data = {}
    
    for abbreviation, proteome_id in input_list:
        query = 'proteome:'+proteome_id       
        success, data, column_names = uniprot_query(query, query_columns)        
        if success == False:
            continue
        else:
            entries = []
            for s in data[1:]:
                fields = s.split('\t')
                if len(fields) < len(column_names):
                    continue
                entries.append(fields)        
        
        initial_df = pandas.DataFrame(data=entries, columns=column_names, index=None)
        initial_df.columns = df_columns
        
        f = correct_uniprot_dataframe
        final_df = f(initial_df, ['primary_gene_names', 'GO_terms'], 'primary_gene_names')
        final_df = final_df.loc[:, ['primary_gene_names', 'GO_terms']]
        
        tmp_dict = {}
        for i, row in final_df.iterrows():
            primary_gene = row.primary_gene_names
            terms = row.GO_terms.split('; ')
            for term in terms:            
                tmp_dict.setdefault(primary_gene, []).append(term)
            tmp_dict[primary_gene] = list(set(tmp_dict[primary_gene]))
        annotation_data.update({abbreviation:tmp_dict})

    G_GO_P = construct_graph_from_mongo('GO_P', mongo_database='BIM_background')
    G_GO_C = construct_graph_from_mongo('GO_C', mongo_database='BIM_background')
    G_GO_F = construct_graph_from_mongo('GO_F', mongo_database='BIM_background')  

    for i, (organism, organism_mapping) in enumerate(annotation_data.items()):
        tmp_data = []
        for gene, terms in organism_mapping.items():
            go_p_terms = list(set(G_GO_P.entries.keys()).intersection(terms))
            for term in go_p_terms:
                tmp_data.append([gene, term, 'Unknown', 'biological_process']) 
            go_c_terms = list(set(G_GO_C.entries.keys()).intersection(terms))
            for term in go_c_terms:
                tmp_data.append([gene, term, 'Unknown', 'cellular_component']) 
            go_f_terms = list(set(G_GO_F.entries.keys()).intersection(terms))
            for term in go_f_terms:
                tmp_data.append([gene, term, 'Unknown', 'molecular_function']) 
        df = pandas.DataFrame(tmp_data, columns=['external_gene_name', 'go_id', 
                                                 'go_linkage_type', 'namespace_1003'])
        df.to_csv(annotation_dir+organism+'_mapping_GO.tsv', sep='\t', index=False)

    queue.put({_id:'ok'})


'''
DESCRIPTION:

That script retrieves species' GO annotation from UniProt.

GO mapping: contains the GO terms - gene symbols (or ensembl ids) mapping
Output directory: /data/annotation_data//
'''

if __name__ == "__main__":
    
    for c in prokaryotic_classes: 

        main_dir = define_main_dir(c)
        annotation_dir = main_dir+'data/annotation_data/'
        if not os.path.exists(annotation_dir):
            os.makedirs(annotation_dir)    

        species_df = pandas.read_csv(main_dir+'/preprocessing/03_species_with_rRNA_and_HSPs_annotation.tsv', 
                                     sep='\t', index_col=0)

        pairs = zip(species_df.abbreviation, species_df.proteome_id)       
        queue = multiprocessing.Queue()
        procs = []
        nprocs = 16
        batch = int(numpy.math.floor(len(pairs)/float(nprocs)))
        for i in range(nprocs):
            if i == 15:
                tmp = pairs[i*batch:]
            else:
                tmp = pairs[i*batch:(i+1)*batch]
            p = multiprocessing.Process(target=worker, args=(str(i).zfill(2), tmp, ))
            p.start()
            procs.append(p)
        for p in procs:
            p.join()
