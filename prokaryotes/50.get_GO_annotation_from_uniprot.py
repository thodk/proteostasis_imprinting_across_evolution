#!/usr/bin/python
# -*- coding: utf-8 -*-
import pandas
import sys
import os
import multiprocessing
import numpy
sys.path.append('../')
from core_functions import mongodb
from core_functions import construct_graph_from_mongo
from core_functions import uniprot_query
from core_functions import correct_uniprot_dataframe
from constant_variables import prokaryotic_classes
from constant_variables import define_main_dir


def worker(_id, input_list, GO_P_terms):

    # Query attributes
    query_columns = ['id','organism', 'genes(PREFERRED)', 'genes(ALTERNATIVE)', 
                     'go-id']
    # AAlternative names for the query attributes
    df_columns = ['entry', 'organism', 'primary_gene_names', 'synonym_gene_names', 
                  'GO_terms']

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



    # Cross-check the valid GO terms for each species
    for i, (species, species_mapping) in enumerate(annotation_data.items()):
        tmp_data = []
        for gene, terms in species_mapping.items():
            go_p_terms = list(set(GO_P_terms).intersection(terms))
            for term in go_p_terms:
                tmp_data.append([gene, term, 'Unknown', 'biological_process']) 
            #go_c_terms = list(set(G_GO_C.entries.keys()).intersection(terms))
            #for term in go_c_terms:
            #    tmp_data.append([gene, term, 'Unknown', 'cellular_component']) 
            #go_f_terms = list(set(G_GO_F.entries.keys()).intersection(terms))
            #for term in go_f_terms:
            #    tmp_data.append([gene, term, 'Unknown', 'molecular_function']) 
        df = pandas.DataFrame(tmp_data, columns=['external_gene_name', 'go_id', 
                                                 'go_linkage_type', 'namespace_1003'])
        df.to_csv(annotation_dir+species+'_mapping_GO.tsv', sep='\t', index=False)

    queue.put({_id:'ok'})




if __name__ == "__main__":

    '''
    DESCRIPTION:
    
    This script download species' GO annotation from the UniProt database.
    
    GO mapping: contains the GO terms - gene symbols (or ensembl ids) mapping
    
    Output directory: '/data/class/data/annotation_data/'
    '''
    
    # Load GO graphs from MongoDB
    G_GO_P = construct_graph_from_mongo('GO_P', mongo_database=mongodb)
    GO_P_terms = list(G_GO_P.entries.keys())[:]
    #G_GO_C = construct_graph_from_mongo('GO_C', mongo_database=mongodb)
    #G_GO_F = construct_graph_from_mongo('GO_F', mongo_database=mongodb)  
    
    for c in prokaryotic_classes: 

        main_dir = define_main_dir(c)
        annotation_dir = './tsv/'
        if not os.path.exists(annotation_dir):
            os.makedirs(annotation_dir)    

        species_df = pandas.read_csv(main_dir+'/preprocessing/03_species_with_rRNA_and_HSPs_annotation.tsv', 
                                     sep='\t', index_col=0)

        pairs = list(zip(species_df.abbreviation, species_df.proteome_id))       
        queue = multiprocessing.Queue()
        procs = []
        nprocs = 32
        batch = int(numpy.math.floor(len(pairs)/float(nprocs)))
        for i in range(nprocs):
            if i == 31:
                tmp = pairs[i*batch:]
            else:
                tmp = pairs[i*batch:(i+1)*batch]
            p = multiprocessing.Process(target=worker, args=(str(i).zfill(2), tmp, GO_P_terms,))
            p.start()
            procs.append(p)
        for p in procs:
            p.join()
