# -*- coding: utf-8 -*-

import pandas
import numpy
import os
import multiprocessing
import sys
sys.path.append('../')
from core_functions import construct_graph_from_mongo
from core_functions import invert_mapping
from core_functions import uniprot_query
from core_functions import correct_uniprot_dataframe
from core_functions import filter_hsp_data
from constant_variables import prokaryotic_classes
from constant_variables import define_main_dir


def worker(_id, ref_proteomes_list):

    # load GO graph to count the annotated terms for each proteome
    GOP_instance = construct_graph_from_mongo(ontology='GO_P', mongo_database='background')
    GOP_terms = GOP_instance.entries.keys()[:]

    query_columns = ['id', 'reviewed', 'organism', 'genes(PREFERRED)', 
                     'genes(ALTERNATIVE)', 'protein names', 'go-id', 
                     'annotation score']
    df_columns = ['entry', 'status', 'species', 'primary_gene_names', 
                  'synonym_gene_names', 'protein_name', 'GO_terms', 
                  'annotation_score']
    
    # check if that ref_proteomes_list has already been scanned    
    if not os.path.isfile(main_dir+'scanning/details_'+_id+'.txt'):
        F_details = open(main_dir+'scanning/details_'+_id+'.txt', 'a')
        details = ['proteome_id', 'organism', 'genes_amount', 'terms_amount', 
                   'hsp70_bool', 'hsp70', 'hsp40_bool', 'hsp40']
        F_details.write('\t'.join(details)+'\n')
        F_details.close()
        proteomes = []
    else:
        df = pandas.read_csv(main_dir+'scanning/details_'+_id+'.txt', sep='\t')
        proteomes = df.proteome_id.tolist()

    for i, ref_proteome in enumerate(ref_proteomes_list):
        if ref_proteome in proteomes:
            continue
        # if the ref_proteome has not been scanned, then send a query to UniProt
        query = 'proteome:'+ref_proteome
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
        final_df = correct_uniprot_dataframe(initial_df, ['primary_gene_names', 'GO_terms'],
                                             'primary_gene_names')
        species = list(set(final_df.organism))[0]
            
        #======================================================================
        # Examine species annotation
        #======================================================================
        species_gop_annotation = []
        gop_mapping = {}
        for gene_symbol, terms_string in zip(final_df.primary_gene_names, final_df.GO_terms):
            terms = terms_string.split('; ')
            gop_terms = list(set(GOP_terms).intersection(terms))
            for term in gop_terms:
                gop_mapping.setdefault(term, []).append(gene_symbol)
            species_gop_annotation.extend(gop_terms)
        species_gop_annotation = list(set(species_gop_annotation))
        
        extended_gop_annotation = species_gop_annotation[:]
        for term in species_gop_annotation:
            extended_gop_annotation.extend(GOP_instance.get_entry_obj(term).ancestors)
        extended_gop_annotation = list(set(extended_gop_annotation))        
        
        genes_mapping = invert_mapping(gop_mapping)

        #======================================================================
        # Examine HSPs
        #======================================================================
        hsp40_df = filter_hsp_data(final_df, 'hsp40')
        hsp40 = hsp40_df.entry.tolist()
        hsp40 = ', '.join(hsp40) if len(hsp40)>0 else ''
        hsp40_bool = False if hsp40 == '' else True
        
        hsp70_df = filter_hsp_data(final_df, 'hsp70')
        hsp70 = hsp70_df.entry.tolist()
        hsp70 = ', '.join(hsp70) if len(hsp70)>0 else ''
        hsp70_bool = False if hsp70 == '' else True

        #======================================================================
        # Write details
        #======================================================================
        details = [ref_proteome, species, str(len(genes_mapping)), 
                   str(len(extended_gop_annotation)), str(hsp70_bool), hsp70, 
                   str(hsp40_bool), hsp40]
        F_details = open(main_dir+'scanning/details_'+_id+'.txt', 'a')
        F_details.write('\t'.join(details)+'\n')
        F_details.close() 


    queue.put({_id:'ok'})



'''
DESCRIPTION:

That script reads the reference_proteomes.tab file, which should have been
downloaded manually from UniProt. For each reference proteome with CPD value
equal to 'Standard', a uniprot query is performed to some basic annotation 
('entry', 'status', 'species', 'primary_gene_names', 'synonym_gene_names',
'protein_name', 'GO_terms', 'annotation_score').

Then the amount of GO-BP annotated genes is count, using the custom version of
GO-BO database, loaded from MongoDB (database name = 'background'). Finally,
a search for annotated HSPs is perfomed and all these data are recorded in the
output data frame.

The whole process is parallelized, using the multiprocessing package (default
number of processes is 32). So 32 different output files are created in the
'scanning' directory. All these will be concatenated in the following steps
to clarify the species that will be used in the analysis. 



A detail file (output):
proteome_id	organism	genes_amount	terms_amount	hsp70_bool	hsp70	hsp40_bool	hsp40
UP000019772	Paenibacillus sabinae T27	486	1028	True	X4ZMJ1	True	X4ZG03
UP000284395	Altererythrobacter sp. HN-Y73	530	1084	True	A0A420ES37	True	A0A420ES75
UP000009877	Kocuria palustris PEL	370	888	True	M2YFL3	True	M2WH59
'''

if __name__ == "__main__":

    for c in prokaryotic_classes: 

        main_dir = define_main_dir(c)+'preprocessing/'
        if not os.path.exists(main_dir+'scanning/'):
            os.makedirs(main_dir+'scanning/')
    
        # reference_proteomes.tab contains some features to describe the
        # reference proteomes of Uniprot database
        ref_df = pandas.read_csv(main_dir+'reference_proteomes.tab',
                                           sep='\t', index_col=0)
        # keep only those with CPD == Standard
        # CPD: Complete Proteome Detector ---> proteome's protein count vs. the 
        # standard distribution of protein count expected for completeness
        ref_df = ref_df.loc[ref_df.CPD == 'Standard', :]
        
        # in order to accelerate the process, the multiprocessing package is used
        # to parallelize the worker function. By default number of processes is
        # defined to 32.
        queue = multiprocessing.Queue()
        procs = []
        nprocs = 32
        batch = int(numpy.math.floor(len(ref_df.index.tolist())/float(nprocs)))
        for i in range(nprocs):
            if i == 31:
                tmp = ref_df.index.tolist()[i*batch:]
            else:
                tmp = ref_df.index.tolist()[i*batch:(i+1)*batch]
            p = multiprocessing.Process(target=worker, args=(str(i).zfill(2), tmp,))
            p.start()
            procs.append(p)        
        for p in procs:
            p.join()















