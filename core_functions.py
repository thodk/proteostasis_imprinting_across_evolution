#!/usr/bin/python
# -*- coding: utf-8 -*-

import re
import sys
import copy
import numpy
import pymongo
import pandas
import requests
import core_classes


valid_keywords_for_hsp40 = ['heat shock 40', 'heat shock protein 40', 
                            'heat shock cognate 40', 'shock 40kDa', 
                            'shock 40 kDa', 'chaperone 40', 'hsp40', 'dnaj']
                            
valid_keywords_for_hsp70 = ['heat shock 70', 'heat shock protein 70', 
                            'heat shock cognate 70', 'heat shock cognate 71',
                            'shock 70kDa', 'shock 70 kDa', 'shock 71 kDa',
                            'chaperone 70', 'hsp70', 'dnak']

prohibited_keywords = ['interact', 'cooporate', 'stimulate', 'resemble', 
                       'contain', 'unknown function', 'with the hsp70', 
                       'with the hsp40', 'escort', 'fragment', 'pseudogene', 
                       'complex', 'sense overlapping', 'sense intronic', 
                       'antisense', 'binding', 'exchange factor',
                       'co-chaperone', 'suppressor', 'cooperates with',
                       'organizing protein', 'readthrough', 'divergent transcript']


def str_modification(s):
    special_chrs = list('-_!@#$%^&*(),.?":{}|<>]/ ')
    output = []
    for l in list(s):
        if l in special_chrs:
            output.append('\\'+l)
            continue
        else:
            pass
        try:
            int(l)
            output.append(l)
        except ValueError:
            output.append('['+l.lower()+l.upper()+']')
    output = ''.join(output)
    return output


def uniprot_query(query, columns, output_format='tab'):
    # 'proteome:proteome_id'

    columns_string = ','.join(columns)
    url = 'https://www.uniprot.org/uniprot/?query='+query+'&columns='+columns_string+'&format='+output_format

    request_times = 0
    success = False
    data = []
    
    while success==False and request_times<10:    
        try:
            response = requests.get(url=url)
            data = response.text.split('\n')
            success = True
            if len(data) == 0:
                success = False
                request_times += 1
            else:
                column_names = [str(i) for i in data[0].split('\t')]
                if len(column_names) != len(columns):
                    success = False
                    request_times+=1
        except requests.exceptions.ChunkedEncodingError:      
            request_times+=1

    return success, data, column_names


def uniprot_query_with_limit(query, columns, limit=5, output_format='tab'):
    
    columns_string = ','.join(columns)
    url = 'https://www.uniprot.org/uniprot/?query='+query+'&columns='+columns_string+'&limit='+str(limit)+'&format='+output_format

    request_times = 0
    success = False
    data = []
    
    while success==False and request_times<10:    
        try:
            response = requests.get(url=url)
            data = response.text.split('\n')
            success = True
            if len(data) == 0:
                success = False
                request_times += 1
            else:
                column_names = [str(i) for i in data[0].split('\t')]
                if len(column_names) != len(columns):
                    success = False
                    request_times+=1
        except requests.exceptions.ChunkedEncodingError:      
            request_times+=1

    return success, data, column_names


def correct_uniprot_dataframe(initial_df, columns_to_check, column_to_split):
    tmp_df = initial_df.loc[:,columns_to_check]
    tmp_df = tmp_df.replace('', numpy.nan)
    
    contains_nan = list(tmp_df.isnull().any(axis=1))
    not_contains_nan = map(lambda x: not x , contains_nan)
    df_without_nan = initial_df.iloc[not_contains_nan,:]

    pairs = zip(df_without_nan.index, df_without_nan.loc[:,column_to_split])
    multiple_indices_pairs = filter(lambda x: re.search(pattern='; ', string=x[1]), pairs)
    
    splitted_df = pandas.DataFrame(dict((i, {}) for i in initial_df.columns))
    for (index, entries) in multiple_indices_pairs:
        entries = entries.split('; ')
        for entry in entries:
            if entry=='':
                continue
            else:
                row = copy.deepcopy(df_without_nan.loc[index,:])
                row[column_to_split] = entry
                splitted_df = splitted_df.append(row)

    to_remove_indices = [i[1] for i in multiple_indices_pairs]
    valid_indices = list(set(df_without_nan.index).difference(to_remove_indices))
    
    df_without_multiple_genes = df_without_nan.loc[valid_indices]
    final_df = pandas.concat([df_without_multiple_genes, splitted_df], 
                             ignore_index=True, sort=False)
    final_df = final_df.reset_index(drop=True)

    return final_df


def filter_hsp_data(df, hsp):
    
    if hsp == 'hsp70':
        keywords = valid_keywords_for_hsp70
    else:
        keywords = valid_keywords_for_hsp70

    keywords_string = []
    for s in keywords:
        keywords_string.append(str_modification(s))
    keywords_string = '|'.join(keywords_string)
    
    prohibited_keywords_string = []
    for s in prohibited_keywords:
        prohibited_keywords_string.append(str_modification(s))                               
    prohibited_keywords_string = '|'.join(prohibited_keywords_string)
    
    pos_bool1 = df.protein_name.str.contains(keywords_string, regex=True)
    pos_bool2 = df.primary_gene_names.str.contains(keywords_string, regex=True)
    neg_bool = df.protein_name.str.contains(prohibited_keywords_string, regex=True)
    output_df = df.loc[((pos_bool1) | (pos_bool2)) & (~neg_bool), :]
    
    return output_df



def find_unannotated(G):
    to_remove = []
    for term in G.entries.keys():
        try:
            G.reference_mapping[term]
        except KeyError:
            mapping_len = 0
            terms_for_while = G.get_entry_obj(term).descendants[:]
            while (mapping_len == 0 and len(terms_for_while) > 0):
                terms_for_loop = terms_for_while[:]
                for term_desc in terms_for_loop:
                    try:
                        mapping = G.reference_mapping[term_desc]
                        mapping_len = mapping_len + len(mapping)
                        if mapping_len > 0:
                            break
                    except KeyError:
                        terms_for_while.remove(term_desc)
                        term_desc_obj = G.get_entry_obj(term_desc)
                        terms_for_while = list( set(terms_for_while) | set(term_desc_obj.descendants[:]) )
            if mapping_len == 0:
                to_remove.append(term)
    return to_remove



def remove_unannotated(G, to_remove_list):
    for term in G.entries.keys():
        term_obj = G.get_entry_obj(term)
        children = list(set(term_obj.children).difference(to_remove_list))
        descendants = list(set(term_obj.descendants).difference(to_remove_list))
        term_obj.children = children
        term_obj.descendants = descendants
        G.entries.update({term:term_obj})
    
    for term in to_remove_list:
        del G.entries[term]
        
    return G


def invert_mapping(mapping):
    output = {}
    for key, values in mapping.items():
        if type(values) == str:
            output.setdefault(values, []).append(key)
        else:
            for value in values:
                output.setdefault(value, []).append(key)
    return output


def get_mapping_from_mongo(ontology="GO_P", organism='hsapiens', keys="terms", 
                           corrected_mapping=True, specific_list=[], 
                           specific_list_content=None, mongo_database='Background'):

    client = pymongo.MongoClient()
    db = client[mongo_database]

    if corrected_mapping == True:
        collection = db[organism+'_mapping_'+ontology+'_corrected']
    else:
        collection = db[organism+'_mapping_'+ontology]


    if len(specific_list) == 0:
        if keys == 'genes':
            results = collection.aggregate([{'$group': {'_id': "$gene_symbol", 
                                             'list': {'$addToSet': '$term_accession'}}}
                                          ])
        else:
            results = collection.aggregate([{'$group': {'_id': "$term_accession", 
                                             'list': {'$addToSet': '$gene_symbol'}}}
                                          ])            
    else:  
        if keys == 'genes' and specific_list_content == 'genes':
            results = collection.aggregate([{'$match': {'gene_symbol': {'$in':specific_list}}}, 
                                            {'$group': {'_id': "$gene_symbol", 
                                             'list': {'$addToSet': '$term_accession'}}}
                                            ])
        elif keys == 'genes' and specific_list_content == 'terms':
            results = collection.aggregate([{'$match': {'term_accession': {'$in':specific_list}}}, 
                                            {'$group': {'_id': "$gene_symbol", 
                                             'list': {'$addToSet': '$term_accession'}}}
                                            ])
        elif keys == 'terms' and specific_list_content == 'genes':
            results = collection.aggregate([{'$match': {'gene_symbol': {'$in':specific_list}}}, 
                                            {'$group': {'_id': "$term_accession", 
                                             'list': {'$addToSet': '$gene_symbol'}}}
                                            ])
        elif keys == 'terms' and specific_list_content == 'terms':
            results = collection.aggregate([{'$match': {'term_accession': {'$in':specific_list}}}, 
                                            {'$group': {'_id': "$term_accession", 
                                             'list': {'$addToSet': '$gene_symbol'}}}
                                            ])
        else:
            results = []
    

    output_mapping = {}
    for entry in results:
        values = [i for i in entry['list'] if i != ""]
        output_mapping.update({entry['_id']: values})
    client.close()
    return output_mapping





def construct_graph_from_mongo(ontology, connections="briefly", 
                               mongo_database='Background'):
    G = core_classes.Graph() 
    
    client = pymongo.MongoClient()
    db = client[mongo_database]
    collection = db[ontology+"_relations"]

    query = collection.aggregate([ {'$group':{'_id':'$term_id', 
                                              'parents': {'$push':'$$ROOT.parents'},
                                              'children':{'$push':'$$ROOT.children'}, 
                                              'ancestors':{'$push':'$$ROOT.ancestors'},
                                              'descendants':{'$push':'$$ROOT.descendants'}}
                                   },
                                   {'$lookup':{'from':ontology+'_base', 
                                               'localField':'_id',
                                               'foreignField':'term_id',
                                               'as':'data'}
                                   }], allowDiskUse=True)
    
    for entry in query:
        data = entry['data'][0]
        term = str(data["term_id"])
        try:
            definition = data["definition"]
        except (KeyError):
            definition = None
        try:
            description = data["description"]
        except (KeyError):
            description = None
        try:
            category = str(data["category"])
        except KeyError:
            category = None 
    
        try:
            parents = {}
            term_parents = entry['parents'][0]
            for parent_dict in term_parents:
                parents.update({str(parent_dict['parent']):str(parent_dict['relation'])})
        except IndexError:
            parents = {}
        G.set_entry(term, definition=definition, description=description, 
                    category=category, parents=parents.keys(), parents_relations=parents)
        term_obj = G.get_entry_obj(term)
        
        
        try:
            term_ancestors = entry['ancestors'][0]
            for ancestor_dict in term_ancestors:
                ancestor = str(ancestor_dict["ancestor"])
                relation =  ancestor_dict["relation"]
                paths = 0 if connections == "briefly" else ancestor_dict["paths"]
                path_lengths = [] if connections == "briefly" else ancestor_dict["path_lengths"]
                ancestor_obj = core_classes.AncestorGraphNode(ancestor, relation=relation, 
                                                              paths=paths, path_lengths=path_lengths)
                term_obj.ancestors_objects.update({ancestor:ancestor_obj})
                term_obj.set_relation(ancestor, 'ancestor')            
        except IndexError:
            pass
        
        try:
            term_obj.children = entry['children'][0]
        except IndexError:
            pass
        try:
            term_obj.descendants = entry['descendants'][0]
        except IndexError:
            pass
        
    client.close()
    return G




