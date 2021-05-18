# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 14:00:23 2018

@author: thodoris
"""

import sys
sys.path.append('../../../BioInfoMiner2/bim2/python_scripts/core/')
import core_classes
import core_functions
from constant_variables import mongodb
import pandas
import os
import time
import cPickle
import pymongo
import numpy
import subprocess
import copy
from bson.binary import Binary





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


def graph_dump_to_mongo(graph_instance, ontology, organism, annotation_type, mongo_collection):
    
    cPickle.dump(graph_instance, file("./graph_instances/"+organism+"_"+ontology+"_"+annotation_type+".obj", "wb"))
    size = os.path.getsize("./graph_instances/"+organism+"_"+ontology+"_"+annotation_type+".obj")
    if size < 15*10**6:
        objs = core_functions.split_Graph(graph_instance, 1)
    else:
        N = str(numpy.math.ceil( numpy.math.floor(size/float(10**6)) / 20.) * 20)[:-3]
        objs = core_functions.split_Graph(graph_instance, int(N)+1)
        
    for i, g in objs.items():
        print i,
        obj = cPickle.dumps(g)
        mongo_collection.insert([{"id":organism+"_"+ontology+"_"+str(i), "graph":Binary(obj)}]) 
        
    command = subprocess.Popen("rm ./graph_instances/"+organism+"_"+ontology+"_"+annotation_type+".obj", stdout=subprocess.PIPE, shell=True)
    command.communicate()









prokaryotes_class = sys.argv[1]
df = pandas.read_csv('./'+prokaryotes_class+'/preprocessing/STEP3_ranked_organisms_with_names.tsv', 
                     sep='\t', index_col=0)
organisms = df.abbreviation.tolist()



if not os.path.exists('./graph_instances'):
    os.makedirs('./graph_instances')

mongodb = 'BIM_background'
client = pymongo.MongoClient()
db = client[mongodb]
'''
collection = db["Graphs"]
collection.remove({})

collection = db["Graphs_corrected"]
collection.remove({})
'''
ontologies = ['GO_P']

for organism in organisms:
    for ontology in ontologies:
        tmp_ontology = ontology

        collection = db["_".join([organism, "mapping", ontology, "corrected"])]
        res = collection.count()
        if res > 0:
            continue

        # DIRECT ANNOTATION GRAPHS
        print organism, ontology
        start = time.time()
        G = core_functions.construct_graph_from_mongo(tmp_ontology, connections="briefly")
        print "initial amount of terms: ", len(G.entries)
        mapping = core_functions.get_mapping_from_mongo(ontology, organism, corrected_mapping=False)
        print "length of mapping :", len(mapping)

        G.set_reference_pool(mapping, mapping_correction=False)
        to_remove_list = find_unannotated(G)
        G = remove_unannotated(G, to_remove_list)
        G_entries = copy.deepcopy(G.entries)
        print "direct annotation: amount of terms: ", len(G.entries)      
        print "time: ", time.time() - start
        print "dumping and splitting..."
        collection = db['Graphs']
        #graph_dump_to_mongo(G, ontology, organism, "direct", collection)
        print "\n"
        
        
        # CORRECTED ANNOTATION GRAPHS        
        start = time.time()
        G_corrected = core_classes.Graph(ontology, organism)
        G_corrected.entries = copy.deepcopy(G_entries)
        G_corrected.set_reference_pool(mapping, mapping_correction=True)
        print "corrected annotation: amount of terms: ", len(G_corrected.entries) 
        print "time: ", time.time() - start        

        # INTERCALATED : CORRECTED ANNOTATION DUMPING  
        mapping = []
        for term, genes_list in G_corrected.reference_mapping.items():
            for gene in genes_list:
                entry = {"term_accession":term, "gene_symbol":gene}
                mapping.append(entry)
        print "corrected mapping: ", len(mapping)
        
        collection = db["_".join([organism, "mapping", ontology, "corrected"])]
        collection.remove({})
        collection.insert(mapping)
        #collection.ensure_index("gene_symbol")
        #collection.ensure_index("term_accession") 
        
        collection = db["_".join([organism, "mapping", ontology])]
        collection.remove({})
        collection.drop()
        # CORRECTED ANNOTATION GRAPHS  (CONTINUE)
        print "dumping and splitting..."
        collection = db['Graphs_corrected']
        #graph_dump_to_mongo(G_corrected, ontology, organism, "corrected", collection)
        print "\n"

client.close()
