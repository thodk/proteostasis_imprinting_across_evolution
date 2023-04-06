#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import pymongo
import subprocess
from core_functions import mongodb
import core_classes


def read_obo(obo_file):

    G = core_classes.Graph()
    have_term = False

    F = open(obo_file, "r")
    for line in F:
        if line.startswith("[Term]"):
            entry = None
            definition = None
            description = None
            category = None
            obsolete = False
            parents = []
            relations = {}
            alternative_entries = []
            replace_entries = []
            consider_entries = []
            have_term = True
        if line in ['\n', '\r\n']: # if the line is empty
            if have_term: # check if a term has been already defined
                if (type(description)==str and "OBSOLETE" in description): # obsolete condition
                    G.set_obsolete_entry(entry, replace_entries=replace_entries, 
                                         consider_entries=consider_entries)
                elif obsolete==True: # obsolete condition
                    G.set_obsolete_entry(entry, replace_entries=replace_entries, 
                                         consider_entries=consider_entries) 
                elif (os.path.basename(obo_file)[0:2] in ["GO"] and category==None): # obsolete condition
                    G.set_obsolete_entry(entry, replace_entries=replace_entries, 
                                         consider_entries=consider_entries)
                else:
                    G.set_entry(entry, definition=definition, description=description, 
                                category=category, parents=parents, parents_relations=relations)
                    if len(alternative_entries) != 0: # alternative condition
                        for alternative_entry in alternative_entries:
                            G.set_alternative_entry(alternative_entry, entry)
                have_term = False # current term has been recorded, so start to search for another one
            else:
                pass
        if have_term:
            line = line.replace("\n", "")
            if line.startswith("id:"):
                entry = line.split(" ")[1]
            elif line.startswith("name:"):
                definition = line.split(" ", 1)[1]
            elif line.startswith("def:"):
                description = line.split(" ", 1)[1]
                fields = description.split("\" [")
                if len(fields) > 1:
                    description = fields[0][1:]
            elif line.startswith("namespace:"):
                category = line.split(" ")[1]
                if category == "biological_process":
                    category = "P"
                elif category == "molecular_function":
                    category = "F"
                elif category == "cellular_component":
                    category = "C"                    
            elif (line.startswith("is_a") and not line.startswith("is_anon")):
                parent = line.split(" ")[1]
                if parent.split(':')[0] != entry.split(':')[0]:
                    pass
                else:
                    parents.append(parent)
                    relations.update({parent:"is_a"})
            elif line.startswith("relationship"):
                parent = line.split(" ")[2]
                relationship = line.split(" ", 2)[1]
                if relationship == "part_of":
                    if parent.startswith("NCBI"):
                        pass
                    else:
                        parents.append(parent)
                        relations.update({parent:"is_a"})
                else:
                    pass
            elif line.startswith("alt_id"):
                alternative_entries.append(line.split(" ")[1])
            elif line.startswith("is_obsolete: true"):
                obsolete = True
            elif line.startswith("replaced_by"):
                replace_entries.append(line.split(" ")[1])
            elif line.startswith("consider"):
                consider_entries.append(line.split(" ")[1])   
            else:
                pass
    F.close()
    return G


def reform_graph_and_save_to_mongo(graph_instance, ontology, mongo_database):
    # 1. G : Graph instance, where the construct_connections_extendedly function
    # has beed already executed
    # 2. database: GO, HPO, Reactome ...
    base_dict = {}
    relations_dict = {}

    for entry, entry_obj in graph_instance.entries.items():    
        # 1. obsolete, alternative, valid terms
        if entry_obj.is_obsolete():
            a = {"term_id":entry, "consider":entry_obj.consider_entries, 
                "replaced_by":entry_obj.replace_entries, "obsolete":True}
            base_dict.setdefault(ontology+"_obsoletes", []).append(a)
        elif entry_obj.is_alternative():
            a = {"term_id":entry, "major_id":entry_obj.major_entry, "alternative":True}
            base_dict.setdefault(ontology+"_alternatives", []).append(a)            
        else:
            a = {"term_id":entry, "definition":entry_obj.definition, 
                 "description":entry_obj.description, "category": entry_obj.category}
            if entry_obj.category == None:
                base_dict.setdefault(ontology, []).append(a)
            else:
                base_dict.setdefault(ontology+"_"+entry_obj.category, []).append(a)
        # 2. set relations
        parents = entry_obj.parents
        ancestors = entry_obj.ancestors
        children = entry_obj.children
        descendants = entry_obj.descendants
        if (len(parents) > 0):
            tmp_list = []
            for parent in parents:
                parent_obj = entry_obj.ancestors_objects[parent]
                tmp_list.append({"parent":parent, "relation":parent_obj.relation})
            a = {"term_id":entry, "parents":tmp_list}
            if entry_obj.category == None:
                relations_dict.setdefault(ontology, []).append(a)
            else:    
                relations_dict.setdefault(ontology+"_"+entry_obj.category, []).append(a)
            # import ancestors
            tmp_list = []
            for ancestor in ancestors:
                ancestor_obj = entry_obj.ancestors_objects[ancestor]   
                tmp_list.append({"ancestor":ancestor, "paths":ancestor_obj.paths,
                                 "path_lengths":ancestor_obj.path_lengths,
                                 "relation":ancestor_obj.relation})
                                 
            a = {"term_id":entry, "ancestors":tmp_list}
            if entry_obj.category == None:
                relations_dict.setdefault(ontology, []).append(a)
            else:    
                relations_dict.setdefault(ontology+"_"+entry_obj.category, []).append(a)
            
        if (len(children) > 0):
            a = {"term_id":entry, "children":children}
            if entry_obj.category == None:
                relations_dict.setdefault(ontology, []).append(a)
            else:    
                relations_dict.setdefault(ontology+"_"+entry_obj.category, []).append(a)
            a = {"term_id":entry, "descendants":descendants}
            if entry_obj.category == None:
                relations_dict.setdefault(ontology, []).append(a)
            else:    
                relations_dict.setdefault(ontology+"_"+entry_obj.category, []).append(a)

    client = pymongo.MongoClient()
    db = client[mongo_database]
    for k, base in base_dict.items():
        collection1 = db[k+"_base"]
        collection1.delete_many({})
        collection1.insert_many(base)
        collection1.create_index("term_id")            
    for k, relation in relations_dict.items():
        collection2 = db[k+"_relations"]
        collection2.delete_many({}) 
        collection2.insert_many(relation)
        collection2.create_index("term_id")  
    client.close()


if __name__ == '__main__':

    wget_files_dir = './files/'
    if not os.path.exists(wget_files_dir):
        os.makedirs(wget_files_dir) 

    url = "http://geneontology.org/ontology/go-basic.obo"
    output_file = wget_files_dir+"GO.obo" 
    p = subprocess.Popen("wget -O ./"+output_file+" "+url+"", stdout=subprocess.PIPE, shell=True)
    p.communicate()
    
    graph_instance = read_obo(output_file)
    graph_instance.construct_connections_briefly()
    reform_graph_and_save_to_mongo(graph_instance, 'GO', mongodb)    

