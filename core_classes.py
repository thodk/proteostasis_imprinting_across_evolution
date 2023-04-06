#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy
import pandas
import operator
import copy
import collections
import core_functions
import core_semantics
import string
symbols_for_random_ids = list(string.ascii_lowercase) + list(string.ascii_uppercase) + list(string.digits)


class Pool(object):


    def __init__(self):
        self.populations = {}
        self.population = 0
        self.counting_elements_population = 0

    
    def set_pool_content(self, mapping, pool_content="keys"):
        pool = []
        if pool_content == "values":
            for key, values in mapping.items():
                pool.extend(values)
            self.counting_elements_population = len(list(mapping.keys()))                 
        else:
            tmp = []
            for key, values_list in mapping.items():
                tmp.extend(values_list)
                pool.extend([key]*len(values_list))
            self.counting_elements_population = len(set(tmp))
        self.populations = collections.Counter(pool)
        self.population = len(pool)
    
    
    def get_population(self, entry=None):
        if entry==None:
            return self.population
        else:
            try:
                return self.populations[entry]
            except KeyError:
                return 0

    def get_unique_pool(self):
        return list(self.populations.keys())

    def get_sorted_pool(self):
        return sorted(self.populations.items(), key=operator.itemgetter(1), 
                      reverse=True)


class GraphNode(object):

    def __init__ (self, name, definition=None, description=None, category=None,
                  parents=[], parents_relations={}):
        self.name = name
        self.definition = definition
        self.description = description
        self.category = category
        self.parents = parents
        self.ancestors = parents
        self.children = []
        self.descendants = []
        self.ancestors_objects = {}
        for parent, relation in parents_relations.items():
            self.ancestors_objects.update({parent:AncestorGraphNode(parent, relation)})
        self.correct_relations()

    def set_relation(self, node, relation_type):
        if relation_type == "parent":
            self.parents.append(node)
            self.parents = list(set(self.parents))
        elif relation_type == "ancestor":
            self.ancestors.append(node)
            self.ancestors = list(set(self.ancestors))
        elif relation_type == "child":
            self.children.append(node)
            self.children = list(set(self.children))            
        elif relation_type == "descendant":
            self.descendants.append(node)
            self.descendants = list(set(self.descendants))            
        else:
            pass

    def correct_relations(self):
        self.parents = list(set(self.parents))    
        self.ancestors = list(set(self.ancestors))
        self.children = list(set(self.children))     
        self.descendants = list(set(self.descendants))         
        
    def get_relations(self):
        return {"parents":self.parents, "ancestors":self.ancestors,
                "children":self.children, "descendants":self.descendants}
        
    def is_obsolete(self):
        return False

    def is_alternative(self):
        return False


class AncestorGraphNode(object):

    def __init__ (self, entry, relation=None, paths=0, path_lengths=[]):
        self.name = entry
        self.relation = relation
        self.paths = paths
        self.path_lengths = path_lengths


class ObsoleteGraphNode(GraphNode):
  
    def __init__(self, entry, replace_entries=[], consider_entries=[]):
        super(ObsoleteGraphNode, self).__init__(entry)
        self.replace_entries = replace_entries
        self.consider_entries = consider_entries

    def is_obsolete(self):
        return True


class AlternativeGraphNode(GraphNode):
  
    def __init__(self, entry, major_entry):
        super(AlternativeGraphNode, self).__init__(entry)
        self.major_entry = major_entry

    def is_alternative(self):
        return True



class Collection(object):


    def __init__(self, name = None, entries = None):
        self.name = name
        self.entries = {}
        
    def set_entry(self):
        print("""Collection class is auxiliary. Please define your 
        Collection subclass, which should contain its specific setEntry
        definition.""")

    def has_entry(self, ID):
        try:
            self.entries[ID]
            return True
        except KeyError:
            return False

    def remove_entry(self, ID):
        del self.entries[ID]
        
    def get_entry_obj(self, ID):
        return self.entries[ID]

    def get_entries_dict(self):
        return self.entries

    def get_entries_list(self):
        return list(self.entries.keys())

    def __str__(self):
        return self.name




class Graph(Collection):


    def __init__(self, ontology="GO_P", organism="hsapiens", name="Graph"):        
        super(Graph, self).__init__(name)
        self.ontology = ontology
        self.organism = organism
        self.reference_mapping = {}     
        self.reference_pool = Pool()
        self.reference_mapping_status = "None"
        self.experiment_mapping = {}
        self.experiment_pool = Pool()
        self.experiment_mapping_status = "None"        

    def set_entry(self, entry, definition=None, description=None, category=None,
                  parents=[], parents_relations={}):
        obj = GraphNode(entry, definition, description, category, parents, parents_relations)
        self.entries.update({entry:obj})

    def set_obsolete_entry(self, entry, replace_entries=[], consider_entries=[]):
        obj = ObsoleteGraphNode(entry, replace_entries, consider_entries)
        self.entries.update({entry:obj})

    def set_alternative_entry(self, entry, major_entry):
        obj = AlternativeGraphNode(entry, major_entry)
        self.entries.update({entry:obj})

    def set_reference_pool(self, input_mapping, mapping_correction=True):
        if mapping_correction == True:
            self.reference_mapping = core_functions.correct_mapping(input_mapping, self)
            self.reference_mapping_status = "corrected"
        else:
            self.reference_mapping = copy.deepcopy(input_mapping) 
            self.reference_mapping_status = "direct"
        self.reference_pool.set_pool_content(self.reference_mapping)

    def set_experiment_pool(self, input_mapping, mapping_correction=True):        
        if mapping_correction == True: 
            self.experiment_mapping = core_functions.correct_mapping(input_mapping, self)
            self.experiment_mapping_status = "corrected"
        else:
            self.experiment_mapping = copy.deepcopy(input_mapping)
            self.experiment_mapping_status = "direct"
        self.experiment_pool.set_pool_content(self.experiment_mapping)

    def construct_connections_briefly(self):      
        for entry in list(self.entries.keys()):
            entry_obj = self.get_entry_obj(entry)
            if entry_obj.is_obsolete():
                continue
            parents = entry_obj.parents[:]
            if parents == []: # root node
                continue
            for parent in parents:
                entry_obj.set_relation(parent, "ancestor")
                parent_obj = self.get_entry_obj(parent)
                parent_obj.set_relation(entry, "child")
                parent_obj.set_relation(entry, "descendant")            
            while_loop = parents[:]
            while while_loop != []:
                for_loop = while_loop[:] # loop over a temporary list
                for ancestor in for_loop:
                    ancestor_obj = self.get_entry_obj(ancestor)
                    ancestor_parents = ancestor_obj.parents[:]
                    for new_ancestor in ancestor_parents:
                        entry_obj.set_relation(new_ancestor, "ancestor")
                        entry_obj.ancestors_objects.update({new_ancestor:AncestorGraphNode(new_ancestor)})
                        new_ancestor_obj = self.get_entry_obj(new_ancestor)
                        new_ancestor_obj.set_relation(entry, "descendant")
                        while_loop.append(new_ancestor)
                    while_loop.remove(ancestor)
                    while_loop = list(set(while_loop)) # remove redundancy

    def construct_connections_completely(self):
        for entry in list(self.entries.keys()):
            entry_obj = self.get_entry_obj(entry)
            if entry_obj.is_obsolete():
                continue
            parents = entry_obj.parents[:]
            if parents == []:
                continue
            for parent in parents:
                entry_obj.set_relation(parent, "ancestor")
                parent_obj = self.get_entry_obj(parent)
                parent_obj.set_relation(entry, "child")
                parent_obj.set_relation(entry, "descendant")
                ancestor_obj = entry_obj.ancestors_objects[parent]
                ancestor_obj.paths = 1
                ancestor_obj.path_lengths = [1]
    
            while_loop = zip(entry_obj.parents[:], [1]*len(entry_obj.parents[:]))
            while len(while_loop) > 0:
                for_loop = while_loop[:]
                for ancestor, path_length in for_loop:
                    new_ancestors = self.get_entry_obj(ancestor).parents
                    for new_ancestor in new_ancestors:
                        try:
                            new_ancestor_obj = entry_obj.ancestors_objects[new_ancestor]
                            new_ancestor_obj.paths = new_ancestor_obj.paths + 1
                            new_ancestor_obj.path_lengths.append(path_length+1)
                        except KeyError:
                            new_ancestor_obj = AncestorGraphNode(new_ancestor, paths=1,
                                                                 path_lengths=[path_length+1])
                        entry_obj.ancestors_objects.update({new_ancestor:new_ancestor_obj})
                        entry_obj.set_relation(new_ancestor, "ancestor")
                        new_ancestor_obj = self.get_entry_obj(new_ancestor)
                        new_ancestor_obj.set_relation(entry, "descendant")
                        while_loop.append((new_ancestor, path_length+1))
                    while_loop.remove((ancestor, path_length))


    def get_non_redundant_set(self, terms):
        output_terms = terms[:]
        for term in terms:
            output_terms = list(set(output_terms).difference(self.get_entry_obj(term).ancestors))
        return output_terms


    def get_ancsestry_pool(self, terms):
        output = []
        for term in terms:
            ancestors = self.get_entry_obj(term).ancestors
            output = output.extend([term]+ancestors)
        output = list(set(output))
        return output


    def substitution(self, new_term, obsolete_terms):
        
        new_obj = GraphNode(new_term, definition=" ".join(obsolete_terms))
        ref_mapping = []
        exp_mapping = []
        for obs_term in obsolete_terms:
            
            parents = self.get_entry_obj(obs_term).parents
            for parent in parents:
                obj = self.get_entry_obj(parent)
                obj.children.remove(obs_term)
                obj.set_relation(new_term, 'child')
                self.entries[parent] = obj
                new_obj.set_relation(parent, 'parent')

            ancestors = self.get_entry_obj(obs_term).ancestors
            for ancestor in ancestors:
                obj = self.get_entry_obj(ancestor)
                obj.descendants.remove(obs_term)
                obj.set_relation(new_term, 'descendant')
                self.entries[ancestor] = obj
                new_obj.set_relation(ancestor, 'ancestor')

            children = self.get_entry_obj(obs_term).children
            for child in children:
                obj = self.get_entry_obj(child)
                obj.parents.remove(obs_term)
                obj.set_relation(new_term, 'parent')
                self.entries[child] = obj
                new_obj.set_relation(child, 'child')            
            
            descendants = self.get_entry_obj(obs_term).descendants         
            for descendant in descendants:
                obj = self.get_entry_obj(descendant)
                obj.ancestors.remove(obs_term)
                obj.set_relation(new_term, 'ancestor')
                self.entries[descendant] = obj
                new_obj.set_relation(descendant, 'descendant')                
            
            if len(self.reference_mapping) != 0:
                try:
                    tmp = self.reference_mapping[obs_term]
                    ref_mapping = list(set(ref_mapping).union(tmp))
                    del self.reference_mapping[obs_term]
                except KeyError:
                    pass
            
            if len(self.experiment_mapping) != 0:
                try:
                    tmp = self.experiment_mapping[obs_term]
                    exp_mapping = list(set(exp_mapping).union(tmp))
                    del self.experiment_mapping[obs_term]
                except KeyError:
                    pass
        
            del self.entries[obs_term]
        
        if len(self.reference_mapping) != 0:
            for ancestor in new_obj.ancestors:
                try:                
                    tmp = self.reference_mapping[ancestor]
                    tmp = set(tmp).union(ref_mapping)
                    self.reference_mapping.update({ancestor:tmp})
                    self.reference_pool.populations[ancestor] = len(tmp)
                except KeyError:
                    pass
            self.reference_mapping.update({new_term:ref_mapping})            
            self.reference_pool.populations.update({new_term:len(ref_mapping)})
            self.reference_pool.population = sum(self.reference_pool.populations.values())            
            
        if len(self.experiment_mapping) != 0:
            for ancestor in new_obj.ancestors:
                try:                
                    tmp = self.experiment_mapping[ancestor]
                    tmp = set(tmp).union(ref_mapping)
                    self.experiment_mapping.update({ancestor:tmp})
                    self.experiment_pool.populations[ancestor] = len(tmp)
                except KeyError:
                    pass              
            self.experiment_mapping.update({new_term:exp_mapping})
            self.experiment_pool.populations.update({new_term:len(exp_mapping)})
            self.experiment_pool.population = sum(self.experiment_pool.populations.values())
            
        self.entries.update({new_term:new_obj})          


          


    

    

class Semantics(core_semantics.Semantics):
    
    def __init__(self, graph_instance, store_similarities=True):
        super(Semantics, self).__init__(graph_instance, store_similarities)






class GraphCluster(Collection):

    def __init__(self, name, keep_self=True, graph_node_instance=None, mapping=None):
        super(GraphCluster, self).__init__()
        self.name = name
        self.graph_node_instance = graph_node_instance
        self.mapping = mapping
        self.entries = []
        self.all_members = []
        self.keep_self = keep_self


    def set_entry(self, cluster):
        self.entries.append(cluster)

    
    def flatten_entries(self):
        
        def recursive_worker(graph_cluster_instance):
            if len(graph_cluster_instance.entries) == 0:
                return [graph_cluster_instance.name]
            else:
                return [recursive_worker(gc_instance) for gc_instance in graph_cluster_instance.entries]
        
        tmp = recursive_worker(self)
        all_members = []
        while len(tmp) > 0:
            for_loop = tmp[:]
            for entry in for_loop:
                if type(entry) == list:
                    tmp.extend(entry)
                else:
                    all_members.append(entry)
                tmp.remove(entry)
        if len(all_members) == 0 or self.keep_self==True:
            all_members.append(self.name)
        else:
            pass
        all_members = list(set(all_members))
        return all_members


    def get_intra_distance(self, graph_instance, metric, criterion):
        semantics = Semantics(graph_instance)
        all_members = self.flatten_entries()
        if len(all_members) < 2:
            return 0.
        else:
            if len(all_members) == 2:
                return 1-semantics.get_pairwise_similarity(all_members, metric, criterion)
            else:
                matrix = semantics.get_pairwise_similarity(all_members, metric, criterion)
                indices = numpy.triu_indices(matrix.shape[0], k=1)
                values = matrix[indices]
                s = numpy.sum(1-values)            
                return s


    def get_all_members(self):
        return self.flatten_entries()






class GraphClusteringProcessMICA(Collection):
    
    
    def __init__(self, graph_instance, entries, metric='resnik', criterion='graph_corpus', 
                 ancestors_set='mica', replace_method='mica'):
        super(GraphClusteringProcessMICA, self).__init__()
        for entry in entries:
            GC = GraphCluster(str(entry), keep_self=True)
            GC.set_entry(GraphCluster(str(entry), keep_self=True))
            self.entries.update({GC.name:GC})
        self.graph_instance = graph_instance
        self.criterion = criterion
        self.metric = metric
        self.ancestors_set = ancestors_set
        self.replace_method = replace_method
        self.semantics = Semantics(graph_instance)
        self.rounds = {}


    def set_entry(self, entry, keep_self=True):
        if type(entry).__name__ != 'GraphCluster':
            GC = GraphCluster(str(entry), keep_self)
            GC.set_entry(GraphCluster(str(entry), keep_self))
        else:
            GC = copy.deepcopy(entry)
        try:
            self.entries[GC.name]
        except KeyError:
            self.entries.update({GC.name:GC})


    def find_cluster3(self, cluster1, cluster2):
        mica = self.semantics.get_mica(cluster1, cluster2, self.criterion)
        return mica


    def construct_similarity_matrix(self):
        terms = sorted(list(self.entries.keys()))
        similarity_matrix = self.semantics.get_pairwise_similarity(terms, metric=self.metric, 
                                                                   criterion=self.criterion,
                                                                   ancestors_set=self.ancestors_set)
        numpy.fill_diagonal(similarity_matrix, 0)
        self.similarity_matrix_df = pandas.DataFrame(similarity_matrix, index=terms, columns=terms)


    def find_closest(self, threshold=0):
        max_value = self.similarity_matrix_df.values.max()
        if threshold > max_value:
            pairs = []
        else:
            indices = numpy.where(self.similarity_matrix_df.values == max_value)
            pairs = zip(indices[0], indices[1])
            pairs = list(tuple(sorted(pair)) for pair in pairs)
            pairs = list(set(pairs))
        return pairs



    def merge_one(self, cluster1, cluster2, cluster3):
        try:
            cluster3_obj = self.get_entry_obj(cluster3)
        except KeyError:
            cluster3_obj = GraphCluster(cluster3, keep_self=False)
            self.set_entry(cluster3_obj)

        cluster1_obj = self.get_entry_obj(cluster1)
        if cluster1 != cluster3:
            cluster3_obj.set_entry(cluster1_obj)
            del self.entries[cluster1]
        cluster2_obj = self.get_entry_obj(cluster2) 
        if cluster2 != cluster3:        
            cluster3_obj.set_entry(cluster2_obj)
            del self.entries[cluster2]
        self.entries[cluster3] = cluster3_obj

                
    def clustering_round(self, threshold=0):
        if len(self.rounds) == 0:
            self.construct_similarity_matrix()
        else:
            pass
        pairs = self.find_closest(threshold)
        if len(pairs) == 0:
            return False
        else:
            pass
        M = self.similarity_matrix_df
        pair = pairs[0]
        cluster1 = M.index[pair[0]]
        cluster2 = M.index[pair[1]]
        
        cluster3 = self.find_cluster3(cluster1, cluster2)
        indexes_to_remove = []
        if cluster1 != cluster3:
            indexes_to_remove.append(cluster1)
        if cluster2 != cluster3:
            indexes_to_remove.append(cluster2)
   
        M = M.drop(index=indexes_to_remove)
        M = M.drop(columns=indexes_to_remove)
        

        if len(indexes_to_remove) == 2:
            try:
                M.loc[cluster3]
            except KeyError:
                new_vector = self.semantics.get_pairwise_similarity([cluster3, M.index.tolist()], 
                                                                    metric=self.metric, 
                                                                    criterion=self.criterion,
                                                                    ancestors_set=self.ancestors_set)
                new_vector = new_vector.round(3) # shape (1,6)
                # append column of cluster3
                cluster3_df = pandas.DataFrame(new_vector.transpose(), columns=[cluster3], index=M.index)
                M = pandas.concat([M, cluster3_df], axis=1)
    
                # append row of cluster3
                new_vector = numpy.append(new_vector, 0) # shape (6,)
                cluster3_df = pandas.DataFrame(new_vector, columns=[cluster3], index=M.index.tolist()+[cluster3]).transpose()
                M = pandas.concat([M, cluster3_df], axis=0)

        self.rounds.update({len(self.rounds) + 1 : [cluster1, cluster2, cluster3]})
        self.merge_one(cluster1, cluster2, cluster3)
        self.similarity_matrix_df = M
        return True

            
    def clusteringK(self, k=1):
        while (len(self.entries) > k):
            if self.clustering_round():
                pass
            else:
                break
            

    def clusteringT(self, threshold=0):
        while (len(self.entries) > 1):
            if self.clustering_round(threshold=threshold):
                pass
            else:
                break






