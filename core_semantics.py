#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy
from numpy import math
import pandas
import operator
import copy
import networkx as nx
from scipy.spatial import distance
from sklearn.cluster import AffinityPropagation as AP
import core_functions



class SemanticsBase(object):
    
    def __init__(self, graph_instance):
        self.graph_instance = graph_instance
        self.probabilities = {}
        self.information_contents = {}
        self.semantic_values = {}

        
    def get_lengths_of_elements(self, term, criterion):
        if criterion=="reference_mapping":
            numerator = int(self.graph_instance.reference_pool.get_population(term))
            denominator = float(self.graph_instance.reference_pool.counting_elements_population)
        elif criterion=="experiment_mapping":
            numerator = int(self.graph_instance.experiment_pool.get_population(term))
            denominator = float(self.graph_instance.experiment_pool.counting_elements_population)
        elif criterion=="graph_corpus":
            term_obj = self.graph_instance.get_entry_obj(term)
            numerator = len(term_obj.descendants) + 1
            denominator = float(len(self.graph_instance.entries))
        else:
            numerator = 0
            denominator = 1
        return numerator, denominator


    def get_elements(self, term, criterion):
        if criterion == "reference_mapping":
            try:
                elements = self.graph_instance.reference_mapping[term]
            except KeyError:
                 elements = []              
        elif criterion == "experiment_mapping":
            try:
                elements = self.graph_instance.experiment_mapping[term]
            except KeyError:
                 elements = []
        elif criterion == "graph_corpus":
            try:
                elements = self.graph_instance.get_entry_obj(term).ancestors
            except KeyError:
                 elements = []
        else:
            elements = []     
        return elements


    def get_probability(self, term, criterion):
        '''
        Relative Frequency of a term in the graph corpus or an annotation dataset

        Paper: A new measure for functional similarity of gene products based on
        Gene Ontology, Schlicker, 2006

        https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-7-302
        '''
        def worker(tmp_term):
            numerator, denominator = self.get_lengths_of_elements(tmp_term, criterion)
            prob = numerator/float(denominator)
            self.probabilities.setdefault(criterion, {}).update({tmp_term:prob})
            return prob
        
        if type(term) == list: # if term is a list of terms
            prob = 0
            for tmp_term in term:
                try:
                    prob += self.probabilities[criterion][tmp_term]
                except KeyError:
                    prob += worker(tmp_term)
            prob = prob / float(len(term))
        else: # single term
            try:
                prob = self.probabilities[criterion][term]
            except KeyError:
                prob = worker(term)   
        return prob


    def get_information_content(self, term, criterion):
        '''
        The fundamental idea of information content

        A Mathematical Theory of Communication, Shannon, 1948
        http://math.harvard.edu/~ctm/home/text/others/shannon/entropy/entropy.pdf
        '''
        def worker(tmp_term):
            numerator, denominator = self.get_lengths_of_elements(tmp_term, criterion)
            if numerator == 0 or numerator==denominator:
                ic = 0
            else:
                p = numerator/float(denominator)
                ic = -numpy.log2(p)
                max_value = -numpy.log2(1./denominator)
                ic = round(ic/max_value, 3) # normalization
            self.information_contents.setdefault(criterion, {}).update({tmp_term:ic})
            return ic

        if type(term) == list: # if term is a list of terms
            ic = 0
            for tmp_term in term:
                try:
                    ic += self.information_contents[criterion][tmp_term]
                except KeyError:
                    ic += worker(tmp_term)
            ic = round(ic / float(len(term)), 3)
        else: # single term
            try:
                ic = self.information_contents[criterion][term]
            except KeyError:
                ic = worker(term)
        return ic


    def get_semantic_weight(self, term, criterion):
        try:
            ic = self.information_contents[criterion][term]
        except KeyError: 
            ic = self.get_information_content(term, criterion)
        if ic==0:
            knowledge = 9999999.
        else:
            knowledge = 1/float(ic)
        norm_knowledge = round(1/(1 + numpy.exp(-knowledge)), 3)
        return norm_knowledge


    def get_semantic_value(self, term, criterion):
        try:
            return self.semantic_values[criterion][term]
        except KeyError:
            pass
        term_obj = self.graph_instance.get_entry_obj(term)
        ancestors = term_obj.ancestors+[term]
        sem_value = 0
        for ancestor in ancestors:
            sem_value += self.get_semantic_weight(ancestor, criterion)
        sem_value = round(sem_value, 3)
        self.semantic_values.setdefault(criterion, {}).update({term:sem_value})
        return sem_value


    def get_topological_position(self, term):
        def worker(obj):
            if len(obj.parents) == 0:
                return 1
            else:
                p = []
                for parent in obj.parents:
                    den = float(len(self.graph_instance.get_entry_obj(parent).children))
                    p.append(worker(self.graph_instance.get_entry_obj(parent))/den)
                return numpy.prod(p)
        return self.graph_instance.get_entry_obj(term)


    def get_topological_information(self, term):
        tp = self.get_topological_position(term)
        return round(-numpy.log2(tp),3)






class SemanticsAncestorsBase(SemanticsBase):
    
    def __init__(self, graph_instance):
        super(SemanticsAncestorsBase, self).__init__(graph_instance)


    def get_mica(self, term1, term2, criterion):
        ancestors1 = self.graph_instance.get_entry_obj(term1).ancestors
        ancestors2 = self.graph_instance.get_entry_obj(term2).ancestors
        common_ancestors = list(set(ancestors1).intersection(ancestors2))
        mica = None
        try:
            ancestors2.index(term1)
            mica = term1 # term1 is ancestor of term 2
        except (IndexError, ValueError):
            pass
        try:
            ancestors1.index(term2)
            mica = term2 # term 2 is ancestor of term 1
        except (IndexError, ValueError):
            pass
        if term1 == term2:
            mica = term1
        if mica == None:
            mica_ic = -1
            for ancestor in common_ancestors:
                ic = self.get_information_content(ancestor, criterion)
                if ic > mica_ic:
                    mica_ic = ic
                    mica = ancestor
        else:
            pass
        return mica


    def get_dca(self, term1, term2, criterion=None):
        ancestors1 = self.get_elements(term1, criterion="graph_corpus")
        ancestors2 = self.get_elements(term2, criterion="graph_corpus")
        common_ancestors = list(set(ancestors1).intersection(ancestors2))
        dca = None
        if len(ancestors1) == 0:
            dca = [term1]
        elif len(ancestors2) == 0:
            dca = [term2]
        else:
            pass
        try:
            ancestors2.index(term1)
            dca = [term1]
        except (IndexError, ValueError):
            pass
        try:
            ancestors1.index(term2)
            dca = [term2]
        except (IndexError, ValueError):
            pass
        if dca == None:
            for_loop = common_ancestors[:]
            dca = set(common_ancestors)
            for anc in for_loop:
                anc_obj = self.graph_instance.get_entry_obj(anc)
                dca = list(set(dca).difference(anc_obj.ancestors[:]))
        else:
            pass
        return dca


    def get_da(self, term, criterion=None):
        ancestors = self.get_elements(term, criterion="graph_corpus")
        da = None
        for_loop = ancestors[:]
        da = set(ancestors)
        for anc in for_loop:
            anc_obj = self.graph_instance.get_entry_obj(anc)
            da = list(set(da).difference(anc_obj.ancestors[:]))
        return da


    def get_dishin(self, term1, term2, criterion=None):
        '''
        dishin needs data for paths among terms, so the graph_instance
        should have een constructed with connections="complete" and not
        "briefly"
        '''
        ancestors1 = self.get_elements(term1, criterion="graph_corpus")
        ancestors2 = self.get_elements(term2, criterion="graph_corpus")
        common_ancestors = list(set(ancestors1).intersection(ancestors2))
        dishin = None
        if len(ancestors1) == 0:
            dishin = [term1]
        elif len(ancestors2) == 0:
            dishin = [term2]
        else:
            pass
        try:
            ancestors2.index(term1)
            dishin = [term1]
        except (IndexError, ValueError):
            pass
        try:
            ancestors1.index(term2)
            dishin = [term2]
        except (IndexError, ValueError):
            pass 
        if dishin == None:
            term1_obj = self.graph_instance.get_entry_obj(term1)
            term2_obj = self.graph_instance.get_entry_obj(term2)
            data = {}
            for ancestor in common_ancestors:
                paths1 = term1_obj.ancestors_objects[ancestor].paths
                paths2 = term2_obj.ancestors_objects[ancestor].paths
                pd = abs(paths1-paths2)
                ic = self.get_information_content(ancestor, 'graph_corpus')
                data.setdefault(pd, []).append([ancestor, ic])
            pds = sorted(list(data.keys()))
            dishin = []
            for pd in pds:
                pd_ancestors = sorted(data[pd], key=operator.itemgetter(1), reverse=True)
                pd_dishin = pd_ancestors[0][0]
                dishin.append(pd_dishin)
        else:
            pass
        return dishin


    def get_xgrasm(self, term1, term2, criterion):
        ancestors1 = self.get_elements(term1, criterion="graph_corpus")
        ancestors2 = self.get_elements(term2, criterion="graph_corpus")
        if len(ancestors1) == 0:
            xgrasm = [term1]
        elif len(ancestors2) == 0:
            xgrasm = [term2]
        else:
            xgrasm = list(set(ancestors1).intersection(ancestors2))
            ranking = list(map(lambda a: [a,self.get_information_content(a, criterion)], xgrasm))
            ranking = sorted(ranking, key=operator.itemgetter(1), reverse=True)
            xgrasm = [x[0] for x in ranking]
        return xgrasm

    
    def get_xgrasm_factor(self, ancestors_set, criterion):
        ranking = list(map(lambda a: self.get_information_content(a, criterion), ancestors_set))
        mica_ic = ranking[0]
        if mica_ic == 0:
            return 1.
        else:
            value = (1+(sum(ranking)-mica_ic)/mica_ic)/float(len(ancestors_set))
            return round(value, 3)
    








class SemanticsAncestorsOperator(SemanticsAncestorsBase):

    def __init__(self, graph_instance):
        super(SemanticsAncestorsOperator, self).__init__(graph_instance)
        self.ancestors_sets = {}
        self.ancestors_type_methods = {
            'mica': self.get_mica,
            'dca': self.get_dca,
            'dishin': self.get_dishin,
            'xgrasm': self.get_xgrasm
        }
        self.ancestors_type_factors = {
            'mica': lambda x,y: 1,
            'dca':   lambda x,y: 1,
            'dishin': lambda x,y: 1,
            'xgrasm': self.get_xgrasm_factor
        }
        self.ancestors_for_ic_calculation = {
            'mica': lambda x: x,
            'dca':   lambda x: x,
            'dishin': lambda x: x,
            'xgrasm': lambda x: x[0]              
        }

        
    def get_ancestors_set(self, term1, term2, ancestors_type, criterion='graph_corpus'):        
        key = (ancestors_type, criterion)
        pair = tuple(sorted([term1, term2]))
        try:
            return self.ancestors_sets[key][pair]
        except KeyError:
            f = self.ancestors_type_methods[ancestors_type]
            ancestors = f(term1, term2, criterion)
            self.ancestors_sets.setdefault(key, {}).update({pair:ancestors})
            return ancestors


    def get_ancestors_set_for_ic_calculation(self, ancestors_type, ancestors_set):
        f = self.ancestors_for_ic_calculation[ancestors_type]
        return f(ancestors_set)
        

    def get_ancestors_factor(self, ancestors_type, ancestors_set, criterion='graph_corpus'):
        f = self.ancestors_type_factors[ancestors_type]
        return f(ancestors_set, criterion)

    '''
    def get_ancestors_from_bounded_graph(self, term, lower_threshold=0.7, 
                                         upper_threshold=0.35, 
                                         metric='information_content',
                                         criterion='graph_corpus'):
        terms = [term]
        ancestors = []
        if metric == 'information_content':
            f = self.get_information_content
        else:
            f = self.get_semantic_value
        while len(terms) > 0:
            loop_terms = terms[:]
            for term in loop_terms:
                ic = f(term, criterion)
                if ic > lower_threshold:
                    parents_ics = numpy.array([f(t,'graph_corpus') for t in self.graph_instance.get_entry_obj(term).parents])
                    if any(parents_ics < upper_threshold):
                        ancestors.append(term)
                    else:
                        terms.extend(self.graph_instance.get_entry_obj(term).parents)
                elif ic <= lower_threshold and ic >= upper_threshold:
                    ancestors.append(term)
                else:
                    pass
                terms.remove(term)
        ancestors = list(set(ancestors))
        return ancestors
    '''

    def get_ancestors_from_bounded_graph(self, term, high_ic=1, low_ic=0,
                                         high_sem_value=100, low_sem_value=0,
                                         criterion='graph_corpus'):

        def is_below_lower_tree_thresholds(term):
            ic = self.get_information_content(term, criterion=criterion)
            sem_value = self.get_semantic_value(term, criterion=criterion)
            if ic < high_ic and sem_value < high_sem_value:
                return False
            else:
                return True

        def is_above_higher_tree_thresholds(term):
            ic = self.get_information_content(term, criterion=criterion)
            sem_value = self.get_semantic_value(term, criterion=criterion)            
            if ic > low_ic and sem_value > low_sem_value:
                return False
            else:
                return True

        terms = [term]
        ancestors = []        
        while len(terms) > 0:
            loop_terms = terms[:]
            for term in loop_terms:
                if is_below_lower_tree_thresholds(term):
                    parents = self.graph_instance.get_entry_obj(term).parents
                    if any([is_above_higher_tree_thresholds(p) for p in parents]):
                        ancestors.append(term)
                    else:
                        terms.extend(self.graph_instance.get_entry_obj(term).parents)
                elif is_above_higher_tree_thresholds(term):
                    ancestors.append(term)
                else:
                    ancestors.append(term)
                terms.remove(term)
        ancestors = list(set(ancestors))
        return ancestors        






class BinaryMetricsOperator(object):
    
    def __init__(self):
        '''
        a: intersection
        b: difference between set 1 and set 2
        c: difference between set 2 and set 1
        n: union
        '''
        self.coefficients = {
            'jaccard': lambda a,b,c: float(a)/(a+b+c),
            'dice': lambda a,b,c: 2.*a/(2.*a+b+c),
            'sokal_and_sneath': lambda a,b,c: float(a)/(a+2.*(b+c)),
            'gower_and_legendre': lambda a,b,c: a/(a+0.5*(b+c)),
            'mean_manhattan': lambda a,b,c: 1 - float(b+c)/(a+b+c),
            'vari': lambda a,b,c: 1 - (b+c)/(4.*(a+b+c)),
            'lance_and_williams': lambda a,b,c: 1 - (b+c)/(2.*a+b+c),
            'sized_difference': lambda a,b,c: 1 - float(math.pow(b+c,2))/math.pow(a+b+c,2),
            'cosine': lambda a,b,c: float(a)/math.sqrt((a+b)*(a+c)),
            'forbesi': lambda a,b,c: float(a)*(a+b+c)/((a+b)*(a+c)),
            'sorgenfrei': lambda a,b,c: float(math.pow(a,2))/((a+b)*(a+c)),
            'mountford': lambda a,b,c: a/(0.5*(a*b+a*c)+b*c),
            'mcconnaughey': lambda a,b,c: float(math.pow(a,2)-b*c)/((a+b)*(a+c)),
            'kulczynski': lambda a,b,c: 0.5*a*(2.*a+b+c)/((a+b)*(a+c)),
            'driver_and_kroeber': lambda a,b,c: 0.5*a*(1./(a+b) + 1./(a+c)),
            'johnson': lambda a,b,c: float(a)/(a+b) + float(a)/(a+c),
            'simpson': lambda a,b,c: float(a)/min(a+b, a+c),
            'braun_and_banquet': lambda a,b,c: float(a)/max(a+b, a+c),
        }


    def get_binary_similarity(self, metric, a, b, c):
        f = self.coefficients.get(metric)
        if a+b+c == 0:
            similarity = 0
        else:
            similarity = round(f(a,b,c), 3)
        return similarity








class AggregateMetricsOperator(object):
    
    def __init__(self):
        self.aggregate_metrics = {
            'average': self.get_average_score,
            'AVE': self.get_average_score,
            'maximum': self.get_max_score,
            'MAX': self.get_max_score,
            'best_match_average': self.get_best_match_average,
            'BMA': self.get_best_match_average,
            'best_match_maximum': self.get_best_match_maximum,
            'BMM': self.get_best_match_maximum,
            'average_best_matches': self.get_average_best_matches,
            'ABM': self.get_average_best_matches
        }
    
    
    def get_average_score(self, similarity_data):
        if type(similarity_data) == list:
            values1 = similarity_data[0].flatten()
            values2 = similarity_data[0].flatten()
            values = numpy.concatenate((values1, values2))
            similarity = round(values.mean(), 3)
        elif similarity_data.shape[0] == similarity_data.shape[1] and all(similarity_data.diagonal()) == 1:
            numpy.fill_diagonal(similarity_data, 0)
            condensed_matrix = distance.squareform(similarity_data, force='to_vector')
            similarity = round(condensed_matrix.mean(), 3)
        else:    
            similarity = round(similarity_data.mean(),3)
        return similarity


    def get_max_score(self, similarity_data):
        if type(similarity_data) == list:
            values1 = similarity_data[0].flatten()
            values2 = similarity_data[0].flatten()
            values = numpy.concatenate((values1, values2))
            similarity = round(values.max(), 3)
        elif similarity_data.shape[0] == similarity_data.shape[1] and all(similarity_data.diagonal()) == 1:
            numpy.fill_diagonal(similarity_data, 0)
            condensed_matrix = distance.squareform(similarity_data, force='to_vector')
            similarity = round(condensed_matrix.max(), 3)
        else:    
            similarity = round(similarity_data.max(),3)
        return similarity
    
    
    def get_best_match_average(self, similarity_data): # BMA
        if type(similarity_data) == list:
            matrix1 = similarity_data[0]
            mean_of_row_max1 = matrix1.max(axis=1).mean()
            matrix2 = similarity_data[1]
            mean_of_row_max2 = matrix2.max(axis=1).mean()
            similarity = (mean_of_row_max1+mean_of_row_max2)/2.
        elif similarity_data.shape[0] == similarity_data.shape[1] and all(similarity_data.diagonal()) == 1:
            numpy.fill_diagonal(similarity_data, 0)
            similarity = similarity_data.max(axis=1).mean()
        else:
            mean_of_row_max = similarity_data.max(axis=1).mean()
            mean_of_col_max = similarity_data.max(axis=0).mean()
            similarity = (mean_of_row_max+mean_of_col_max)/2.
        similarity = round(similarity,3)
        return similarity
        
    
    def get_best_match_maximum(self, similarity_data): # BMM
        if type(similarity_data) == list:
            matrix1 = similarity_data[0]
            mean_of_row_max1 = matrix1.max(axis=1).mean()
            matrix2 = similarity_data[1]
            mean_of_row_max2 = matrix2.max(axis=1).mean()
            similarity = max([mean_of_row_max1, mean_of_row_max2])
        elif similarity_data.shape[0] == similarity_data.shape[1] and all(similarity_data.diagonal()) == 1:
            numpy.fill_diagonal(similarity_data, 0)
            similarity = similarity_data.max()
        else:
            mean_of_row_max = similarity_data.max(axis=1).mean()
            mean_of_col_max = similarity_data.max(axis=0).mean()
            similarity = max([mean_of_row_max, mean_of_col_max])
        similarity = round(similarity,3)
        return similarity


    def get_average_best_matches(self, similarity_data): # ABM
        if type(similarity_data) == list:
            N = similarity_data[0].shape[0]          
            matrix1 = similarity_data[0]
            sum_of_row_max1 = matrix1.max(axis=1).sum()
            M = similarity_data[1].shape[0]  
            matrix2 = similarity_data[1]
            sum_of_row_max2 = matrix2.max(axis=1).sum()
            similarity = (sum_of_row_max1+sum_of_row_max2)/float(N+M)    
        elif similarity_data.shape[0] == similarity_data.shape[1] and all(similarity_data.diagonal()) == 1:
            numpy.fill_diagonal(similarity_data, 0)
            similarity = similarity_data.max(axis=1).mean()   
        else:
            N = similarity_data.shape[0]
            M = similarity_data.shape[1]
            sum_of_row_max = similarity_data.max(axis=1).sum()
            sum_of_col_max = similarity_data.max(axis=0).sum()
            similarity = (sum_of_row_max+sum_of_col_max)/float(N+M)
        similarity = round(similarity,3)
        return similarity


    def get_aggregate_similarity(self, similarity_matrix, metric):
        f = self.aggregate_metrics.get(metric)
        similarity = f(similarity_matrix)
        return similarity












class SemanticsPairwiseMetricsBase(SemanticsAncestorsOperator, BinaryMetricsOperator):

    
    def __init__(self, graph_instance):
        SemanticsAncestorsOperator.__init__(self, graph_instance)
        BinaryMetricsOperator.__init__(self)


    def get_resnik(self, term1, term2, criterion='graph_corpus', ancestors_type='mica'):
        ancestors = self.get_ancestors_set(term1, term2, ancestors_type)
        ancestors_for_ic = self.get_ancestors_set_for_ic_calculation(ancestors_type, ancestors)
        ancestors_factor = self.get_ancestors_factor(ancestors_type, ancestors, criterion)
        similarity = ancestors_factor * self.get_information_content(ancestors_for_ic, criterion)
        return similarity


    def get_lin(self, term1, term2, criterion='graph_corpus', ancestors_type='mica'):
        ancestors = self.get_ancestors_set(term1, term2, ancestors_type)
        ancestors_for_ic = self.get_ancestors_set_for_ic_calculation(ancestors_type, ancestors)
        ancestors_factor = self.get_ancestors_factor(ancestors_type, ancestors, criterion)
        anc_ic = self.get_information_content(ancestors_for_ic, criterion)
        ic1 = self.get_information_content(term1, criterion)
        ic2 = self.get_information_content(term2, criterion)       
        if (ic1+ic2) == 0:
            similarity = 0.
        else:
            similarity = ancestors_factor*2*anc_ic/(ic1+ic2)
        return round(similarity,3)

    
    def get_jiang_and_conrath(self, term1, term2, criterion='graph_corpus', ancestors_type='mica'):
        ancestors = self.get_ancestors_set(term1, term2, ancestors_type)
        ancestors_for_ic = self.get_ancestors_set_for_ic_calculation(ancestors_type, ancestors)
        ancestors_factor = self.get_ancestors_factor(ancestors_type, ancestors, criterion)
        anc_ic = self.get_information_content(ancestors_for_ic, criterion)
        ic1 = self.get_information_content(term1, criterion)
        ic2 = self.get_information_content(term2, criterion)
        distance = ic1+ic2 - 2.*anc_ic
        similarity = ancestors_factor*(1-distance/2.)
        return round(similarity,3)
    

    def get_pirro_and_euzenat(self, term1, term2, criterion='graph_corpus', ancestors_type='mica'):
        ancestors = self.get_ancestors_set(term1, term2, ancestors_type)
        ancestors_for_ic = self.get_ancestors_set_for_ic_calculation(ancestors_type, ancestors)
        ancestors_factor = self.get_ancestors_factor(ancestors_type, ancestors, criterion)
        anc_ic = self.get_information_content(ancestors_for_ic, criterion)
        ic1 = self.get_information_content(term1, criterion)
        ic2 = self.get_information_content(term2, criterion)
        if (ic1+ic2-anc_ic) == 0 or ic1 == 0 or ic2 == 0:
            similarity = 0.
        else:
            similarity = ancestors_factor*anc_ic/(ic1+ic2-anc_ic)
        return round(similarity,3)

    
    def get_schlicker(self, term1, term2, criterion='graph_corpus', ancestors_type='mica'):
        ancestors = self.get_ancestors_set(term1, term2, ancestors_type)
        anc_prob = self.get_probability(ancestors, criterion)
        sim_lin = self.get_lin(term1, term2, criterion, ancestors_type)
        similarity = sim_lin*(1-anc_prob)
        return round(similarity,3)


    def get_nunivers(self, term1, term2, criterion='graph_corpus', ancestors_type='mica'):
        ancestors = self.get_ancestors_set(term1, term2, ancestors_type)
        ancestors_for_ic = self.get_ancestors_set_for_ic_calculation(ancestors_type, ancestors)
        ancestors_factor = self.get_ancestors_factor(ancestors_type, ancestors, criterion)
        anc_ic = self.get_information_content(ancestors_for_ic, criterion)
        ic1 = self.get_information_content(term1, criterion)
        ic2 = self.get_information_content(term2, criterion)        
        similarity = ancestors_factor*anc_ic/max(ic1,ic2)
        return round(similarity,3)


    def get_simGIC(self, term1, term2, criterion=None, ancestors_type=None):
        # GIC or Jaccard with ICs
        ancestors1 = self.get_elements(term1, criterion="graph_corpus")
        ancestors2 = self.get_elements(term2, criterion="graph_corpus")
        common = list(set(ancestors1).intersection(ancestors2))
        diff1 = list(set(ancestors1).difference(ancestors2))
        diff2 = list(set(ancestors2).difference(ancestors1))        
        f = self.get_information_content
        a = sum([f(t, criterion) for t in common])        
        b = sum([f(t, criterion) for t in diff1])
        c = sum([f(t, criterion) for t in diff2])
        similarity = self.get_binary_similarity('jaccard', a,b,c)     
        return similarity


    def get_simDIC(self, term1, term2, criterion='graph_corpus', ancestors_type=None):
        # DIC or Dice coefficient with ICs
        ancestors1 = self.get_elements(term1, criterion="graph_corpus")
        ancestors2 = self.get_elements(term2, criterion="graph_corpus")
        common = list(set(ancestors1).intersection(ancestors2))
        diff1 = list(set(ancestors1).difference(ancestors2))
        diff2 = list(set(ancestors2).difference(ancestors1))        
        f = self.get_information_content
        a = sum([f(t, criterion) for t in common])        
        b = sum([f(t, criterion) for t in diff1])
        c = sum([f(t, criterion) for t in diff2])
        similarity = self.get_binary_similarity('dice', a, b, c)
        return similarity


    def get_simUIC(self, term1, term2, criterion='graph_corpus', ancestors_type='mica'):
        # UIC or Braun and Blanquet with ICs
        ancestors1 = self.get_elements(term1, criterion="graph_corpus")
        ancestors2 = self.get_elements(term2, criterion="graph_corpus")
        common = list(set(ancestors1).intersection(ancestors2))
        diff1 = list(set(ancestors1).difference(ancestors2))
        diff2 = list(set(ancestors2).difference(ancestors1))        
        f = self.get_information_content
        a = sum([f(t, criterion) for t in common])        
        b = sum([f(t, criterion) for t in diff1])
        c = sum([f(t, criterion) for t in diff2])
        similarity = self.get_binary_similarity('braun_and_banquet', a, b, c)     
        return similarity


    def get_simUI(self, term1, term2, criterion=None, ancestors_type=None):
        # Jaccard
        ancestors1 = self.get_elements(term1, criterion="graph_corpus")
        ancestors2 = self.get_elements(term2, criterion="graph_corpus")
        common = list(set(ancestors1).intersection(ancestors2))
        diff1 = list(set(ancestors1).difference(ancestors2))
        diff2 = list(set(ancestors2).difference(ancestors1))
        a = len(common)     
        b = len(diff1)
        c = len(diff2)
        similarity = self.get_binary_similarity('jaccard', a, b, c)
        return similarity


    def get_simDB(self, term1, term2, criterion=None, ancestors_type=None):
        # Dice
        ancestors1 = self.get_elements(term1, criterion="graph_corpus")
        ancestors2 = self.get_elements(term2, criterion="graph_corpus")
        common = list(set(ancestors1).intersection(ancestors2))
        diff1 = list(set(ancestors1).difference(ancestors2))
        diff2 = list(set(ancestors2).difference(ancestors1))
        a = len(common)     
        b = len(diff1)
        c = len(diff2)
        similarity = self.get_binary_similarity('dice', a, b, c)
        return similarity

        
    def get_simUB(self, term1, term2, criterion=None, ancestors_type=None):
        # Braun and Blanquet
        ancestors1 = self.get_elements(term1, criterion="graph_corpus")
        ancestors2 = self.get_elements(term2, criterion="graph_corpus")
        common = list(set(ancestors1).intersection(ancestors2))
        diff1 = list(set(ancestors1).difference(ancestors2))
        diff2 = list(set(ancestors2).difference(ancestors1))
        a = len(common)     
        b = len(diff1)
        c = len(diff2)
        similarity = self.get_binary_similarity('braun_and_banquet', a, b, c)
        return similarity      


    def get_simNTO(self, term1, term2, criterion=None, ancestors_type=None):
        # Simpson
        ancestors1 = self.get_elements(term1, criterion="graph_corpus")
        ancestors2 = self.get_elements(term2, criterion="graph_corpus")
        common = list(set(ancestors1).intersection(ancestors2))
        diff1 = list(set(ancestors1).difference(ancestors2))
        diff2 = list(set(ancestors2).difference(ancestors1))
        a = len(common)     
        b = len(diff1)
        c = len(diff2)
        similarity = self.get_binary_similarity('simpson', a, b, c)
        return similarity  


    def get_aggregate_ic(self, term1, term2, criterion='graph_corpus', ancestors_type=None):
        sv_term1 = self.get_semantic_value(term1, criterion)
        sv_term2 = self.get_semantic_value(term2, criterion)
        ancestors1 = self.get_elements(term1, "graph_corpus")
        ancestors2 = self.get_elements(term2, "graph_corpus")
        common_ancestors = list(set(ancestors1).intersection(ancestors2))
        sv_common = sum([self.get_semantic_weight(a, criterion) for a in common_ancestors])
        similarity = 2*sv_common / (sv_term1 + sv_term2)
        return round(similarity,3)


    def get_custom_aggregate_ic(self, term1, term2, criterion="graph_corpus", ancestors_type='mica'):
        sv_term1 = self.get_semantic_value(term1, criterion)
        sv_term2 = self.get_semantic_value(term2, criterion)
        ancestors1 = self.get_elements(term1, "graph_corpus")
        ancestors2 = self.get_elements(term2, "graph_corpus")
        common_ancestors = list(set(ancestors1).intersection(ancestors2))
        sv_common = sum([self.get_semantic_weight(a, criterion) for a in common_ancestors])
        resnik_sim = self.get_resnik(term1, term2, criterion, ancestors_type)
        similarity = (1+resnik_sim)*sv_common / (sv_term1 + sv_term2)
        return round(similarity,3)


    def get_reverse_eng_metric(self, term1, term2, criterion="graph_corpus", ancestors_type='mica'):
        resnik_sim = self.get_resnik(term1, term2, criterion, ancestors_type)
        
        genes1 = self.get_elements(term1, criterion='reference_mapping')
        genes2 = self.get_elements(term2, criterion='reference_mapping')
        a = len(list(set(genes1).intersection(genes2)))
        b = len(list(set(genes1).difference(genes2)))
        c = len(list(set(genes2).difference(genes1)))
        coeff = self.get_binary_similarity('dice', a, b, c)
        
        max_value = max([resnik_sim, coeff])
        min_value = min([resnik_sim, coeff])
        similarity = max_value + (1-max_value)*min_value
        return similarity


           
class SemanticsPairwiseMetricsOperator(SemanticsPairwiseMetricsBase):
    
    def __init__(self, graph_instance, store_similarities=True):
        super(SemanticsPairwiseMetricsOperator, self).__init__(graph_instance)
        self.similarities = {}
        self.store_similarities=store_similarities
        self.simple_metrics = {
            'resnik': self.get_resnik,
            'lin': self.get_lin,
            'jiang_and_conrath': self.get_jiang_and_conrath,
            'pirro_and_euzenat': self.get_pirro_and_euzenat,
            'schlicker': self.get_schlicker,
            'simGIC': self.get_simGIC,
            'simDIC': self.get_simDIC, 
            'simUIC': self.get_simUIC,            
            'nunivers': self.get_nunivers,
            'simUI': self.get_simUI,
            'simDB': self.get_simDB,
            'simUB': self.get_simUB,
            'simNTO': self.get_simNTO,
            'aggregate_ic': self.get_aggregate_ic,
            'custom_aggregate_ic': self.get_custom_aggregate_ic,
            'reverse_eng_metric': self.get_reverse_eng_metric
        }
        

    def get_pairwise_similarity(self, list_of_terms, metric, criterion='graph_corpus', 
                                ancestors_set='mica'):
        
        def worker(term1, term2):        
            if term1 == term2:
                return 1.
            #else:
            pair = tuple(sorted([term1, term2])) 
            key = (metric, '', ancestors_set)
            try:
                return self.similarities[key][pair]
            except KeyError:
                f = self.simple_metrics.get(metric)
                similarity = f(term1, term2, criterion, ancestors_set)
                if self.store_similarities == True:
                    self.similarities.setdefault(key, {}).update({pair:similarity})
                else:
                    pass
                return similarity

        output = None
        if len(list_of_terms) == 2:
            if type(list_of_terms[0]) == type(list_of_terms[1]) == str: # both are single terms
                similarity = worker(list_of_terms[0], list_of_terms[1])
                output = similarity
            else:
                if type(list_of_terms[0]) != list: # first element is single term
                    list_of_terms = [[list_of_terms[0]], list_of_terms[1]] # put the element in list
                if type(list_of_terms[1]) != list:
                    list_of_terms = [list_of_terms[0], [list_of_terms[1]]]
                N = len(list_of_terms[0])
                M = len(list_of_terms[1])
                matrix = numpy.ones((N,M))
                for n in range(N):
                    for m in range(M):
                        similarity = worker(list_of_terms[0][n], list_of_terms[1][m])
                        matrix[n,m] = similarity
                output = matrix
        elif len(list_of_terms) == 1:
            output = 1.
        else:
            N = len(list_of_terms)
            matrix = numpy.ones((N,N))
            for i in range(N):
                for j in range(i+1, N):
                    similarity = worker(list_of_terms[i], list_of_terms[j])
                    matrix[i,j] = matrix [j,i] = similarity
            output = matrix        
        return output
    






        
class SemanticsGroupwiseMetricsOperator(SemanticsPairwiseMetricsOperator, AggregateMetricsOperator):

    
    def __init__(self, graph_instance):
        SemanticsPairwiseMetricsOperator.__init__(self, graph_instance)
        AggregateMetricsOperator.__init__(self)


    def get_groupwise_intra_similarity(self, list_of_terms, pairwise_metric='resnik', 
                                       groupwise_metric='ABM', criterion='graph_corpus',
                                       ancestors_set='mica'):
                                       
        similarity_matrix = self.get_pairwise_similarity(list_of_terms, 
                                                         pairwise_metric,
                                                         criterion,
                                                         ancestors_set)
    
        similarity = self.get_aggregate_similarity(similarity_matrix, groupwise_metric)
        
        return similarity
    
    
    def get_groupwise_inter_similarity(self, pair_of_term_lists, pairwise_metric='resnik', 
                                       groupwise_metric='ABM', criterion='graph_corpus',
                                       ancestors_set='mica'):

        N = len(pair_of_term_lists[0])
        M = len(pair_of_term_lists[1])
        similarity_matrix = numpy.ones((N,M))
        for n in range(N):
            term1 = pair_of_term_lists[0][n]
            similarity_vector = self.get_pairwise_similarity([term1, pair_of_term_lists[1]],
                                                             pairwise_metric,
                                                             criterion, 
                                                             ancestors_set)
            similarity_matrix[n,] = similarity_vector

        similarity = self.get_aggregate_similarity(similarity_matrix, groupwise_metric)

        return similarity




class Semantics(SemanticsGroupwiseMetricsOperator):


    def __init__(self, graph_instance, store_similarities=True):
        super(Semantics, self).__init__(graph_instance)
        self.store_similarities=store_similarities


    def get_unique_terms(self, terms_list):
        initial_terms = terms_list[:]
        
        for term in initial_terms:
            terms_list = list(set(terms_list).difference(self.graph_instance.get_entry_obj(term).ancestors))    
          
        return terms_list       



    def get_term_lists_similarity(self, terms_dict, pairwise_metric='resnik', 
                                  groupwise_metric='ABM', criterion='graph_corpus',
                                  ancestors_set='mica', uniqueness=False,
                                  print_process=False):
        
        '''
        terms_dict: dictionary which contains list of terms and their ids.
                    For example: {'list1': [t1, t2], 'list2':[t1, t2, t3]}
        
        pairwise_metric: the basic semantic similarity metric that will be used 
                         in order to compare a pair of terms, default: resnik.
                         SemanticsPairwiseMetricsOperator class contains all the
                         available metrics in self.simple_metrics variable
        
        groupwise_metric: the metric that will be used to aggregate all the
                          pairwise similarity between two lists of terms, default: ABM.
                          ABM means Average Best Matches. AggregateMetricsOperator 
                          class contains all the available groupwise metrics.
                          
        criterion: graph_corpus means all the pairwise similarities will be 
                   calculated on the ontological graph (semantics). Otherwise
                   the reference_mapping could be used.
                   
        ancestors_set: a parameter which is used from some pairwise metric,
                      in order to define the common ancestors set between two
                      terms. Mica means that the metric will take into account
                      only the most informative common ancestor. A stricter 
                      method is xgrasm, which is based on the whole common 
                      ancestral tree.
        
        uniqueness: in order to compare two lists of terms, one could use those
                    lists, or their non-redundant sets. As non-redunant is defined
                    a list where only the most specific terms are included, as
                    their ancestors have been removed. 
        
        print_process: print in the console some elements about the procedure,
                       as it might take a long time to complete.
        '''
        
    
        def worker(terms_dict1, terms_dict2):

            list_of_terms1 = [terms_dict1['to_analyze'], terms_dict2['to_compare']]
            sim_matrix1 = self.get_pairwise_similarity(list_of_terms1,
                                                       metric=pairwise_metric, 
                                                       criterion='graph_corpus', 
                                                       ancestors_set='mica')
                                                       
            list_of_terms2 = [terms_dict2['to_analyze'], terms_dict1['to_compare']]                                         
            sim_matrix2 = self.get_pairwise_similarity(list_of_terms2,
                                                       metric=pairwise_metric, 
                                                       criterion='graph_corpus', 
                                                       ancestors_set='mica')
            
            similarity = self.get_aggregate_similarity([sim_matrix1, sim_matrix2],
                                                       groupwise_metric)
            return similarity


        terms_dict = copy.deepcopy(terms_dict)    
    
        if uniqueness == False:
            for key, terms_list in terms_dict.items():
                terms_dict.update({key: {'to_analyze':terms_list,'to_compare':terms_list}})
        else:
            for key, terms_list in terms_dict.items():
                to_analyze = self.get_unique_terms(terms_list)
                to_compare = terms_list
                terms_dict.update({key: {'to_analyze':to_analyze,'to_compare':to_compare}})

        keys = sorted(list(terms_dict.keys()))
        matrix = numpy.ones((len(keys), len(keys)))

        for i in range(len(keys)):
            key1 = keys[i]
            tmp_terms1 = terms_dict[key1]
            for j in range(i+1, len(keys)):
                if print_process == True:
                    pass
                else:
                    pass
                key2 = keys[j]
                tmp_terms2 = terms_dict[key2]
                similarity = worker(tmp_terms1, tmp_terms2)
                matrix[i,j] = matrix[j,i] = similarity                
        
        matrix_df = pandas.DataFrame(matrix, columns=keys, index=keys)
        return matrix_df



    def get_term_lists_distance(self, terms_dict, pairwise_metric='resnik', 
                                  groupwise_metric='ABM', criterion='graph_corpus',
                                  ancestors_set='mica', uniqueness=False,
                                  print_process=False):
        
        sim_matrix_df = self.get_term_lists_similarity(terms_dict, pairwise_metric,
                                                       groupwise_metric, criterion,
                                                       ancestors_set, uniqueness,
                                                       print_process)
        dis_matrix = 1. - sim_matrix_df

        return dis_matrix
        






class GenesSemanticsOperator(Semantics):
    
    
    def __init__(self, graph_instance):
        super(Semantics, self).__init__(graph_instance)
        self.genes_similarities = {}
        self.genes_mapping = core_functions.invert_mapping(graph_instance.reference_mapping)
        self.tmp_genes_list = []
        self.tmp_unannotated_genes = []
        self.tmp_genes_mapping = {}


    def upload_gene_list(self, genes_list, detect_specific_terms=True):
        self.tmp_genes_list = genes_list[:]
        self.tmp_unannotated_genes = []
        self.tmp_genes_mapping = {}
        for gene in self.tmp_genes_list:
            try:
                terms = self.genes_mapping[gene]
                if detect_specific_terms == True:
                    tmp_values = {
                        'specific_terms':self.graph_instance.get_non_redundant_set(terms),
                        'entire_corpus_terms':terms[:]
                    }
                else:
                    tmp_values = terms[:]
                self.tmp_genes_mapping.update({gene: tmp_values})
            except KeyError:
                self.tmp_unannotated_genes.append(gene)


    def get_genes_similarity(self, genes_dict, pairwise_metric, groupwise_metric='BMA', 
                             criterion='graph_corpus', ancestors_set='mica'):
    
        def worker(tmp_genes_dict):
            try:
                terms1_to_analyze = tmp_genes_dict.values()[0]['specific_terms']
                terms1_to_compare = tmp_genes_dict.values()[0]['entire_corpus_terms']
                terms2_to_analyze = tmp_genes_dict.values()[1]['specific_terms']
                terms2_to_compare = tmp_genes_dict.values()[1]['entire_corpus_terms']            
            except (TypeError, KeyError):
                terms1_to_analyze = tmp_genes_dict.values()[0]
                terms1_to_compare = tmp_genes_dict.values()[0]            
                terms2_to_analyze = tmp_genes_dict.values()[1]
                terms2_to_compare = tmp_genes_dict.values()[1]

            sim_matrix1 = self.get_pairwise_similarity([terms1_to_analyze, terms2_to_compare],
                                                       pairwise_metric, criterion='graph_corpus', 
                                                       ancestors_set='mica')
            sim_matrix2 = self.get_pairwise_similarity([terms2_to_analyze, terms1_to_compare],
                                                       pairwise_metric, criterion='graph_corpus', 
                                                       ancestors_set='mica')
            
            similarity = self.get_aggregate_similarity([sim_matrix1, sim_matrix2],
                                                       groupwise_metric)
            return similarity
        
        key = (pairwise_metric, criterion, ancestors_set, groupwise_metric)
        
        genes = sorted(list(genes_dict.keys()))
        matrix = numpy.ones((len(genes), len(genes)))

        for i in range(len(genes)):
            gene1 = genes[i]
            gene1_mapping = genes_dict[gene1]
            for j in range(i+1, len(genes)):
                gene2 = genes[j]
                gene2_mapping = genes_dict[gene2]

                pair = sorted([gene1, gene2])                
                try:
                    similarity = self.genes_similarities[key][pair]
                except KeyError:
                    similarity = worker({gene1:gene1_mapping, gene2:gene2_mapping})
                matrix[i,j] = matrix[j,i] = similarity                
        
        return matrix




    def construct_similarity_matrix(self, pairwise_metric, groupwise_metric='BMA', 
                                    criterion='graph_corpus', ancestors_set='mica'):
        
        similarity_matrix = self.get_genes_similarity(self.tmp_genes_mapping, 
                                                      pairwise_metric, 
                                                      groupwise_metric='BMA', 
                                                      criterion='graph_corpus', 
                                                      ancestors_set='mica')
        indices = sorted(list(self.tmp_genes_mapping.keys()))
        df = pandas.DataFrame(similarity_matrix, columns=indices, index=indices)
        return df


    def construct_semantic_similarity_network(self, similarity_matrix, threshold):
        values = similarity_matrix.values
        numpy.fill_diagonal(values, 0)
        indices = numpy.where(values >= threshold)
        network = nx.Graph()
        
        for i,gene in enumerate(similarity_matrix.index.tolist()):
            network.add_node(i, label=gene)
    
        for pair in zip(indices[0], indices[1]):
            network.add_edge(pair[0], pair[1])

        return network
        
        
    def export_hubs_from_network(self, network, small_network_threshold=5):        
        output = {}        
        for i, component in enumerate(nx.connected_components(network)):
            if len(component) == 1:
                gene = network.node[list(component)[0]]['label']
                output.setdefault(str(i), {}).update({'genes': [gene]})
                output.setdefault(str(i), {}).update({'hubs':[gene]})
                output.setdefault(str(i), {}).update({'network_nodes': list(component)})
            elif len(component) <= small_network_threshold:
                max_length = 0
                hub = None
                for node in component:
                    gene = network.node[node]['label']
                    try:
                        terms = self.tmp_genes_mapping[gene]['entire_corpus_terms']
                    except KeyError:
                        terms = self.tmp_genes_mapping[gene]
                    length = len(terms)
                    if length > max_length:
                        max_length = length
                        hub = gene
                    else:
                        pass
                    output.setdefault(str(i), {}).setdefault('genes', {}).update({gene:length})
                output.setdefault(str(i), {}).setdefault('hubs', []).append(hub)
                output.setdefault(str(i), {}).update({'network_nodes':list(component)})
            else:
                # 1. construct sub network
                to_remove_nodes = set(network.nodes).difference(component)
                sub_network = nx.restricted_view(network, to_remove_nodes, [])
                sub_network_nodes = sorted(sub_network.nodes())
                # 2 calculate shortest paths
                distances = nx.floyd_warshall_numpy(sub_network, nodelist=sub_network_nodes)
                similarities = 1 - distances / distances.max()
                # 3. calculate betweenness centrality
                betweenness_centrality = nx.betweenness_centrality(sub_network)
                # 4. apply affinity propagation
                preferences = [round(betweenness_centrality[j], 3) for j in sub_network_nodes]               
                ap_operator = AP(preference=numpy.array(preferences), affinity='precomputed')                
                labels = ap_operator.fit_predict(similarities)
                # 5. collect clusters
                clusters = sorted(list(set(labels)))
                for cluster in clusters:
                    indices = numpy.where(labels == cluster)
                    cluster_nodes = [sub_network_nodes[j] for j in list(indices[0])]
                    to_remove_nodes = set(network.nodes).difference(cluster_nodes)
                    cluster_network = nx.restricted_view(network, to_remove_nodes, [])
                    cluster_bc_dict = nx.degree_centrality(cluster_network)
                    cluster_bc_dict = dict( (network.node[j]['label'], v) for j,v in cluster_bc_dict.items())
                    threshold = numpy.percentile(cluster_bc_dict.values(), 80)
                    hubs = [j[0] for j in filter(lambda d: d[1] > threshold, cluster_bc_dict.items())]
                    exemplar_index = ap_operator.cluster_centers_indices_[cluster]
                    exemplar_node = sub_network_nodes[exemplar_index]
                    exemplar = cluster_network.node[exemplar_node]['label']                    
                    hubs = list(set(hubs).union([exemplar]))                    
                    k = str(i) + '_' + str(cluster)                    
                    output.setdefault(k, {}).update({'genes':cluster_bc_dict})
                    output.setdefault(k, {}).update({'hubs':hubs})
                    output.setdefault(k, {}).update({'network_nodes':list(cluster_network.nodes)})
        return output