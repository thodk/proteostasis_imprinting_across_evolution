#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import pandas
import os
import operator
import decimal
import pymongo
sys.path.append('../')
from core_classes import EnrichmentAnalysis
import core_functions


def core_process(gene_ids, ontology, organism, corrected_mapping, n_boot, 
                 hyper_pvalue_threshold, corrected_pvalue_threshold, mongo_database):


    P = EnrichmentAnalysis()
    
    # 1. background mapping & pool creation
    P.set_background_pool(ontology, organism, corrected_mapping=corrected_mapping, 
                          mongo_database=mongo_database)

    # 2. sample mapping & pool creation
    P.set_sample_pool(gene_ids)

    # 3. terms mapping creation
    P.set_terms_mapping()

    # 4. generate the enrichments from the above mappings & pools
    P.set_enrichments()
    P.set_elements_pool() # generate the elements - a list with enrichments
    unique_elements = P.get_elements()  # unique elements
    elements_counts = P.get_elements_populations()

    # hypergeometric test on elements
    Pt = P.get_background_population()
    St = P.get_sample_population()
    hyper_pvalues = core_functions.hypergeometric_test_for_enrichments(Pt, St, unique_elements)
    
    
    # rank the elements based on their count and their hyperg pvalue
    data_for_ranking = [[es, elements_counts[es], hyper_pvalues[es]] for es in hyper_pvalues.keys()]
    ranked_data = sorted(data_for_ranking, key=operator.itemgetter(1,2))
    # rankedData: [ [es1,count,hp], [es2,count,hp], [es3,count,hp] ]

    if len(ranked_data) == 0:
        return []

    # bootstrapping test on indexes
    corrected_pvalues = core_functions.bootstrapping_test_for_enrichments(range(len(ranked_data)))
    # substitute the counts in rankedData with the respective bootstrapping pvalues
    for index, pvalue in corrected_pvalues.items():
        ranked_data[index][1] = pvalue
        

    # rankedData: [ [es1,bp,hp], [es2,bp,hp], [es3,bp,hp] ]
    ranked_data = sorted(ranked_data, key=operator.itemgetter(1,2))
    data = filter(lambda x: x[2]<=hyper_pvalue_threshold, ranked_data)
    data = filter(lambda x: x[1]<=corrected_pvalue_threshold, data)    
    # data: [ [es1,bp,hp], [es2,bp,hp], [es3,bp,hp] ]

    # create the final output:
    client = pymongo.MongoClient()
    db = client[mongo_database]
    collection = db[ontology+"_base"]
    output = []
    i = 0
    for entry in data:
        terms = P.get_terms_of_enrichment(entry[0])
        terms_to_genes = P.get_terms_mapping(terms)
        for term in terms:
            i += 1
            try:
                query = list(collection.find({"term_id":term}, {'definition':1}))[0]
                definition = str(query['definition'])
                d = {"rank":i, "term_id" : term, 
                "term_definition" : definition, 
                "enrichment" : str(entry[0]), 
                "hyper_pvalue" : "{:.3E}".format(decimal.Decimal(entry[2])), 
                "corrected_pvalue" : entry[1], 
                "genes" : list(set(i for i in terms_to_genes[term]))
                }
                output.append(d)
            except IndexError:
                pass
    return {ontology: output}




if __name__ == '__main__':

    '''
    DESCRIPTION:

    Execution of the enrichment analysis for all the gene sets of the 20
    conserved mechanisms, as they have been defined in the previous steps. 
    '''
    
    mongodb = 'background'
    main_dir = './other_terms_analysis/'
    terms = os.listdir(main_dir+'pathway_analysis_inputs/')
    
    for term in terms[15:]:

        output_dir = main_dir+'pathway_analysis_outputs/'+term+'/'
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        files = os.listdir(main_dir+'pathway_analysis_inputs/'+term+'/')
        ontology = 'GO_P'
    
        for en, f in enumerate(files):
            print(term, en, len(files))
            label = os.path.basename(f).split('.')[0]
            tmp_species = '_'.join(label.split('_')[0:2])
           
            F = open(main_dir+'pathway_analysis_inputs/'+term+'/'+f, 'r')
            input_list = [i.replace('\n', '').replace('\r', '') for i in F.readlines()]

            ea_results = core_process(gene_ids=input_list, ontology=ontology, 
                                      organism=tmp_species, corrected_mapping=True, 
                                      n_boot=1000, hyper_pvalue_threshold=0.1, 
                                      corrected_pvalue_threshold=1, mongo_database=mongodb)

            cols = ['rank', 'term_id', 'term_definition', 'enrichment', 
                    'hyper_pvalue', 'corrected_pvalue', 'genes']    
            df = pandas.DataFrame(ea_results['GO_P'], columns=cols)
            df.to_csv(output_dir + label + '_GO_P.tsv', index=None, sep='\t')


