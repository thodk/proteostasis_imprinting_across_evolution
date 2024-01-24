#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import pandas
import numpy
import pymongo
import sys
from scipy import stats
sys.path.append('../')
from core_functions import mongodb
from core_functions import construct_graph_from_mongo
from core_functions import get_mapping_from_mongo
import core_classes


if __name__ == '__main__':
    
    '''
    DESCRIPTION:

    This script generates the inputs gene sets for the pathway analysis of the
    20 conserved biological processes. The whole workflow is described in the 
    section 2.6.1 of the manuscript.
    '''
    
    def get_pvalue(size):
        gene_set = numpy.random.choice(ref_gene_set, size=size, replace=False)
        sample_mapping = {}
        for gene in gene_set:
            terms = []
            try:
                terms = background_mapping[gene]
            except KeyError:
                pass
            sample_mapping.update({gene:terms})
        sample_pool = core_classes.Pool()
        sample_pool.set_pool_content(sample_mapping, pool_content="values")
        Pt = background_pool.get_population()
        Ps = len(ref_gene_set)
        S = sample_pool.get_population()
        x = len(gene_set)
        distribution = stats.hypergeom(Pt, Ps, S)
        pvalue = distribution.sf(x) + distribution.pmf(x)
        return pvalue
    
    G = construct_graph_from_mongo('GO_P', mongo_database=mongodb)
    selected_terms = ['GO:0006457', 'GO:0051641', 'GO:0008033', 'GO:0006974', 
                      'GO:0009260', 'GO:0008104', 'GO:0051252', 'GO:0032259',
                      'GO:0006629', 'GO:0006355', 'GO:0022607', 'GO:0033554', 
                      'GO:0015031', 'GO:0006096', 'GO:0060255', 'GO:0016310',
                      'GO:0006605', 'GO:0006310', 'GO:0006281', 'GO:0046034']

    files = os.listdir('./PN_analysis/pathway_analysis_outputs/')
    statistics_dict = {}
    for f in files:
        df = pandas.read_csv('./PN_analysis/pathway_analysis_outputs/'+f, sep='\t')
        tmp_species = '_'.join(f.split('_')[0:2])
        power = numpy.mean([-numpy.log10(p) for p in sorted(df.hyper_pvalue.tolist())[0:10]])
        statistic = numpy.power(10, -power)
        statistics_dict.update({tmp_species:statistic})


    for tmp_species in list(statistics_dict.keys()):
        for term in selected_terms[:10]:
            client = pymongo.MongoClient()
            db = client['background']
            collection = db[tmp_species+'_mapping_GO_P_corrected']
            results = collection.aggregate([
                {'$match': {'term_accession': {'$in': [term]}}},
                {'$group': {'_id': "$term_accession", 
                            'list': {'$addToSet': '$gene_symbol'}}}
                ])
            for entry in results:
                gene_set = entry['list']
            
            background_mapping = get_mapping_from_mongo('GO_P', tmp_species,
                                                        keys="genes",
                                                        mongo_database=mongodb)
            background_pool = core_classes.Pool()
            background_pool.set_pool_content(background_mapping, pool_content="keys")
            thr = statistics_dict[tmp_species]
            ref_gene_set = gene_set[:]
            sizes = []
            if len(ref_gene_set) <= 10:
                sizes.append(len(ref_gene_set))
            else:
                for n in range(30):
                    pvalue = 0
                    size = len(ref_gene_set)
                    step = int(numpy.ceil(len(ref_gene_set)/2))
                    while step > 1:
                        if pvalue > thr:
                            step = int(numpy.ceil(step/2))
                            tmp_pvalue = get_pvalue(size+step)
                            if tmp_pvalue < thr:
                                size = size+step
                            else:
                                size = size+2*step
                            pvalue=0
        
                        while pvalue < thr:
                            tmp_size = max([size - step, 1])
                            pvalue = get_pvalue(tmp_size)
                            if pvalue < thr:
                                size = tmp_size
                                step = max([int(numpy.ceil(step/2)), 1])
                            else:
                                size = tmp_size
                    sizes.append(size)

            if numpy.mean(sizes) <= 10 and len(ref_gene_set) > 10:
                sizes = [10] * len(sizes)
            
            if len(ref_gene_set) > 2*numpy.mean(sizes):
                N = 10
            elif len(ref_gene_set) > numpy.mean(sizes):
                N = 5
            else:
                N = 1
            term_definition = G.get_entry_obj(term).definition.replace(' ', '_')
            tmp_dir = './other_terms_analysis/pathway_analysis_inputs/'+term_definition+'/'
            if not os.path.exists(tmp_dir):
                os.makedirs(tmp_dir)
            for n in range(N):
                if len(sizes) > 1:
                    s = int(numpy.random.normal(loc=numpy.mean(sizes), scale=numpy.std(sizes), size=1)[0])
                    tmp_gene_set = numpy.random.choice(ref_gene_set, size=s, replace=False)
                else:
                    tmp_gene_set = ref_gene_set[:]
                F = open(tmp_dir+tmp_species+'_'+str(n)+'.txt', 'w')
                F.write('\n'.join(tmp_gene_set))
                F.close()
                
        
