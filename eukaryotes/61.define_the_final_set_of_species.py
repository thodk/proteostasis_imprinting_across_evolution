#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas
import sys
sys.path.append('../')
from constant_variables import define_main_dir
from constant_variables import ensembl_classes

'''
def worker(ensembl_class):
    
    preprocessing_dir = define_main_dir(ensembl_class)
    pn_lists_dir = os.path.abspath('./data/'+ensembl_class+'/data/PN_gene_lists')+'/'
    
    df = pandas.read_csv(preprocessing_dir+'04_species_with_genomic_annotation.tsv', sep='\t', index_col=0) 
    species_dict = dict(zip(df.abbreviation, df.species_abbreviation))
    files = [pn_lists_dir+i for i in os.listdir(pn_lists_dir)]
    sizes = []
    for f in files:
        F = open(f, 'r')
        species = f.split('/')[-1].replace('.txt', '')
        sizes.append([species, len(F.readlines())])
        F.close()
    
    sizes_df = pandas.DataFrame(sizes)

    sizes_df.loc[:,'logsizes'] = sizes_df.iloc[:,1]
    q25 = numpy.percentile(sizes_df.loc[:,'logsizes'], 25)
    q75 = numpy.percentile(sizes_df.loc[:,'logsizes'], 75)
    median = numpy.percentile(sizes_df.loc[:,'logsizes'], 50)
    mad = numpy.median(numpy.absolute(sizes_df.loc[:,'logsizes'] - median))

    i1 = sizes_df.loc[(sizes_df.loc[:,'logsizes'] > q75+3*mad) |
                      (sizes_df.loc[:,'logsizes'] < q25-3*mad)]
    #for i, row in i1.iterrows():
    #    print(species_dict[row[0]])
    #print(i1)
    print(i1)    
'''
    
if __name__ == '__main__':

    dfs = []
    for ensembl_class in ensembl_classes:    
        f = define_main_dir(ensembl_class)+'04_species_with_genomic_annotation.tsv'
        df = pandas.read_csv(f, sep='\t', index_col=0)
        df.loc[:, 'ensembl_class'] = [ensembl_class]*len(df)
        dfs.append(df)
    global_df = pandas.concat(dfs, axis=0, ignore_index=True, sort=True)
    global_df.to_csv('./files/final_species_set_for_analysis.tsv', sep='\t')