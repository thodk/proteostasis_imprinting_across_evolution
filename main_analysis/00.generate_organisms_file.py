# -*- coding: utf-8 -*-


import pandas
import os


def worker_for_prokaryotes(category):
    filename = '../prokaryotes/'+category+'/preprocessing/03_species_with_rRNA_and_HSPs_annotation.tsv'
    global_df = pandas.read_csv(filename, sep='\t', index_col=0)
    global_df.loc[:,'taxonomy'] = [category]*len(global_df)

    return global_df



def worker_for_eukaryotes(category):
    filename = '../eukaryotes/data/'+category+'/preprocessing/04_species_with_genomic_annotation.tsv'
    global_df = pandas.read_csv(filename, sep='\t', index_col=0)
    global_df = global_df.loc[:, ['species_name', 'species_abbreviation', 'abbreviation', 'rRNA']]
    
    global_df.loc[:,'taxonomy'] = ['eukaryotes']*len(global_df)
    global_df.loc[:,'ensembl_mart'] = [category]*len(global_df)

    return global_df


ensembl_parameters_df = pandas.read_csv('../eukaryotes/files/ensembl_parameters.tsv',
                                        sep='\t', index_col=0)

data_frames = []
for c in ensembl_parameters_df.loc[:,'class'].tolist():
    data_frames.append(worker_for_eukaryotes(c))

data_frames.append(worker_for_prokaryotes('bacteria'))
data_frames.append(worker_for_prokaryotes('archaea'))

final_df = pandas.concat(data_frames, axis=0, sort=True, ignore_index=True)
columns = ['species_name', 'species_abbreviation', 'abbreviation',
           'taxonomy', 'ensembl_mart', 'rRNA']
final_df = final_df[columns]
final_df = final_df.drop_duplicates(subset='abbreviation', keep='last')
final_df = final_df.sort_values(by=['taxonomy', 'ensembl_mart'], 
                                ascending=[False, True])

if not os.path.exists('./files'):
    os.makedirs('./files')
final_df.to_csv('./files/species.tsv', sep='\t')


