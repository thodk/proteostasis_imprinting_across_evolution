# -*- coding: utf-8 -*-
import pandas
import os


mongodb = 'background'
parameters_df = pandas.read_csv('./files/ensembl_parameters.tsv', sep='\t',
                                index_col=0)
ensembl_classes = parameters_df.loc[:, 'class'].tolist()


model_species_df = pandas.read_csv('./files/model_species.tsv', sep='\t', 
                                     index_col=0)
model_species = model_species_df.species.tolist()
model_species_classes = dict(zip(model_species_df.species, 
                                   model_species_df.loc[:,'class']))

                                 
define_main_dir = lambda x: os.path.abspath('./data/'+x+'/preprocessing/')+'/'
