# -*- coding: utf-8 -*-

import pandas
import re
import os
from constant_variables import define_main_dir
from constant_variables import model_species_df, model_species


def worker(ensembl_class, tmp_df):
    
    no_data = {}
    main_dir = define_main_dir(ensembl_class)+'annotation/'
    
    for reference_species in tmp_df.species.tolist():

        print ensembl_class, reference_species
        
        ref_df = pandas.read_csv(main_dir+'model_species/corrected_gene_lists/'+reference_species+'.tsv', 
                                 sep='\t')
        bool1 = ref_df.ensembl_ids.isna()
        bool2 = ref_df.loc[:,'gene_symbols_(primaries)'].isna()
        ref_df = ref_df.loc[(~bool1) | (~bool2)]
        
        ref_ensembl_ids = ref_df.ensembl_ids.tolist()
        ref_ensembl_ids = [re.split('\s+', i) for i in ref_ensembl_ids if not pandas.isnull(i)]
        ref_ensembl_ids = [i for j in ref_ensembl_ids for i in j]

        ref_gene_symbols = ref_df.loc[:,'gene_symbols_(primaries)']
        ref_gene_symbols = ref_gene_symbols.tolist()
        ref_gene_symbols = [i for i in ref_gene_symbols if not pandas.isnull(i)]
        
        homologies_files = os.listdir(main_dir+'homologies/')
        homologies_files = filter(lambda x: re.search(reference_species, x), homologies_files)
        
        for f in homologies_files:
            fields = f.split('_')
            if fields[0] == reference_species:
                continue
            if fields[1] != reference_species:
                continue
            if fields[0] in model_species:
                continue
            f = main_dir+'homologies/'+f
            tmp_species = fields[0]
            try:
                tmp_df = pandas.read_csv(f, sep='\t', index_col=0)
                bool1 = tmp_df.ensembl_gene_id.isin(ref_ensembl_ids) 
                bool2 = tmp_df.external_gene_name.isin(ref_gene_symbols)
                bool3 = tmp_df.iloc[:,2].isna()
                tmp_df = tmp_df.loc[((bool1) | (bool2)) & (~bool3),:]
                
                tmp_gene_symbols = tmp_df.external_gene_name.tolist()
                tmp_ensembl_ids = tmp_df.ensembl_gene_id.tolist()
                tmp_ref_df = ref_df.loc[(ref_df.loc[:,'gene_symbols_(primaries)'].isin(tmp_gene_symbols)) | 
                                        (ref_df.ensembl_ids.str.contains('|'.join(tmp_ensembl_ids) )),]
                
                ratio = round(len(tmp_ref_df) / float(len(ref_df)), 3)
                if ratio < 0.5:
                    no_data.setdefault(tmp_species, []).append(reference_species)
                
            except IOError:
                no_data.setdefault(tmp_species, []).append(reference_species)
    
    return no_data

'''
DESCRIPTION:

That script reads the homologies between the reference species and all the others.
Then for each new species, it calculates the ratio of proteostasis genes, 
defined for the reference species, which have homologies. If the ratio is lower
than 50% then the examined species is labeled with 'lack of annotation' and it 
will be excluded from the analysis in the following steps.
'''

if __name__ == "__main__":

	groups = model_species_df.groupby(by='class')

	for c, df in groups:
	    data = worker(c, df)
	    F = open('./data/'+c+'/preprocessing/files/species_with_lack_of_annotation.tsv', 'w')
	    F.write('\tspecies\treference_species\tpass\n')
	    i = 0    
	    for o, v in data.items():
		for ref_o in v:
		    F.write(str(i) + '\t' + o + '\t' + ref_o + '\t' + str(len(v) != len(df))+'\n')
		    i+=1
	    F.close()

