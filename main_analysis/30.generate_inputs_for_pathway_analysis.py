# -*- coding: utf-8 -*-

import pandas
import pymongo
import os


def worker_for_prokaryotes(category):
    
    tmp_species_df = species_df.loc[species_df.taxonomy == category, :]
    species_dict = dict(zip(tmp_species_df.species_name, tmp_species_df.abbreviation))
    
    client = pymongo.MongoClient()
    db = client['background']    

    details = []    
    for species, abbreviation in species_dict.items():
        collection = db[abbreviation + '_mapping_GO_P_corrected']
        results = collection.find({'gene_symbol':{'$in':prokaryotic_genes}})
        tmp_genes = []
        for i in results:
            tmp_genes.append(i['gene_symbol'])
        tmp_genes = list(set(tmp_genes))
        F = open('./PN_analysis/pathway_analysis_inputs/'+abbreviation+'.txt', 'w')
        F.write('\n'.join(tmp_genes))
        F.close()
        
        details.append([species, abbreviation, category, len(tmp_genes)])
        
    return details
    

def worker_for_eukaryotes(category):
    
    tmp_species_df = species_df.loc[species_df.ensembl_mart == category, :]
    species_dict = dict(zip(tmp_species_df.species_abbreviation, tmp_species_df.abbreviation))
    details = []
    
    ref_dir = '../eukaryotes/data/'+category+'/preprocessing/annotation/model_species/corrected_gene_lists/'
    reference_files = os.listdir(ref_dir)
    reference_mapping = {}
    for reference_file in reference_files:
        reference_species = reference_file.split('.')[0]        
        df = pandas.read_csv(ref_dir+reference_file, sep='\t')
        if reference_species == 'xtropicalis':
            gene_symbols = df.iloc[:,2].tolist()
            gene_symbols = [j for j in gene_symbols if not pandas.isnull(j)]
            gene_symbols = [i for j in gene_symbols for i in str(j).split(' ')]
        else:
            gene_symbols = df.iloc[:,1].tolist()
        gene_symbols = filter(lambda x: not pandas.isnull(x), gene_symbols)
        reference_mapping.update({reference_species:gene_symbols})
        print reference_species, len(gene_symbols)
        
    main_dir = '../eukaryotes/data/'+category+'/preprocessing/annotation/homologies/'
    for species, abbreviation in species_dict.items():
        if species in reference_eukaryotes:
            continue
        final_ensembl_ids = []
        for ref_species, gene_symbols in reference_mapping.items():
            tmp_df = pandas.read_csv(main_dir+species+'_'+ref_species+'_homologs.tsv',
                                     sep='\t', index_col=0)
            tmp_df = tmp_df.loc[tmp_df.external_gene_name.isin(gene_symbols), :]
            if category == 'ensembl':
                homologies = tmp_df.loc[:,species+'_homolog_ensembl_gene'].tolist()
            else:
                try:
                    homologies = tmp_df.loc[:,species+'_eg_homolog_ensembl_gene'].tolist()
                except KeyError:
                    homologies = tmp_df.loc[:,species+'_eg_homolog_associated_gene_name'].tolist()
            homologies = filter(lambda x: not pandas.isnull(x), homologies)
            final_ensembl_ids = final_ensembl_ids + homologies
        final_ensembl_ids = list(set(final_ensembl_ids))
        
        F = open('./PN_analysis/pathway_analysis_inputs/'+abbreviation+'.txt', 'w')
        F.write('\n'.join(final_ensembl_ids))
        F.close()
        
        details.append([species, abbreviation, 'eukaryotes', len(final_ensembl_ids)])

    for reference_species, gene_symbols in reference_mapping.items():
        F = open('./PN_analysis/pathway_analysis_inputs/'+species_dict[reference_species]+'.txt', 'w')
        F.write('\n'.join(gene_symbols))
        F.close()
        details.append([reference_species, species_dict[reference_species], 
                        'eukaryotes', len(gene_symbols)])

    return details
    



if __name__ == '__main__':

    species_df = pandas.read_csv('./files/species.tsv', sep='\t')
    prokaryotic_genes = pandas.read_csv('./files/prokaryotic_genes_for_proteostasis.tsv',
                                        sep='\t', index_col=0).gene_symbol.tolist()
    if not os.path.exists('./PN_analysis/pathway_analysis_inputs'):
        os.makedirs('./PN_analysis/pathway_analysis_inputs')
    
    ensembl_parameters_df = pandas.read_csv('../eukaryotes/files/ensembl_parameters.tsv',
                                            sep='\t', index_col=0)
    reference_eukaryotes = []
    for c in ensembl_parameters_df.loc[:,'class'].tolist():
        ref_dir = '../eukaryotes/data/'+c+'/preprocessing/annotation/model_species/corrected_gene_lists/'
        ref_files = os.listdir(ref_dir)
        reference_eukaryotes = reference_eukaryotes + [i.split('.')[0] for i in ref_files]
    
    details = []
    for c in ensembl_parameters_df.loc[:,'class'].tolist():
        details.append(worker_for_eukaryotes(c))
    details.append(worker_for_prokaryotes('bacteria'))
    details.append(worker_for_prokaryotes('archaea'))
    details = [i for j in details for i in j]
    details_df = pandas.DataFrame(details, columns=['species_name', 'abbreviation', 
                                                    'taxonomy', 'gene_list_length'])
    details_df.to_csv('./files/pathway_analysis_inputs_details.tsv', sep='\t')


