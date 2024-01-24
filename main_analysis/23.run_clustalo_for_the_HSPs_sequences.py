# -*- coding: utf-8 -*-

import os
import subprocess
import shutil
import pandas

def worker(hsp):
    
    main_dir = './HSPs_analysis/'+hsp+'/consensus_sequences/'
    main_dir2 = './HSPs_analysis/'+hsp+'/'
    fasta_files = [main_dir+i for i in os.listdir(main_dir)]
    
    with open(main_dir2+'final_seq.fa', 'wb') as wfd:
        for f in fasta_files:
            with open(f, 'rb') as fd:
                shutil.copyfileobj(fd, wfd)
    
    i = main_dir2+'final_seq.fa'
    o = main_dir2+'distance_matrix.txt'
    p = subprocess.Popen('clustalo -i '+i+' --infmt fa --full --force --distmat-out '+o+' -o '+main_dir2+'msa.txt', 
                         shell=True, stdout=subprocess.PIPE) 
    p.communicate(o)
    
    raw_distance_matrix = open(o, 'r')
    next(raw_distance_matrix)
    data = []
    for line in raw_distance_matrix:
        line = line.replace('\n', '')
        fields = list(filter(lambda x: x != '', line.split(' ')))    
        data.append(fields)
    raw_distance_matrix.close()
    
    df = pandas.read_csv('./files/species.tsv', sep='\t')
    species_mapping = dict(zip(df.abbreviation, df.species_name))    
    
    names = [species_mapping[k[0]] for k in data]
    FO = open('./HSPs_analysis/'+hsp+'/final_distance_matrix.tsv', 'w')
    FO.write('\t' + '\t'.join(names)+'\n')
    for i in data:
        FO.write('\t'.join([species_mapping[i[0]]]+i[1:])+'\n')
    FO.close()



if __name__ == '__main__':
    
    '''
    DESCRIPTION:
        
    This scripts calls the 'clustalo' program to run the analysis and creates a 
    distance matrix file for the HSP sequences of the same family (hsp40 and
    hsp70).
    
    The final outputs ('final_distance_matrix.tsv'), as well intermediate files, 
    are stored in the respective subfolders in 'HSPs_analysis'.
    '''
    
    worker('hsp40')
    worker('hsp70')