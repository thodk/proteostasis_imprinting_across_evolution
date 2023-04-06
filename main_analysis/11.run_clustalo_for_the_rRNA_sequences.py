#!/usr/bin/python
# -*- coding: utf-8 -*-
import pandas
import subprocess
import os
import shutil


if __name__ == '__main__':
    
    '''
    DESCRIPTION:
        
    This scripts calls the 'clustalo' program to run the analysis and creates a 
    distance matrix file for the rRNA sequences ('final_distance_matrix.tsv').
    
    The final output as well intermediate files are stored in the new 'rRNA_analysis'
    directory.
    '''
    
    df = pandas.read_csv('./files/species.tsv', sep='\t')
    mapping = dict(zip(df.rRNA.tolist(), df.species_name.tolist()))

    fasta_files = ['./rRNAs_analysis/sequences/'+i for i in os.listdir('./rRNAs_analysis/sequences/')]
    with open('./rRNAs_analysis/final_seq.fasta', 'wb') as wfd:
        for f in fasta_files:
            with open(f,'rb') as fd:
                shutil.copyfileobj(fd, wfd)
    i = './rRNAs_analysis/final_seq.fasta'
    o = './rRNAs_analysis/distance_matrix.txt'

    p = subprocess.Popen('clustalo -i '+i+' --infmt fa --full --force --distmat-out '+o+' -o ./rRNAs_analysis/msa.txt', 
                         shell=True, stdout=subprocess.PIPE) 
    p.communicate(o)

    raw_distance_matrix = open(o, 'r')
    next(raw_distance_matrix)
    data = []
    for line in raw_distance_matrix:
        line = line.replace('\n', '')
        if line.startswith('ENA'):
            seq_id = line.split('|')[1].split('.')[0]
        else:
            seq_id = line.split('.')[0]
        fields = list(filter(lambda x: x != '', line.split(' ')))
        fields = [mapping[seq_id]] + fields[1:]    
        data.append(fields)
    raw_distance_matrix.close()
    
    names = [k[0] for k in data]
    FO = open('./rRNAs_analysis/final_distance_matrix.tsv', 'w')
    FO.write('\t' + '\t'.join(names)+'\n')
    for i in data:
        FO.write('\t'.join(i)+'\n')
    FO.close()
