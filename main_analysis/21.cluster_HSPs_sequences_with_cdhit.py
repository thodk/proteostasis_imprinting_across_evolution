#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
import subprocess


def worker(hsp):
    files = os.listdir('./HSPs_analysis/'+hsp+'/sequences/')
    initial_fasta_files = ['./HSPs_analysis/'+hsp+'/sequences/'+i for i in files]
    main_dir = './HSPs_analysis/'+hsp+'/after_clustering/'
    
    if not os.path.exists(main_dir):
        os.makedirs(main_dir)
        
    for f in initial_fasta_files:
        organism = os.path.basename(f).split('.')[0]
        F = open(f, 'r')
        c=0
        for line in F:
            if line.startswith('>'):
                c+=1
            else:
                pass
        o = main_dir + organism + '.fasta'
        if c > 1:
            c = '0.95'
            p = subprocess.Popen('cd-hit -i '+f+' -o '+o+' -c '+c , shell=True, 
                                 stdout=subprocess.PIPE)
        else:
            p = subprocess.Popen('cp '+f+' '+o, shell=True, stdout=subprocess.PIPE)
        p.communicate(o)

    files = os.listdir('./HSPs_analysis/'+hsp+'/after_clustering/')
    to_remove = ['./HSPs_analysis/'+hsp+'/after_clustering/'+i for i in files if i.endswith('.clstr')]
    for f in to_remove:
        os.remove(f)


if __name__ == '__main__':
    
    '''
    DESCRIPTION:
        
    Many species have plenty of hsp40 and hsp70 sequences. However, a single 
    consensus sequence should be created for each species. In order to
    reduce the redundancy of these sequences, this script calls the 'CDHIT'
    software to cluster all the HSP of the same family for each species, with 
    similarity threshold at 95%.
    
    If the output clusters are more than one, then they will be used to create 
    a consensus sequence, based on an HMM model of these clusters (see the next
    script '22.generate_consensus_HSPs_sequences.py'). The final goal is to 
    create a single hsp40 and hsp70 sequence for each species.
    '''
    
    worker('hsp40')
    worker('hsp70')