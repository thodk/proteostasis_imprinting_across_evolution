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
        o = main_dir + organism + '.fa'
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
    worker('hsp40')
    worker('hsp70')