# -*- coding: utf-8 -*-

import subprocess
import os


def worker(hsp):

    files = os.listdir('./HSPs_analysis/'+hsp+'/after_clustering/')
    clustering_fasta_files = ['./HSPs_analysis/'+hsp+'/after_clustering/'+i for i in files]
    
    main_dir = './HSPs_analysis/'+hsp+'/consensus_sequences/'
    hmm_dir = './HSPs_analysis/'+hsp+'/hmm_draft/'

    if not os.path.exists(main_dir):
        os.makedirs(main_dir)

    if not os.path.exists(hmm_dir):
        os.makedirs(hmm_dir)
        
    for f in clustering_fasta_files:
        species = os.path.basename(f).split('.')[0]
        F = open(f, 'r')
        c=0
        for line in F:
            if line.startswith('>'):
                c+=1
            else:
                pass
 
        if c>1:
            i = f
            o = hmm_dir+species+'_msa.txt'
            p = subprocess.Popen('clustalo -i '+i+' -t Protein --infmt fa -o '+o, 
                                 shell=True, stdout=subprocess.PIPE) 
            p.communicate(o)
        
            i = hmm_dir+species+'_msa.txt'
            o = hmm_dir+species+'_hmm.txt'
            p = subprocess.Popen('hmmbuild --amino '+o+' '+i, shell=True, 
                                 stdout=subprocess.PIPE) 
            p.communicate(o) 
        
            i = hmm_dir+species+'_hmm.txt'
            o = hmm_dir+species+'_consensus.fasta'
            p = subprocess.Popen('hmmemit -c -o '+o+' '+i, shell=True, 
                                 stdout=subprocess.PIPE) 
            p.communicate(o) 
        else:
            i = f
            o = hmm_dir+species+'_consensus.fasta'
            p = subprocess.Popen('cp '+i+' '+o, shell=True, stdout=subprocess.PIPE)
            p.communicate(o)       

        FR = open(hmm_dir+species+'_consensus.fasta', 'r') 
        FO = open(main_dir+species+'.fasta', 'w') 
        entry = None
        seq = str()
        for line in FR:
            if line.startswith('>'):
                entry = line.replace('\n', '')
            else:
                seq = seq + line.replace('\n', '')   
        FO.write('>'+species+'\n'+seq+'\n')
        FO.close()
        FR.close()



if __name__ == '__main__':
    
    '''
    DESCRIPTION:
        
    If the output of CDHIT contains more than one clusters, this script creates
    a multiple sequence alignment (MSA) and the subsequent hidden markov model
    (HMM) for these sequences (using the 'clustalo' and 'hmmbuild' functions).
    Then, another function of HMMER3 is used ('hmmemit') to produce a consensus 
    sequence for each HMM. In this way, all the species will be annotated with
    a single sequence for each HSP family.
    '''
    
    worker('hsp40')
    worker('hsp70')