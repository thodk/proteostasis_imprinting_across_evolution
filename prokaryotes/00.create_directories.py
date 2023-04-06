#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
from constant_variables import prokaryotic_classes



if __name__ == '__main__':

    os.makedirs('./files/', exist_ok=True)

    for c in prokaryotic_classes:
        os.makedirs('./data/'+c+'/data/rRNA_fasta', exist_ok=True)
        os.makedirs('./data/'+c+'/data/hsp40_fasta', exist_ok=True)
        os.makedirs('./data/'+c+'/data/hsp70_fasta', exist_ok=True)
        os.makedirs('./data/'+c+'/preprocessing/files', exist_ok=True)
