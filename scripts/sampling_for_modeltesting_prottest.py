#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 16:53:35 2020

@author: meike
"""

'''
Searching for the best model for maximum-likelihood based phylogenetic tree. Test sets: 5 fasta files with 100 random seqs.
'''

import os
from pathlib import Path
from Bio import SeqIO
from random import sample

def random_pick(inputfile, outputfile, number=100):
    '''
    Takes random lines of a fastafile (default 100) id and seq and writes it into a new file
    '''
    with open(inputfile) as f:
        with open(outputfile, 'w') as f2:
            seqs = SeqIO.parse(f,"fasta")
            samps = ((seq.name, seq.seq) for seq in sample(list(seqs), number))
            for samp in samps:
                f2.write(">{}\n{}\n".format(*samp))            
               
path = os.getcwd()
p = Path(path)

file_path = '/home/meike/tests/Files/concat_alignments'

for i in range(1,6):
    random_pick(file_path,'/home/meike/tests/Files/prottest_testset'+ str(i)+'.fa' )

random_pick(file_path,'/home/meike/tests/Files/prottest_testset_larger.fa', number=1000)

'''
run following line on newly generated files to check for the best model:
    
java -jar /home/meiker/software/prottest-3.4.2/prottest-3.4.2.jar -i /home/meiker/phylo_tree/model_testing/prottest_testset1.fa -o /home/meiker/phylo_tree/model_testing/prottest_prediction_1.txt -all -threads 8
'''
