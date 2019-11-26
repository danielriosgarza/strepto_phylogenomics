#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 15:37:56 2019

@author: meike
"""
'''Parameters for genome quality: 
    1. genome quality = good
    2. Completness of more than or equal to 90%
    3. Genome length at least 1000
    4. Status not plasmid
    5. if status nothing than cds of at least 700 
    6. Consistencies above or equal to 95% '''

import csv
import numpy as np


count=[]
genomes =[]
with open('/home/meike/strepto_phylogenomics/files/strepto_all_genomes.tsv') as f:
    f_reader = csv.reader(f, delimiter="\t")
    for line in f_reader:
        if line[0].startswith('genome.genome_id'):
            genomes.append(line)
        else:
            info = line
            count.append(line)
            if info[2] != 'Poor' and info[3] != 'Plasmid':      #filter out all genomes with poor genome quality and plasmids
                if int(info[4]) >= 700 and info[5] != '' and int(info[5]) >= 1000: #filter cds of at least 700 and genome size of at least 1000
                    if info[6] is not "" and float(info[6]) >= 90 and float(info[7]) <= 10: #filter on completness more or equal 90% and contamination less than 10
                        if float(info[11]) >= 95 and float(info[12]): # filter on consistencies
                            genomes.append(info)
with open ('/home/meike/strepto_phylogenomics/files/strepto_genomes.tsv', 'w') as f:
    for row in genomes:
        f.write('\t'.join(row) + '\n')
print(len(count), len(genomes))