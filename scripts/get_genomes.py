#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 14:03:50 2019

@author: meike
"""

'''Get patric genomes of interest. Make a file containing patric_ids of interest.
Loop through ids and download .fna files in chosen directory'''

import csv
import os

def get_patric_genomes(patric_ids_file, save_directory):
    '''Needs a file with patric_ids of interest and a saving directory. Downloads the genomes (fna files)
    of the ids.'''
    if not os.path.exists(save_directory):
        os.makedirs(save_directory)
    with open (patric_ids_file) as f: 
        for identifier in f:
            identifier = identifier.strip()
            url = ''.join('ftp://ftp.patricbrc.org/genomes/'+identifier+'/' +identifier+'.fna')
            os.system("wget " + url + " -P" + save_directory)
        return save_directory

ids=[]
with open ('/home/meike/strepto_phylogenomics/files/strepto_genomes_quality.tsv') as f:
    f_reader = csv.reader(f, delimiter="\t")
    for line in f_reader:
        if not line[0].startswith('genome.genome_id'):
            ids.append(line[0])
with open ('/home/meike/strepto_phylogenomics/files/strepto_genomes_ids.txt', 'w') as f:
    for identifier in ids:
        f.write(identifier+'\n')