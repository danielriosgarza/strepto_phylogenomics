#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 15:37:56 2019

@author: meike
"""
'''Parameters for genome quality: 
    1. genome quality = good
    2. Completness of more than or equal to 90%
    3. Genome length of at least 1000
    4. Status not plasmid
    5. If status is nothing ('') than cds of at least 700
    6. It should have less or equal to 500 contigs
    7. Consistencies above or equal to 95% '''

import csv
import os

def patric_genome_quality_fields(genus, save_directory):
    '''Needs genus of interest and makes file with all species found on patric ftp server with quality related
    fields'''
#    if not os.path.exists(save_directory):
#        os.makedirs(save_directory)
    info =''.join('p3-all-genomes --eq genus,"' +genus+ '" | p3-get-genome-data --attr genome_name --attr genome_quality --attr genome_status --attr patric_cds --attr genome_length --attr checkm_completeness --attr checkm_contamination --attr contigs --attr contig_n50 --attr contig_l50 --attr coarse_consistency --attr fine_consistency > '+ save_directory)
    os.system(info)
    return save_directory

original_count=[]
genomes =[]
with open('/home/meike/strepto_phylogenomics/files/strepto_all_genomes_quality.tsv') as f:
    f_reader = csv.reader(f, delimiter="\t")
    for line in f_reader: 
        if line[0].startswith('genome.genome_id'):
            genomes.append(line)
        else:
            original_count.append(line)
            if line[2] != 'Poor' and line[3] != 'Plasmid':      #filter out all genomes with poor genome quality and plasmids
                if int(line[4]) >= 700 and line[5] != '' and int(line[5]) >= 1000: #filter cds of at least 700 and genome size of at least 1000
                    if line[6] != "" and float(line[6]) >= 90 and float(line[7]) <= 10: #filter on completness more or equal 90% and contamination less than 10
                        if int(line[8]) <= 500 and float(line[11]) >= 95 and float(line[12]): # filter on contigs and consistencies
                            genomes.append(line)

with open ('/home/meike/strepto_phylogenomics/files/strepto_genomes_quality.tsv', 'w') as f:
    for row in genomes:
        f.write('\t'.join(row) + '\n')
print(len(original_count), len(genomes))

patric_genome_quality_fields('Streptococcus', '/home/meike/phylogenomics/files/strepto_all_genomes_quality.tsv')