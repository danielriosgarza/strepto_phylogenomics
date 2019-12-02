#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 09:24:02 2019

@author: meike
"""
'''Makes tsv file of all patric genome fields of interest'''
import os
import csv
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def patric_fields_str(fields_file):
    '''Needs file with fields of interest. Generates string that can be used in os system to retrieve fields.'''
    fields_l =[]
    fields =''
    with open (fields_file) as f:
        for line in f:
            line = line.strip()
            line = ''.join('--attr ' + line)
            fields_l.append(line)
    fields = ' '.join(fields_l)
    return fields      
def patric_genome_fields(genus, fields_file, save_directory):
    '''Needs genus and fields (list,  field of interest) 
    and saving directory and makes file with all species found \
    on patric ftp server'''
#    if not os.path.exists(save_directory):
#        os.makedirs(save_directory)
    fields = patric_fields_str(fields_file)
    info = ''.join('p3-all-genomes --eq genus,"' +genus+ '" | p3-get-genome-data ' + fields + ' > ' + save_directory)
    os.system(info)
    return save_directory

patric_genome_fields("Streptococcus", '/home/meike/strepto_phylogenomics/files/patric_fields', '/home/meike/strepto_phylogenomics/files/strepto_all_genome_fields.tsv') 


def thresholds (parameter_l):
    '''Uses a list to make a boxplot and returns the whisker ends that can be used threshold. Default settings'''
    parameter_plot = plt.boxplot(parameter_l)
    parameter_whis = [item.get_ydata() for item in parameter_plot['whiskers']]
    parameter_thresholds = [parameter_whis[0][1], parameter_whis[1][1]]
    return parameter_thresholds

lengths = []
contigs =[]
with open ('/home/meike/strepto_phylogenomics/files/strepto_all_genome_fields.tsv') as f:
    f_reader = csv.reader(f, delimiter="\t")
    for line in f_reader:
        if line[0].startswith('genome.genome_id'):  #find indexes to use columns as filters
            quality_index = line.index('genome.genome_quality')
            status_index = line.index('genome.genome_status')
            cds_index = line.index('genome.patric_cds')
            length_index = line.index('genome.genome_length')
            completeness_index = line.index('genome.checkm_completeness')
            contamination_index = line.index('genome.checkm_contamination')
            contigs_index= line.index('genome.contigs')
            coarse_con_index = line.index('genome.coarse_consistency')
            fine_con_index = line.index('genome.fine_consistency')
        else:
            if line[length_index] != '':
                lengths.append(int(line[length_index]))
thresholds(lengths)
