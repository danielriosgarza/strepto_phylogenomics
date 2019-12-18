#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 14:48:37 2019

@author: meike
"""

'''Making a File with all genomes that are of good quality'''
import csv
#import matplotlib.pyplot as plt

def thresholds (parameter_l):
    '''Uses a list of numbers to determine thresholds for outliers. 
    Gives lower bound with q1-1.5*iqr and upper bound q3+1.5*iqr'''
    q1 = sorted(parameter_l)[int(len(parameter_l) * .25)]
    q3 = sorted(parameter_l)[int(len(parameter_l) * .75)]
    iqr = q3 - q1
    lower_bound = q1 -(1.5 * iqr) 
    upper_bound = q3 +(1.5 * iqr)
    parameter_thresholds = [lower_bound, upper_bound]
    return parameter_thresholds

'''Find appropiate thresholds for genome quality check'''
lengths = []
contigs =[]
contaminations = []
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
            contigs.append(int(line[contigs_index]))
            if line[contamination_index] != '':
                contaminations.append(float(line[contamination_index]))
print(thresholds(lengths), thresholds(contigs),thresholds(contaminations))
#plt.boxplot(lengths)                
'''Parameters for genome quality: 
    1. genome quality = good
    2. Completness of more than or equal to 90%
    3. Status not plasmid
    /4. If status is nothing ('') than cds of at least 700
    /5. It should have less or equal to 150 contigs
    6. Consistencies above or equal to 95%'''
original_count=[]
genomes =[]
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
            genomes.append(line)
        else:
            original_count.append(line)
            if line[quality_index] == 'Good' and line[status_index] != 'Plasmid':     #filter out all genomes with poor genome quality and plasmids
                if line[completeness_index]!= '' and float(line[completeness_index]) >= 90:  #filter on completness
                    if line[contamination_index] != "" and float(line[contamination_index]) <= 2: #filter on contamination
                        if float(line[coarse_con_index]) >= 95 and float(line[fine_con_index]) >= 95: # filter on consistencies
                            if int(line[cds_index]) >= 700:    #filter out probably wrong anntoated plasmids
                                genomes.append(line)
                                
print(len(original_count), len(genomes))

additional_ids = []
with open ("/home/meike/strepto_phylogenomics/files/strepto_genomes_missing_species_check.tsv") as file:
    f_reader = csv.reader(file, delimiter="\t")
    for line in f_reader:
        if line[0].startswith('genome.genome_id'):
            cds_index = line.index('genome.patric_cds')
            status_index = line.index('genome.genome_status')
            quality_index = line.index('genome.genome_quality')
            completeness_index = line.index('genome.checkm_completeness')
        else:
            if int(line[cds_index]) >= 700 and line[status_index] != 'Plasmid':
                if line[quality_index] == 'Good' or line[quality_index] == '':
                    if line[completeness_index] != "" and float(line[completeness_index]) >= 90:
                        additional_ids.append(line[0])
                    elif line[completeness_index] == "":
                        additional_ids.append(line[0])
print(len(original_count), len(genomes))
original_genomes = len(original_count)-1
first_selection_count = len(genomes)-1

#Browse through all_genomes_file to get all information of the missing species                
with open ('/home/meike/strepto_phylogenomics/files/strepto_all_genome_fields.tsv') as f:
    for line in f:
          line = line.strip().split('\t')
          for identifier in additional_ids:
              if line[0] == identifier:
                  genomes.append(line)


final_genome_count =len(genomes)-1

with open ('/home/meike/strepto_phylogenomics/files/strepto_genomes_quality.tsv', 'w') as f:
    for row in genomes:
        f.write('\t'.join(row) + '\n')

#summarize table of genome counts
cols = ['species','original_genomes', 'first_selection', 'final_genomes']
genome_numbers = ["streptococcus", str(original_genomes), str(first_selection_count), str(final_genome_count)]
with open ('/home/meike/strepto_phylogenomics/files/summary_table.tsv' , 'w') as f:
        f.write('\t'.join(cols) + '\n' + '\t'.join(genome_numbers) + '\n')