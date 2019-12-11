#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 09:29:24 2019

@author: meike
"""

'''Getting fields of the Outgenomes: Floricoccus, Lactococcus, Lactovum, Okadaella. 
--> Patric no genomes of Okadaella, Lactovum. Floricoccus two, Lactococcus with quality check'''

import os
import csv
import io
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
    and saving directory and makes file with all species found
    on patric ftp server'''
#    if not os.path.exists(save_directory):
#        os.makedirs(save_directory)
    fields = patric_fields_str(fields_file)
    info = ''.join('p3-all-genomes --eq genus,"' +genus+ '" | p3-get-genome-data ' + fields + ' > ' + save_directory)
    os.system(info)
    return save_directory

def convert2utf8(filename):
    '''Downloaded file from Patric was not utf-8 type. Function converts files to utf-8 encoding files.
    Replaces errors.'''
    with io.open(filename, 'r', encoding='utf-8', errors= 'replace') as f:
        text = f.read()
    # process Unicode text
    with io.open(filename, 'w', encoding='utf-8') as f:
        f.write(text)
    return filename

def thresholds (parameter_l):
    '''Uses a list to make a boxplot and returns the whisker ends that can be used as threshold. Default settings'''
    parameter_plot = plt.boxplot(parameter_l, showfliers=False)
    parameter_whis = [item.get_ydata() for item in parameter_plot['whiskers']]
    parameter_thresholds = [parameter_whis[0][1], parameter_whis[1][1]]
    return parameter_thresholds

#Look for averages of contamination to determine thresholds
contaminations =[]
with open ('/home/meike/strepto_phylogenomics/files/lactococcus_all_genome_fields.tsv') as f: 
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
            if line[contamination_index] != '':
                contaminations.append(float(line[contamination_index]))

#Outgenomes, where fields will be downloaded
genera = {'Floricoccus':'/home/meike/strepto_phylogenomics/files/floricoccus_all_genome_fields.tsv', 
          'Lactococcus':'/home/meike/strepto_phylogenomics/files/lactococcus_all_genome_fields.tsv',
          }
#file with all patric genome fields (except comments field)
fieldsfile = '/home/meike/strepto_phylogenomics/files/patric_fields'

#Get field information
for k,v in genera.items():
    patric_genome_fields(k, fieldsfile, v)
print(thresholds(contaminations))

#convert file to make it readable
convert2utf8('/home/meike/strepto_phylogenomics/files/lactococcus_all_genome_fields.tsv')


''''Quality check'''

#Same check as for Streptococcus for the Lactococcus genus
original_count=[]
genomes =[]
with open ('/home/meike/strepto_phylogenomics/files/lactococcus_all_genome_fields.tsv') as f:
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
original_genomes = len(original_count)
first_selection_count = len(genomes)
with open ('/home/meike/strepto_phylogenomics/files/lactococcus_genomes_quality.tsv', 'w') as f:
    for row in genomes:
        f.write('\t'.join(row) + '\n')

############ Missing species########
        
#Count how many species you loose after quality check
species_counter = {'sorts_afterwards' : 0, 'sorts_before' : 0}
with open ('/home/meike/strepto_phylogenomics/files/lactococcus_genomes_quality.tsv') as f:
    species_after_l =[]
    for line in f:
        line = line.strip().split('\t')
        if line[0].startswith('genome.genome_id'):
            field_i = line.index('genome.species')
        else:
            if line[field_i] not in species_after_l:
                species_after_l.append(line[field_i])
                species_counter['sorts_afterwards'] += 1
with open ('/home/meike/strepto_phylogenomics/files/lactococcus_all_genome_fields.tsv') as file:
    species_before_l =[]
    for line in file:
        line = line.strip().split('\t')
        if line[0].startswith('genome.genome_id'):
            field_i = line.index('genome.species')
        else:
            if line[field_i] not in species_before_l:
                species_before_l.append(line[field_i])
                species_counter['sorts_before'] += 1
print(len(species_before_l) - len(species_after_l))
missing_specs =[]
for species in species_before_l:
    if species not in species_after_l:
        missing_specs.append(species)
print(missing_specs)

#Create new file to check if species are reasonable left out
with open ('/home/meike/strepto_phylogenomics/files/lactococcus_all_genome_fields.tsv') as f:
    with open ('/home/meike/strepto_phylogenomics/files/lactococcus_genomes_missing_species_check.tsv', 'w') as file:       
        for line in f:
            row=[]
            line = line.strip().split('\t')
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
                species_i = line.index('genome.species')
                row.extend((line[0], line[species_i], line[quality_index], line[status_index] , line[cds_index] , line[completeness_index] , line[contamination_index] , line[coarse_con_index] , line[fine_con_index]))
                file.writelines('\t'.join(row) + '\n')
            else:
                for species in missing_specs:
                    if line[species_i] == species:
                        row.extend((line[0], line[species_i], line[quality_index], line[status_index] , line[cds_index] , line[completeness_index] , line[contamination_index] , line[coarse_con_index] , line[fine_con_index]))
                        file.write('\t'.join(row) + '\n')
#New parameters to include some species and increase diversity                      
additional_ids = []
with open ("/home/meike/strepto_phylogenomics/files/lactococcus_genomes_missing_species_check.tsv") as file:
    for line in file:
        line = line.strip().split('\t')
        if line[0].startswith('genome.genome_id'):
            quality_index = line.index('genome.genome_quality')
        else:
            if line[quality_index] == 'Good':
                additional_ids.append(line[0])
#Browse through all_genomes_file to get all information of the missing species                
with open ('/home/meike/strepto_phylogenomics/files/lactococcus_all_genome_fields.tsv') as f:
    for line in f:
          line = line.strip().split('\t')
          for identifier in additional_ids:
              if line[0] == identifier:
                  genomes.append(line)

final_genome_count =len(genomes)

#Make final file with all good (enough) quality species
with open ('/home/meike/strepto_phylogenomics/files/lactococcus_genomes_quality.tsv', 'w') as f:
    for row in genomes:
        f.write('\t'.join(row) + '\n')    

#summarize table of genome counts
lines = []       
with open('/home/meike/strepto_phylogenomics/files/summary_table.tsv') as f:
    for line in f:
        line = line.strip().split('\t')
        lines.append(line)        
        
cols = ['species','original_genomes', 'first_selection', 'final_genomes']
genome_numbers = ["lactococcus", str(original_genomes), str(first_selection_count), str(final_genome_count)]
flori = ['floricoccus', str(2), str(2), str(2)]
with open ('/home/meike/strepto_phylogenomics/files/summary_table.tsv' , 'w') as f:
    for line in lines:
        f.write('\t'.join(line) + '\n')
    f.write('\t'.join(genome_numbers) + '\n')
    f.write('\t'.join(flori) + '\n')                