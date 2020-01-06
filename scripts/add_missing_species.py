#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 13:43:32 2019

@author: meike
"""

'''Species comparison to check if one misses and check if it is reasonably left out.'''
                
'''Comparison of species of all genomes to filtered genomes'''

import csv

species_counter = {'sorts_afterwards' : 0, 'sorts_before' : 0}
with open ('/home/meike/strepto_phylogenomics/files/strepto_genomes_quality.tsv') as f:
    species_after_l =[]
    for line in f:
        line = line.strip().split('\t')
        if line[0].startswith('genome.genome_id'):
            field_i = line.index('genome.species')
        else:
            if line[field_i] not in species_after_l:
                species_after_l.append(line[field_i])
                species_counter['sorts_afterwards'] += 1
with open ('/home/meike/strepto_phylogenomics/files/strepto_all_genome_fields.tsv') as file:
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

with open ('/home/meike/strepto_phylogenomics/files/strepto_all_genome_fields.tsv') as f:
    with open ('/home/meike/strepto_phylogenomics/files/strepto_genomes_missing_species_check.tsv', 'w') as file:
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
         