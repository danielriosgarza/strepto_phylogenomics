#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 09:41:25 2019

@author: meike
"""

'''Makes bash file to retrieve genomes from Patric. Saving directory is octarine server'''
#wget ftp://ftp.patricbrc.org/genomes/1074061.3/1074061.3.fna -O /home/meiker/git/genomes/streptocuccus_000001.fna
import csv

def identifier2bash (patric_ids_file):
    '''Needs file with identifiers and returns list with wget for all identifiers'''
    bash_lines = []
    with open (patric_ids_file) as f:
        for identifier in f:
            identifier = identifier.strip()
            bash_lines.append(''.join('wget ftp://ftp.patricbrc.org/genomes/'+identifier+'/' +identifier+'.fna'))
    return bash_lines

#Make file with all Patric ids
ids =[]
with open ('/home/meike/strepto_phylogenomics/files/strepto_genomes_quality.tsv') as f:
    f_reader = csv.reader(f, delimiter="\t")
    for line in f_reader:
        if not line[0].startswith('genome.genome_id'):
            ids.append(line[0])        
with open ('/home/meike/strepto_phylogenomics/files/strepto_genomes_ids.txt', 'w') as f:
    for identifier in ids:
        f.write(identifier+'\n')

bash_lines=identifier2bash('/home/meike/strepto_phylogenomics/files/strepto_genomes_ids.txt')

#generate ids for own database later. Containing in total 6 digits
db_id = [] 
strepto_db_id =[]
for i in range(1, len(bash_lines)+1):
    db_id.append(' -O /home/meiker/git/genomes/streptococcus_' + "%05d" % i + '.fna') #id for bash file
    strepto_db_id.append("streptococcus_"+"%05d" % i) #id for databse table
#make bash file to get partic genomes
with open ('/home/meike/strepto_phylogenomics/scripts/get_strepto_genomes_patric.sh', 'w') as f: 
    for i in range(len(db_id)):
        line = bash_lines[i]
        f.writelines(line + db_id[i] + '\n') #combines bash lines to retrieve information from patric and the saving directory in a file

#add database_identifier as first column in strepto_genomes_quality.tsv'        
lines =[]
with open('/home/meike/strepto_phylogenomics/files/strepto_genomes_quality.tsv') as f:
    for line in f:
        line = line.strip().split('\t')
        lines.append(line)
with open('/home/meike/strepto_phylogenomics/files/streptococcus_genomes_quality.tsv', 'w') as f:
    for i, line in enumerate(lines):
        if i==0:
            f.writelines('\t'.join(["database_id"] + line) + '\n')
        else:
            f.writelines('\t'.join([strepto_db_id[i-1]] + line) + '\n')