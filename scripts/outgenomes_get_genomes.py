#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 14:28:03 2019

@author: meike
"""

'''Outgenomes bash file to retrieve genomes from Patric. Saving directory is octarine server'''
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

########## Lactococcus ##########
#Make file containing all ids of interest
ids =[]
with open ('/home/meike/strepto_phylogenomics/files/lactococcus_genomes_quality.tsv') as f:
    f_reader = csv.reader(f, delimiter="\t")
    for line in f_reader:
        if not line[0].startswith('genome.genome_id'):
            ids.append(line[0])                 
with open ('/home/meike/strepto_phylogenomics/files/lactococcus_genomes_ids.txt', 'w') as f:
    for identifier in ids:
        f.write(identifier+'\n')
#Create bash line to retrieve genomes from patric server
bash_lines=identifier2bash('/home/meike/strepto_phylogenomics/files/lactococcus_genomes_ids.txt')

#Generate own database ids
db_id = []
lacto_db_id = []
for i in range(1, len(bash_lines)+1):
    db_id.append(' -O /home/meiker/git/genomes/lactococcus_' + "%05d" % i + '.fna')
    lacto_db_id.append("lactococcus_"+"%05d" % i)
#combines bash lines to retrieve information from patric and the saving directory in a file
with open ('/home/meike/strepto_phylogenomics/scripts/get_lactococcus_genomes_patric.sh', 'w') as f: 
    for i in range(len(db_id)):
        line = bash_lines[i]
        f.writelines(line + db_id[i] + '\n') 
        
#add database_identifier as first column in lactococcus_genomes_quality.tsv'   
lines =[]
with open('/home/meike/strepto_phylogenomics/files/lactococcus_genomes_quality.tsv') as f:
    for line in f:
        line = line.strip().split('\t')
        lines.append(line)
with open('/home/meike/strepto_phylogenomics/files/lactococcus_genomes_quality.tsv', 'w') as f:
    for i, line in enumerate(lines):
        if i==0:
            f.writelines('\t'.join(["database_id"] + line) + '\n')
        else:
            f.writelines('\t'.join([lacto_db_id[i-1]] + line) + '\n')

########## Floricoccus ##########

ids =[]
with open ('/home/meike/strepto_phylogenomics/files/floricoccus_all_genome_fields.tsv') as f:
    f_reader = csv.reader(f, delimiter="\t")
    for line in f_reader:
        if not line[0].startswith('genome.genome_id'):
            ids.append(line[0])                 
with open ('/home/meike/strepto_phylogenomics/files/floricoccus_genomes_ids.txt', 'w') as f:
    for identifier in ids:
        f.write(identifier+'\n')
#Create bash line to retrieve genomes from patric server
bash_lines=identifier2bash('/home/meike/strepto_phylogenomics/files/floricoccus_genomes_ids.txt')

#Generate own database ids
db_id = [] #generate ids for own database later. Containing in total 6 digits
flori_db_id = []
for i in range(1, len(bash_lines)+1):
    db_id.append(' -O /home/meiker/git/genomes/floricoccus_' + "%05d" % i + '.fna')
    flori_db_id.append("floricoccus_"+"%05d" % i)
#combines bash lines to retrieve information from patric and the saving directory in a file
with open ('/home/meike/strepto_phylogenomics/scripts/get_floricoccus_genomes_patric.sh', 'w') as f: 
    for i in range(len(db_id)):
        line = bash_lines[i]
        f.writelines(line + db_id[i] + '\n')
        
#add database_identifier as first column in floricoccus_genomes_quality.tsv'
lines =[]
with open('/home/meike/strepto_phylogenomics/files/floricoccus_all_genome_fields.tsv') as f:
    for line in f:
        line = line.strip().split('\t')
        lines.append(line)
with open('/home/meike/strepto_phylogenomics/files/floricoccus_genomes_quality.tsv', 'w') as f:
    for i, line in enumerate(lines):
        if i==0:
            f.writelines('\t'.join(["database_id"] + line) + '\n')
        else:
            f.writelines('\t'.join([flori_db_id[i-1]] + line) + '\n')



            
            