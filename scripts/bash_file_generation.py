#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 09:41:25 2019

@author: meike
"""

'''Makes bash file to retrieve genomes from Patric. Saving directory is octarine server'''
#wget ftp://ftp.patricbrc.org/genomes/1074061.3/1074061.3.fna -O /home/meiker/git/genomes/streptocuccus_000001.fna

def identifier2bash (patric_ids_file):
    '''Needs file with identifiers and returns list with wget for all identifiers'''
    bash_lines = []
    with open (patric_ids_file) as f:
        for identifier in f:
            identifier = identifier.strip()
            bash_lines.append(''.join('wget ftp://ftp.patricbrc.org/genomes/'+identifier+'/' +identifier+'.fna'))
    return bash_lines

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

db_id = [] #generate ids for own database later. Containing in total 6 digits
for i in range(1, len(bash_lines)+1):
    db_id.append(' -O /home/meiker/git/genomes/streptocuccus_' + "%05d" % i + '.fna')

with open ('/home/meike/strepto_phylogenomics/scripts/get_strepto_genomes_patric.sh', 'w') as f: 
    for i in range(len(db_id)):
        line = bash_lines[i]
        f.writelines(line + db_id[i] + '\n') #combines bash lines to retrieve information from patric and the saving directory in a file
