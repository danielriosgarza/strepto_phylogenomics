#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 11:16:48 2019

@author: meike
"""

'''Gene annotation using prokka'''
#prokka db_id.fna --/home/meiker/git/prokka_annotation  --prefix db_id.fna --genus Streptococcus

#function to generate command to run prokka (after genome is linked to directory)
def identifier2bash (database_ids, genus):
    '''Needs list with identifiers and returns list with prokka line for all ids'''
    bash_lines = []
    for identifier in database_ids:
            bash_lines.append(''.join('prokka /home/meiker/git/genomes/'+identifier+'.fna --outdir /home/meiker/git/prokka_annotation --prefix ' +identifier+' --genus ' +genus))
    return bash_lines
#retrieve selfmade database identifier
def get_db_ids(file):
    db_ids =[]
    with open(file) as f:
        for line in f:
            line= line.strip().split('/')
            identifier = line[-1]
            db_ids.append(identifier[0:-4])
    return db_ids

########## Floricoccus link to directory and running prokka ##########
flori_ids = get_db_ids('/home/meike/strepto_phylogenomics/scripts/get_floricoccus_genomes_patric.sh')
#get bash lines for running prokka
floricoccus_prokka = identifier2bash(flori_ids, 'Floricoccus') 

with open('/home/meike/strepto_phylogenomics/scripts/floricoccus_annotation_prokka.sh', 'w') as f:
    for identifier in floricoccus_prokka:
        f.write(identifier+ '\n')

########## Lactococcus annotation #########
lacto_ids = get_db_ids('/home/meike/strepto_phylogenomics/scripts/get_lactococcus_genomes_patric.sh')
lactococcus_prokka = identifier2bash(lacto_ids, 'Lactococcus')
with open ('/home/meike/strepto_phylogenomics/scripts/lactococcus_annotation_prokka.sh', 'w') as f:
    for identifier in lactococcus_prokka:
        f.write(identifier+ '\n')

########## Streptococcus annotation #########
strepto_ids = get_db_ids('/home/meike/strepto_phylogenomics/scripts/get_strepto_genomes_patric.sh')
streptococcus_prokka = identifier2bash(strepto_ids, 'Streptococcus')
with open ('/home/meike/strepto_phylogenomics/scripts/streptococcus_annotation_prokka.sh', 'w') as f:
   for identifier in streptococcus_prokka:
        f.write(identifier+ '\n')



