#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 11:16:48 2019

@author: meike
"""

'''Gene annotation using prokka'''
#prokka db_id.fna --/home/meiker/git/prokka_annotation  --prefix db_id.fna --genus Streptococcus

#function to Link the metagenome assembly file into desired directory
def link2genome(database_ids):
    bash_lines = []
    for identifier in database_ids:
        bash_lines.append(''.join('ln -fs home/meiker/git/genomes/' +identifier))
    return bash_lines
    
#function to generate command to run prokka (after genome is linked to directory)
def identifier2bash (database_ids, genus):
    '''Needs list with identifiers and returns list with prokka line for all ids'''
    bash_lines = []
    for identifier in database_ids:
            bash_lines.append(''.join('prokka '+identifier+' --/home/meiker/git/prokka_annotation/' +identifier+' --prefix ' +identifier+' --genus ' +genus))
    return bash_lines
#retrieve selfmade database identifier
def get_db_ids(file):
    db_ids =[]
    with open(file) as f:
        for line in f:
            line= line.strip().split('/')
            db_ids.append(line[-1])
    return db_ids

########## Floricoccus link to directory and running prokka ##########
flori_ids = get_db_ids('/home/meike/strepto_phylogenomics/scripts/get_floricoccus_genomes_patric.sh')
#bash to link genomes to annotation dir
floricoccus_link = link2genome(flori_ids)   
#get bash lines for running prokka
floricoccus_prokka = identifier2bash(flori_ids, 'Floricoccus') 

with open('/home/meike/strepto_phylogenomics/scripts/floricoccus_annotation_prokka.sh', 'w') as f:
    for i in range(len(floricoccus_link)):
        f.writelines(floricoccus_link[i] + '\n' + floricoccus_prokka[i] + '\n')

########## Lactococcus annotation #########
lacto_ids = get_db_ids('/home/meike/strepto_phylogenomics/scripts/get_lactococcus_genomes_patric.sh')
lacto_link = link2genome(lacto_ids)
lactococcus_prokka = identifier2bash(lacto_ids, 'Lactococcus')
with open ('/home/meike/strepto_phylogenomics/scripts/lactococcus_annotation_prokka.sh', 'w') as f:
    for i in range(len(lactococcus_prokka)):
        f.writelines(lacto_link[i] + '\n' + lactococcus_prokka[i] + '\n')

########## Streptococcus annotation #########
strepto_ids = get_db_ids('/home/meike/strepto_phylogenomics/scripts/get_strepto_genomes_patric.sh')
strepto_link = link2genome(strepto_ids)
streptococcus_prokka = identifier2bash(strepto_ids, 'Streptococcus')
with open ('/home/meike/strepto_phylogenomics/scripts/streptococcus_annotation_prokka.sh', 'w') as f:
    for i in range(len(streptococcus_prokka)):
        f.writelines(strepto_link[i] + '\n' + streptococcus_prokka[i] + '\n')



