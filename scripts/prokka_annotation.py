#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 11:16:48 2019

@author: meike
"""

'''Gene annotation using prokka'''
#prokka db_id.fna --/home/meiker/git/prokka_annotation  --prefix db_id.fna --genus Streptococcus

def make_prokka_annotation(genus):
    '''Give genus as string and writes bash lines to run prokka annotation. 
    Splits bash lines over 5 different files with each 2400 lines (last one all more than 9000'''
    genus= genus.lower()
    #count amount of genomes that needs to be annotated
    lines = []
    with open('/home/meike/strepto_phylogenomics/files/' + genus + '_genomes_quality.tsv') as f:
        for i, line in enumerate(f):
            lines.append(line)
    #if less than 2400 than generate one file with bash lines
    if len(lines) <= 2400:
        with open('/home/meike/strepto_phylogenomics/scripts/bash_scripts/' + genus + '_prokka_annotate_1.sh','w') as f1:
            for i, line in enumerate(lines):
                if i <= 2400:
                    a=line.strip().split('\t')
                    if not a[0].startswith('database_id'):
                        f1.write('prokka /home/meiker/git/genomes/'+a[0]+'.fna --outdir /home/meiker/git/prokka_annotation/' + a[0] + ' --prefix ' +a[0]+' --genus ' +genus +' \n')
    #else create 5 files with each aroung 2400 genomes
    else:    
        with open('/home/meike/strepto_phylogenomics/scripts/bash_scripts/' + genus + '_prokka_annotate_1.sh','w') as f1:
            with open('/home/meike/strepto_phylogenomics/scripts/bash_scripts/' + genus + '_prokka_annotate_2.sh','w') as f2:
                with open('/home/meike/strepto_phylogenomics/scripts/bash_scripts/' + genus + '_prokka_annotate_3.sh','w') as f3:
                    with open('/home/meike/strepto_phylogenomics/scripts/bash_scripts/' + genus + '_prokka_annotate_4.sh','w') as f4:
                        with open('/home/meike/strepto_phylogenomics/scripts/bash_scripts/' + genus + '_prokka_annotate_5.sh','w') as f5:
                            for i, line in enumerate(lines):
                                if i <= 2400:
                                    a=line.strip().split('\t')
                                    if not a[0].startswith('database_id'):
                                        f1.write('prokka /home/meiker/git/genomes/'+a[0]+'.fna --outdir /home/meiker/git/prokka_annotation/' + a[0] + ' --prefix ' +a[0]+' --genus ' +genus +' \n')
                                elif i >=2400 and i<= 4800:
                                    a=line.strip().split('\t')
                                    f2.write('prokka /home/meike/genomes/'+a[0]+'.fna --outdir /home/meiker/git/prokka_annotation/' + a[0] + ' --prefix ' +a[0]+' --genus ' +genus +' \n')
                                elif i >=4800 and i<= 7200:
                                    a=line.strip().split('\t')
                                    f3.write('prokka /home/meiker/git/genomes/'+a[0]+'.fna --outdir /home/meiker/git/prokka_annotation/' + a[0] + ' --prefix ' +a[0]+' --genus ' +genus +' \n')
                                elif i >=7200 and i<= 9600:
                                    a=line.strip().split('\t')
                                    f4.write('prokka /home/meiker/git/genomes/'+a[0]+'.fna --outdir /home/meiker/git/prokka_annotation/' + a[0] + ' --prefix ' +a[0]+' --genus ' +genus +' \n')
                                elif i >=9600:
                                    a=line.strip().split('\t')
                                    f5.write('prokka /home/meike/genomes/'+a[0]+'.fna --outdir /home/meiker/git/prokka_annotation/' + a[0] + ' --prefix ' +a[0]+' --genus ' +genus +' \n')
    return genus

make_prokka_annotation("streptococcus")
make_prokka_annotation('lactococcus')
make_prokka_annotation('Floricoccus')


##function to generate command to run prokka (after genome is linked to directory)
#def identifier2bash (database_ids, genus):
#    '''Needs list with identifiers and returns list with prokka line for all ids'''
#    bash_lines = []
#    for identifier in database_ids:
#            bash_lines.append(''.join('prokka /home/meiker/git/genomes/'+identifier+'.fna --outdir /home/meiker/git/prokka_annotation --prefix ' +identifier+' --genus ' +genus))
#    return bash_lines
##retrieve selfmade database identifier
#def get_db_ids(file):
#    db_ids =[]
#    with open(file) as f:
#        for line in f:
#            line= line.strip().split('/')
#            identifier = line[-1]
#            db_ids.append(identifier[0:-4])
#    return db_ids
#
########### Floricoccus link to directory and running prokka ##########
#flori_ids = get_db_ids('/home/meike/strepto_phylogenomics/scripts/get_floricoccus_genomes_patric.sh')
##get bash lines for running prokka
#floricoccus_prokka = identifier2bash(flori_ids, 'Floricoccus') 
#
#with open('/home/meike/strepto_phylogenomics/scripts/floricoccus_annotation_prokka.sh', 'w') as f:
#    for identifier in floricoccus_prokka:
#        f.write(identifier+ '\n')
#
########### Lactococcus annotation #########
#lacto_ids = get_db_ids('/home/meike/strepto_phylogenomics/scripts/get_lactococcus_genomes_patric.sh')
#lactococcus_prokka = identifier2bash(lacto_ids, 'Lactococcus')
#with open ('/home/meike/strepto_phylogenomics/scripts/lactococcus_annotation_prokka.sh', 'w') as f:
#    for identifier in lactococcus_prokka:
#        f.write(identifier+ '\n')
#
########### Streptococcus annotation #########
#strepto_ids = get_db_ids('/home/meike/strepto_phylogenomics/scripts/get_strepto_genomes_patric.sh')
#streptococcus_prokka = identifier2bash(strepto_ids, 'Streptococcus')
#with open ('/home/meike/strepto_phylogenomics/scripts/streptococcus_annotation_prokka.sh', 'w') as f:
#   for identifier in streptococcus_prokka:
#        f.write(identifier+ '\n')



