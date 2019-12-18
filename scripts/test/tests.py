#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 09:16:34 2019

@author: meike
"""

#--outdir /home/meiker/git/prokka_annotation/streptococcus_02451
# 'streptococcus_11895', 'streptococcus_11939']
# 'lactococcus_00188'
#wget ftp://ftp.patricbrc.org/genomes/1739284.3/1739284.3.fna -O /home/meiker/git/genomes/streptocuccus_11846.fna

missing_genomes= ['streptocuccus_11895.fna', 'streptocuccus_11939.fna', 'lactococcus_00188.fna']

with open ('/home/meike/strepto_phylogenomics/scripts/get_strepto_genomes_patric.sh') as f1:
    with open ('/home/meike/strepto_phylogenomics/scripts/get_lactococcus_genomes_patric.sh') as f2:
        with open ('/home/meike/strepto_phylogenomics/scripts/get_genomes_missed.sh', 'w') as f3:
            for line in f1:
                a = line.strip().split('/')
                if a[-1] in missing_genomes:
                    f3.write(line)
            for line in f2:
                b = line.strip().split('/')
                if b[-1] in missing_genomes:
                    f3.write(line)