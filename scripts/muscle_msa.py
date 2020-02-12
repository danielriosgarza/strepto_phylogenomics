#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 10:56:16 2020

@author: meike
"""

'''
Bash script for MSA with Muscle (output Fasta format)
'''
import os
from pathlib import Path

path = os.getcwd()
p = Path(path)

genes = ["ADK", 'Exonuc_VII_L', 'RBFA']

with open (os.path.join(p,'bash_scripts', 'phylogenetic_tree', 'muscle_msa.sh'), 'w') as f:
    for gene in genes:
        command = '/home/meiker/software/muscle3.8.31_i86linux64 -in /home/meiker/phylo_tree/mfa/' + gene + ' -out /home/meiker/phylo_tree/msa_test/' + gene + ' -maxiters 1 -diags1 -sv\n'
        f.write(command)