#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 09:56:44 2020

@author: meike
"""

'''
Building phylogenetic tree based on HMM profiles. (Python 3)
'''

import os
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

path = os.getcwd()
p = Path(path)

def find_best_hits(hmm_file_location):
    hits = {}
    evalues = {}
    with open (hmm_file_location) as f:
        for line in f:
            if line[0] == '#':
                pass
            else:
                split_line = line.split()
                target = split_line[0]
                query_name = split_line[2]
    #            print target
    #            print query_name
                if query_name not in hits:
                    hits[query_name] = target
                    evalues[query_name] = split_line[4]
    return hits, evalues

def insert_newlines(string, every=60):
    return '\n'.join(string[i:i+every] for i in range(0, len(string), every))

profiles = '/home/meiker/phylo_tree/genes.hmm'
output_path = '/home/meiker/phylo_tree/'
files_dir = '/home/meiker/git/data/prokka_annotation/'

#%% runcell 1

'''
Homology search for each proteome using 71 HMM profiles.
Use --tblout to save a parseble table which is used to determine the best hit. (output in hmm/)
Use -o to output a file from which the alignment from the protein sequence to the profile can be parsed.(output in hmmalign/)
'''

db_ids = []
with open (os.path.join(p.parents[0], 'files', 'taxon_list')) as f:
    for line in f:
        a = line.strip()
        db_ids.append(a)

#first path give outdir
with open (os.path.join(p.parents[0],'scripts', 'bash_scripts', 'pyhlogenetic_tree', 'get_hmm_search.sh'), 'w') as f:
    for id_ in db_ids:
        f.write('hmmsearch --tblout '+ output_path + 'output_hmm/hmm/'+id_+'.txt -o '+ output_path+'output_hmm/hmmalign/'+id_+' --cut_tc --cpu 8 ' + profiles + ' ' + files_dir + id_ + '/'+ id_ +'.faa\n\n')
                
#%% runcell 2

'''
Copy files to local drive ti check the hits. (change paths)
Find best hit for each gene, save target name.
Parse both output files to get a dict that maps gene name to protein sequence for every strain.
'''

#get all hits and e_values from each proteome
hits = {}
e_values = {}
for id_ in db_ids:
    hits[id_], e_values[id_] = find_best_hits(output_path+ 'hmm/' + id_ + '.txt')


profile_seqs = {id_:{} for id_ in hits}
for id_ in profile_seqs:
    prot_seq = ''
    protein_name = 'tempname'
    gene_map = {v:k for k,v in hits[id_].items()}
    with open (output_path +'/hmmalign/'+id_) as f:
        for line in f:
            if line[0] != '#':
                if line[0] == '>':
                    if protein_name in gene_map:
                        gene_name = gene_map[protein_name]
                        profile_seqs[id_][gene_name] = prot_seq #Store protein sequence to profile_seqs dict under gene key.
                    protein_name = line.strip().split()[1]
                    prot_seq = ''
                elif protein_name in line:
                    prot_seq += line.strip().split()[2].replace('-', '').upper() #Leave out gaps in sequence
                
#Throw out bad hits
bad_hits = []
for id_ in list(profile_seqs):
    for gene in list(profile_seqs[id_]):
        if len(profile_seqs[id_][gene]) < 10:
            del profile_seqs[id_][gene]
            bad_hits.append((id_, gene))
#all_genes_found = [strain for strain in hits.keys() if len(hits[strain]) == nr_of_genes]
nr_genes_found = {id_:len(hits[id_]) for id_ in hits.keys()}

nr_genes_found_vals = [i for i in nr_genes_found.values()]

genes = {}
for id_ in hits:
    for gene in hits[id_]:
        if gene not in genes:
            genes[gene] = {id_}
        else:
            genes[gene].add(id_)

nr_times_found = {gene:len(genes[gene]) for gene in genes.keys()}

nr_times_found_vals = [i for i in nr_times_found.values()]

#Check the hits
nr_genome = np.asarray(nr_genes_found_vals)
plt.hist(nr_genome)
plt.savefig(output_path+'nr_genomes.png', bbox_inches='tight')
nr_gene = np.asarray(nr_times_found_vals)
plt.hist(nr_gene)
plt.savefig(output_path+'nr_genes.png', bbox_inches='tight')

#%% runcell 3

'''
Make multifasta file for each gene
'''

for gene in genes:
    with open (output_path+'mfa/'+gene, 'w') as f:
        for id_ in profile_seqs:
            if gene in profile_seqs[id_]:
                f.write('>'+id_ + '\n')
                sequence_w_newlines = insert_newlines(profile_seqs[id_][gene])
                f.write(sequence_w_newlines.strip()+'\n') #.strip() to aviod double \n\n.
