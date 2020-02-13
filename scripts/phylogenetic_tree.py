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


def check_alignment_len(path_to_file):
    with open (path_to_file) as f:
        lens = []
        counter = 0
        for line in f:
            if '>' in line:
                if counter != 0:
                    lens.append(counter)
                counter = 0
            else:
                counter += len(line.strip())
    return lens

profiles = '/home/meiker/phylo_tree/genes.hmm'
phylo_path = '/home/meiker/phylo_tree/'
files_dir = '/home/meiker/git/data/prokka_annotation/'
bashscript_path = '/home/meiker/git/strepto_phylogenomics/scripts/bash_scripts/'

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

# #first path give outdir
# with open (os.path.join(p.parents[0],'scripts', 'bash_scripts', 'pyhlogenetic_tree', 'get_hmm_search.sh'), 'w') as f:
#     for id_ in db_ids:
#         f.write('hmmsearch --tblout '+ phylo_path + 'output_hmm/hmm/'+id_+'.txt -o '+ phylo_path+'output_hmm/hmmalign/'+id_+' --cut_tc --cpu 8 ' + profiles + ' ' + files_dir + id_ + '/'+ id_ +'.faa\n\n')
                
#%% runcell 2

'''
Find best hit for each gene, save target name.
Parse both output files to get a dict that maps gene name to protein sequence for every strain.
'''

#get all hits and e_values from each proteome
# hits = {}
# e_values = {}
# for id_ in db_ids:
#     hits[id_], e_values[id_] = find_best_hits(phylo_path+ 'output_hmm/hmm/' + id_ + '.txt')


# profile_seqs = {id_:{} for id_ in hits}
# for id_ in profile_seqs:
#     prot_seq = ''
#     protein_name = 'tempname'
#     gene_map = {v:k for k,v in hits[id_].items()}
#     with open (phylo_path +'output_hmm/hmmalign/'+id_) as f:
#         for line in f:
#             if line[0] != '#':
#                 if line[0] == '>':
#                     if protein_name in gene_map:
#                         gene_name = gene_map[protein_name]
#                         profile_seqs[id_][gene_name] = prot_seq #Store protein sequence to profile_seqs dict under gene key.
#                     protein_name = line.strip().split()[1]
#                     prot_seq = ''
#                 elif protein_name in line:
#                     prot_seq += line.strip().split()[2].replace('-', '').upper() #Leave out gaps in sequence
                
# #Throw out bad hits
# bad_hits = []
# for id_ in list(profile_seqs):
#     for gene in list(profile_seqs[id_]):
#         if len(profile_seqs[id_][gene]) < 10:
#             del profile_seqs[id_][gene]
#             bad_hits.append((id_, gene))
# #all_genes_found = [strain for strain in hits.keys() if len(hits[strain]) == nr_of_genes]
# nr_genes_found = {id_:len(hits[id_]) for id_ in hits.keys()}

# nr_genes_found_vals = [i for i in nr_genes_found.values()]

# genes = {}
# for id_ in hits:
#     for gene in hits[id_]:
#         if gene not in genes:
#             genes[gene] = {id_}
#         else:
#             genes[gene].add(id_)

# nr_times_found = {gene:len(genes[gene]) for gene in genes.keys()}

# nr_times_found_vals = [i for i in nr_times_found.values()]

# #Check the hits
# nr_genome = np.asarray(nr_genes_found_vals)
# plt.hist(nr_genome)
# plt.savefig(phylo_path+'nr_genomes.png', bbox_inches='tight')
# nr_gene = np.asarray(nr_times_found_vals)
# plt.hist(nr_gene)
# plt.savefig(phylo_path+'nr_genes.png', bbox_inches='tight')

#%% runcell 3

# '''
# Make multifasta file for each gene
# '''

# for gene in genes:
#     with open (phylo_path+'mfa/'+gene, 'w') as f:
#         for id_ in profile_seqs:
#             if gene in profile_seqs[id_]:
#                 f.write('>'+id_ + '\n')
#                 sequence_w_newlines = insert_newlines(profile_seqs[id_][gene])
#                 f.write(sequence_w_newlines.strip()+'\n') #.strip() to aviod double \n\n.
                
#%% runcell 4

'''
Run a multiple sequence aligner over each mfa file.
Run "run_msa.sh"
'''

# if not os.path.isdir(phylo_path+'msa'):
#     os.mkdir(phylo_path+'msa')
    
#To large dataset
# with open (os.path.join(p,'bash_scripts', 'phylogenetic_tree', 'run_msa.sh'), 'w') as f:
#     for gene in genes:
#         command = '/home/meiker/software/clustalo-1.2.4-Ubuntu-x86_64 -i ' + phylo_path + 'mfa/' + gene + ' -o ' + phylo_path + 'msa/' + gene + ' --dealign --threads 8\n'
#         # f.write(command)

#Muscle as MSA –clwstrict: writes output in ClustalW format
# with open (os.path.join(p,'bash_scripts', 'phylogenetic_tree', 'muscle_msa.sh'), 'w') as f:
#     for gene in genes:
#         command = '/home/meiker/software/muscle3.8.31_i86linux64 -in ' + phylo_path + 'mfa/' + gene + ' -out ' + phylo_path + 'msa/' + gene + ' -maxiters 1 -diags1 -sv -clwstrict\n'
#         f.write(command)

# with open (bashscript_path + 'phylogenetic_tree/muscle_msa.sh', 'w') as f:
#     for gene in genes:
#         command = '/home/meiker/software/muscle3.8.31_i86linux64 -in /home/meiker/phylo_tree/mfa/' + gene + ' -out /home/meiker/phylo_tree/msa/' + gene + ' -maxiters 1 -diags1 -sv\n'
#         f.write(command)

# msa_lens = {}
# for gene in genes:
#     msa_lens[gene] = let(check_alignment_len(phylo_path+'msa_trimmed/'+gene))
    
#%% runcell 5

'''
Write a bash script that trims every gene msa file.
'''

# with open (bashscript_path + 'phylogenetic_tree/run_trimal.sh', 'w') as f:
#     for gene in genes:
#         f.write('trimal -in '+ phylo_path+ 'msa/' + gene+' -out '+phylo_path+'msa_trimmed/'+gene+' -automated1\n')



#%% runcell 6

'''
Concetenate all files together that can be used for distance analysis with iqtree (fasta format)
1. 
a) Check each alignment file (per gene) if all genomes (ids) are present. 
b) If not add the database_id with a gap sequence ('-') with same length as the other alignments.
2. Write concatenated file: '>db_id' and next line all alignments "glued" to each other
(example: gene1gene2---3gen-4g--e5...)
'''
path= phylo_path+'msa_trimmed/'
files = os.listdir(path)

#dict to store all alignments: key=filename(gene), value: dict with key: db_id and value: seq alignment of that gene
all_alignments ={} 

for fl in files:
    sequences ={} #temp dict that is asigned to all_alignments later 
    aligned_ids = [] #to see which ids (genomes) had the gene
    with open(path+fl) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                seq=''
                db_id = line[1::] #get rid of '>'
                if db_id not in aligned_ids:
                    aligned_ids.append(db_id)
                    sequences[db_id]=''
            else:
                seq+=line
            sequences[db_id]=seq
    #Check which genomes (db_ids) miss certain gene and add gaps ('-') for this db_id
    for id_ in db_ids: 
        if id_ not in aligned_ids:
            sequences[id_]=''
    
    #Determine the longest alignment and add gaps to alignment if they are shorter (or miss the gene)       
    longest = len(max(sequences.values())) 
    
    for k, v in sequences.items():
        if len(v) < longest:
            difference = longest-len(v)
            gap = '-'*difference
            sequences[k] = v+gap
    all_alignments[fl] = sequences #add to all alignments to store the alignments per gene     

#make dict with each id (genome) as key and all alignments of each gene glued together as values
alignments_per_id = {}    
for gene in all_alignments:
        seqD = all_alignments[gene] #get dict with alignments per gene
        for id_ in seqD:
            if id_ not in alignments_per_id:
                alignments_per_id[id_] = seqD[id_]
            else:
                alignments_per_id[id_]+=seqD[id_]
                
#write concatenated file          
if not os.path.isdir(phylo_path+'iqtree'):
    os.mkdir(phylo_path+'iqtree')
                
with open(phylo_path+'/iqtree/concat_alignments', 'w') as f:
    for id_, seqs in alignments_per_id.items():
        f.write('>'+id_+'\n'+seqs+'\n')

'''
Run iqtree wtih following command line:
iqtree -s <concat_file> -bb 1000 -alrt 1000 -nt AUTO -ntmax 50

iqtree -s '/home/meiker/phylo_tree/iqtree/concat_alignments -bb 1000 -alrt 1000 -nt AUTO
-alert --> specifies the number of bootstrap replicates for SH-aLRT (1000 is minimum number recommended)
-bb --> number of bootstrap replicates (1000 is minimum number recommended)
-ntmax 8 --> determine max cores that might be used (otherwise all will be used)
-nt AUTO --> determines best number of cores
'''
#%%
# '''
# Write one msa file where the seqs of all genes are concatenated per genome (id).
# '''
# #Test if each genome has a sequence for enough genes to be considered.
# usable_ids = [id_ for id_ in hits if len(hits[id_].values()) > 10]

# for id_ in hits:
#     len_hits = [i for i in hits[id_].values()]
#     if len(hits[id_].values()) < 10:
#         print(id_)
        
# ids_found_per_gene = {gene:[] for gene in genes}
# for gene in genes:
#     with open (phylo_path+'msa_trimmed/'+gene) as f:
#         for line in f:
#             if '>' in line:
#                 ids_found_per_gene[gene].append(line.strip().split()[0][1:])

# '''
# Add gaps for genes that were not found.
# '''

# for gene in genes:
#     gene_len = get_alignment_len(phylo_path+'msa_trimmed/'+gene)
#     gap = insert_newlines('-'*gene_len)
#     for id_ in usable_ids:
#         if id_ not in ids_found_per_gene[gene]:
#             with open(phylo_path+'msa_trimmed/'+gene, 'a') as f:
#                 f.write('>'+id_+'\n')
#                 f.write(gap+'\n')
                
                
# concatenated_seqs = {id_:'' for id_ in usable_ids}
# concatenation_len = {id_:0 for id_ in usable_ids}
# for gene in genes:
#     with open(phylo_path+'msa_trimmed/'+gene) as f:
#         for line in f:
#             if '>' in line:
#                 id_ = line.strip().split()[0][1:]
#             else:
#                 if id_ in concatenated_seqs:
#                     concatenated_seqs[id_] += line.strip()
#                     concatenation_len[id_] += len(line.strip())

# with open(phylo_path+'msa_trimmed/concatenated_msa_trimmed', 'w') as f:
#     for id_ in concatenated_seqs:
#         f.write('>'+id_+'\n')
#         seq = insert_newlines(concatenated_seqs[id_]).strip()
#         f.write(seq+'\n')

# '''
# Run iqtree on concatenated_msa_trimmed.
# With iqtree -s concatenated_msa_trimmed -bb 1000 -alrt 1000 -nt AUTO
# '''
    
    