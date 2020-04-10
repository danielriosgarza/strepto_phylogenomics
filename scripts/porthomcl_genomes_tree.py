#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 12:16:40 2020

@author: meike
"""

'''
Finding orthologs and paralogs using the software Porthomcl.

First get orthomcl: https://orthomcl.org/common/downloads/software/v2.0/UserGuide.txt
Get porthoMCL git clone and place it in orthomcSoftware/bin folder
$cp -a ~/PorthoMCL/. ~/orthomclSoftware-v2.0.9/bin/ --> to run porthomcl py scripts

Lines starting with '$' in this script were running in the terminal in the output dir of the porthomcl analysis.
'''

import os
from pathlib import Path
from datetime import date
from ete3 import Tree


def porthoMCL_prep(ids, savedir):
    '''
    Writes bash lines for PorthoMCL preparation. 
    '''
    with open (savedir, 'w') as f:
        for id_ in ids:
            f.write('orthomclAdjustFasta '+id_+' /home/meiker/git/data/prokka_annotation/'+id_+'/'+id_+'.faa 1\n')
        f.write('mv *.fasta /home/meiker/phylo_tree/orthomcl/compliantFasta')

def blast_run_bash(ids, savedir):
    '''
    writes bash line to run blast for all ids in taxon_list:
    blastp -query blastquery/DB_ID.fasta  -db blastdb/goodProteins.fasta  -seg yes  -dbsize 100000000  
    -evalue 1e-5  -outfmt 6 -num_threads 8 -out blastres/DB_ID.tab
    '''
    homepath = '/home/meiker/phylo_tree/orthomcl/'
    with open (savedir, 'w') as f:
        for id_ in ids:
            f.write("blastp -query " + homepath + "blastquery/"+id_[0]+".fasta  -db " + homepath + "blastdb/goodProteins.fasta  -seg yes  -dbsize 100000000  -evalue 1e-5  -outfmt 6 -num_threads 8 -out " + homepath + "blastres/"+id_[0]+".tab\n")
    
def blast_Parser_bash(ids, savedir):
    '''
    Writes bash for BlastParser for all ids in taxon_list
    porthomclBlastParser blastres/id_.tab compliantFasta >> splitSimSeq/id_.ss.tsv
    '''
    with open (savedir, 'w') as f:
        for id_ in ids:
            f.write("porthomclBlastParser /home/meiker/phylo_tree/orthomcl/blastres/" + id_[0] + ".tab /home/meiker/phylo_tree/orthomcl/compliantFasta >> /home/meiker/phylo_tree/orthomcl/splitSimSeq/" + id_[0] + ".ss.tsv\n")
            
def finding_best_hits(taxons, savedir):
    '''
     paralogs are found, and an unnormalized score is assigned to them. Step 5.3 will normalize 
     this score so that it be comparable among different genomes.
    '''
    with open(savedir, 'w') as f:
        for i in taxons:
            f.write('porthomclPairsBestHit.py -t /home/meiker/phylo_tree/orthomcl/taxon_list -s /home/meiker/phylo_tree/orthomcl/splitSimSeq -b /home/meiker/phylo_tree/orthomcl/besthit -q /home/meiker/phylo_tree/orthomcl/paralogTemp -x ' + str(i[1]) + ' -l /home/meiker/phylo_tree/orthomcl/logs/' + today + '_logfile_besthits.txt\n')
  
def split_files(ids, nsplits = 10):
    '''
    Splits db_list depending on length and makes lists with savdir that can be used for bash file generation.
    Tuples w/ (db_id, original_index)
    '''    
    groupsize = round(int(len(ids)/nsplits))
    
    db_index = [(id_, i + 1) for i, id_ in enumerate(ids)]
    
    id_i = [db_index[i :i + groupsize] for i in range(0, len(db_index), groupsize)]
    
    return id_i

def find_orthologs(taxons, savedir):
    '''
    Output of bash line is all the ortholog genes.
    '''
    with open (savedir, 'w') as f:
        for i in taxons:
            f.write("porthomclPairsOrthologs.py -t /home/meiker/phylo_tree/orthomcl/taxon_list -b /home/meiker/phylo_tree/orthomcl/besthit -o /home/meiker/phylo_tree/orthomcl/orthologs -x " + str(i[1]) + " -l /home/meiker/phylo_tree/orthomcl/logs/" + today + "_logfile_orthologs.txt\n")

def find_paralogs(taxons, savedir):
    '''
    Bash lines for finding paralogs. Uses split lists (db_id, index).
    '''
    with open (savedir, 'w') as f:
        for i in taxons:
            f.write ("porthomclPairsInParalogs.py -t /home/meiker/phylo_tree/orthomcl/taxon_list -q /home/meiker/phylo_tree/orthomcl/paralogTemp -o /home/meiker/phylo_tree/orthomcl/ogenes -p /home/meiker/phylo_tree/orthomcl/paralogs -x " + str(i[1]) +"\n")

    
    
path = os.getcwd()
p = Path(path)

#get the date to keep track of the scripts (added to scriptname)
today = date.today().strftime("%d/%m/%Y")
today = today.split('/')
today = ''.join(today)

#%% runcell 1
'''
Before adjusting fasta run:
$mkdir compliantFasta/   
    
1. Adjust fasta files:
    Get all db_ids and make bash scripts to adjust the fasta files
    Example bash line:
    orthomclAdjustFasta <id> /home/meiker/git/data/prokka_annotation/<id>/<id>.faa 1
    
Command creates output in the same folder it's ran
'''

#get ids that were used for the tree (203 genomes)
tree = Tree(os.path.join(p.parents[0], 'files', 'phylogenetic_tree', '12032020_reduced_concat_alignments.fa.contree'))
            
taxons = sorted(tree.get_leaf_names())

#adjust fasta files
porthoMCL_prep(taxons, os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl', today + '_porthomcl_prep.sh'))


'''
Get taxon list run following lines in the terminal (output dir): 
$ls -1 compliantFasta/ | sed -e 's/\..*$//'  > taxon_list

#place taxon list in git dir to get all ids in a single list
$cp -a ~/orthomcl/taxon_list ~/git/strepto_phylogenomics/files/
'''

#%% runcell 2
'''
2. Filter the input:
$orthomclFilterFasta compliantFasta/ 10 20
$mkdir filteredFasta
$mv goodProteins.fasta filteredFasta/
$mv poorProteins.fasta filteredFasta/

3.1 Create Blast DB:
$makeblastdb -in filteredFasta/goodProteins.fasta  -dbtype prot
$mkdir blastdb
$mv filteredFasta/goodProteins.fasta.* blastdb/

3.2 Split the inputfile:
$mkdir blastquery
$porthomclSplitFasta.py -i filteredFasta/goodProteins.fasta  -o blastquery
$mkdir blastres


3.3 Run Blast:
    Split the ids and generate several bash files to run blast on cluster
    Example bash line:
    blastp -query blastquery/<id>.fasta  -db blastdb/goodProteins.fasta  -seg yes  -dbsize 100000000  -evalue 1e-5  -outfmt 6 -num_threads 8 -out blastres/<id>.tab
'''

#split ids for parallel running of individual steps
splitted_ids = split_files(taxons)

for i, l_ids in enumerate(splitted_ids):
    blast_run_bash(l_ids, os.path.join(p.parents[0], 'scripts', 'bash_scripts' , 'porthomcl', today + '_blastrun' + str(i) + '.sh'))



#%% runcell 3
'''
4. Parse Blast results:
$mkdir splitSimSeq


Example bash line:
porthomclBlastParser blastres/<id>.tab compliantFasta >> splitSimSeq/<id>.ss.tsv
'''

splitted_ids_4 = split_files(taxons, nsplits = 4)

for i, l_ids in enumerate(splitted_ids_4):
    i += 1
    blast_Parser_bash(l_ids, os.path.join(p.parents[0], 'scripts', 'bash_scripts' , 'porthomcl', 'blastparser', today + '_blastparser' + str(i) + '.sh'))
    


#%% runcell 4
'''
5. Finding best hits:
$mkdir paralogTemp
$mkdir besthit
$mkdir logs #to save logfiles

Make bash script to find best hits (-x <number>, index of taxon to work on)
Example bash line:
porthomclPairsBestHit.py -t taxon_list -s splitSimSeq -b /besthit -q /paralogTemp -x <1> -l /logs/besthit_logfile.txt
'''

for i, l_ids in enumerate(splitted_ids_4):
    i += 1
    finding_best_hits(l_ids, os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl', 'besthit', today + '_find_best_hits' + str(i) +'.sh'))

#%% runcell 5
'''
6. Find orthologs
$mkdir orthologs

Example bash line (again -x <number> = taxon):
porthomclPairsOrthologs.py -t taxon_list -b besthit -o orthologs -x <1>

7. Find paralogos
$mkdir ogenes
$awk -F'[|\t]' '{print $4 >> ("ogenes/"$3".og.tsv")}' orthologs/*.ort.tsv
$awk -F'[|\t]' '{print $2 >> ("ogenes/"$1".og.tsv")}' orthologs/*.ort.tsv

Example bash line (again -x <number> = taxon):
porthomclPairsInParalogs.py -t taxon_list -q paralogTemp -o ogenes -p paralogs -x <1>
'''

for i, l_ids_i in enumerate(splitted_ids):
    i += 1
    find_orthologs(l_ids_i, os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl', 'orthologs', today + '_orthologs' + str(i) +'.sh'))


find_paralogs(l_ids_i, os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl', 'paralogs', today + '_paralogs.sh') )


'''
8. Run MCL
$cat orthologs/*.tsv >> all.ort.tsv
$mcl all.ort.tsv  --abc -I 1.5 -t 4 -o all.ort.group

$cat paralogs/*.tsv >> all.par.tsv
$mcl all.par.tsv  --abc -I 1.5 -t 4 -o all.par.group

To convert MCL files to a binary table read the scripts 'binary_table_preparation.py' and 'sorting_binary_table.py'
'''
