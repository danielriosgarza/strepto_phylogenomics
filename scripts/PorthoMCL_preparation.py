#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 12:10:37 2020

@author: meike
"""

'''
Run first the porthoMCL_prep for all species and follow than the steps in the PorthoPrep Steps document in the terminal,
before running the rest.

'''
#import random
from random import shuffle
import os
from pathlib import Path
from datetime import date

def get_ids(file):
    '''
    Gets ids from database-patric_id file and returns a list with all ids
    '''
    with open (file) as f:
        db_ids = []
        for line in f:
            a = line.strip().split('\t')
            if a[0] != 'database_id':
                db_ids.append(a[0])
        db_ids.sort()
    return db_ids

def porthoMCL_prep(db_ids, savedir):
    '''
    Writes bash lines for PorthoMCL preparation. 
    '''
    with open (savedir, 'w') as f:
        for id_ in db_ids:
            f.write('orthomclAdjustFasta '+id_+' /home/meiker/git/data/prokka_annotation/'+id_+'/'+id_+'.faa 1\n')
        f.write('mv *.fasta /home/meiker/orthomcl/compliantFasta')

def get_taxon_list(taxon_list):
    '''
    After creation of taxon_list in terminal, run function to get all ids.
    '''
    db_ids = []
    with open(taxon_list) as f:
        for line in f:
            line = line.strip()
            db_ids.append(line)
    return db_ids


def blast_run_bash(db_ids, savedir):
    '''
    writes bash line to run blast for all ids in taxon_list:
    blastp -query blastquery/DB_ID.fasta  -db blastdb/goodProteins.fasta  -seg yes  -dbsize 100000000  
    -evalue 1e-5  -outfmt 6 -num_threads 8 -out blastres/DB_ID.tab
    '''
    homepath = '/home/meiker/orthomcl/'
    with open (savedir, 'w') as f:
        for id_ in db_ids:
            f.write("blastp -query " + homepath + "blastquery/"+id_[0]+".fasta  -db " + homepath + "blastdb/goodProteins.fasta  -seg yes  -dbsize 100000000  -evalue 1e-5  -outfmt 6 -num_threads 8 -out " + homepath + "blastres/"+id_[0]+".tab\n")
    
def blast_Parser_bash(db_ids, savedir):
    '''
    Writes bash for BlastParser for all ids in taxon_list
    porthomclBlastParser blastres/id_.tab compliantFasta >> splitSimSeq/id_.ss.tsv
    '''
    with open (savedir, 'w') as f:
        for id_ in db_ids:
            f.write("porthomclBlastParser /home/meiker/orthomcl/blastres/"+id_+".tab /home/meiker/orthomcl/compliantFasta >> /home/meiker/orthomcl/splitSimSeq/"+id_+".ss.tsv\n")
            
def finding_best_hits(db_ids, savedir):
    '''
     paralogs are found, and an unnormalized score is assigned to them. Step 5.3 will normalize 
     this score so that it be comparable among different genomes.
    '''
    with open(savedir, 'w') as f:
        for i, id_ in enumerate(db_ids):
            f.write("porthomclPairsBestHit.py -t /home/meiker/orthomcl/taxon_list -s /home/meiker/orthomcl/splitSimSeq -b /home/meiker/orthomcl/besthit -q /home/meiker/orthomcl/paralogTemp -x "+str(i + 1)+"\n")
  
def split_files(db_ids):
    '''
    Splits db_list depending on length and makes lists with savdir that can be used for bash file generation.
    Tuples w/ (db_id, original_index)
    '''    
    #groupsize = int(len(db_ids)/16)
    groupsize = 761
    db_index = [(id_, i + 1) for i, id_ in enumerate(db_ids)]
    
    id_i = [db_index[i :i + groupsize] for i in range(0, len(db_index), groupsize)]
    
    return id_i

def find_orthologs(db_ids, savedir):
    '''
    Output of bash line is all the ortholog genes.
    '''
    with open (savedir, 'w') as f:
        for i, id_ in enumerate(db_ids):
            f.write("porthomclPairsOrthologs.py -t /home/meiker/orthomcl/taxon_list -b /home/meiker/orthomcl/besthit -o /home/meiker/orthomcl/orthologs -x "+ str(i + 1)+"\n")

def find_paralogs(db_ids, savedir):
    '''
    Bash lines for finding paralogs. Uses split lists (db_id, index).
    '''
    with open (savedir, 'w') as f:
        for i, id_ in enumerate(db_ids):
            f.write ("porthomclPairsInParalogs.py -t /home/meiker/orthomcl/taxon_list -q /home/meiker/orthomcl/paralogTemp -o /home/meiker/orthomcl/ogenes -p /home/meiker/orthomcl/paralogs -x "+str(i + 1)+"\n")

def randomizer(infile, outfile):
    '''
    Randomizes order of bash lines
    '''
    lines = []
    with open(infile) as f:
        for line in f:
            lines.append(line.strip())
    shuffle(lines)
    with open(outfile, 'w') as f2:
        for line in lines:
            f2.write(line + '\n')
    
    
    
path = os.getcwd()
p = Path(path)

#%% runcell 1
'''
First get orthomcl: https://orthomcl.org/common/downloads/software/v2.0/UserGuide.txt
Get porthoMCL git clone and place it in orthomcSoftware/bin folder
$cp -a ~/PorthoMCL/. ~/orthomclSoftware-v2.0.9/bin/ --> to run porthomcl py scripts

in the following lines '$' terminal commands were running in the output dir of porthomcl analysis

Before adjusting fasta run:
$mkdir compliantFasta/   
    
1. Adjust fasta files:
    Get all db_ids and make bash scripts to adjust the fasta files
    Example bash line:
    orthomclAdjustFasta <id> /home/meiker/git/data/prokka_annotation/<id>/<id>.faa 1
'''

# #get ids
# strepto_ids = get_ids(os.path.join(p.parents[0], 'files', '20012020streptococcus_patric_id_with_database_id.tsv'))
# flori_ids = get_ids(os.path.join(p.parents[0], 'files', 'floricoccus_patric_id_with_database_id.tsv'))
# lacto_ids = get_ids(os.path.join(p.parents[0], 'files', 'lactococcus_patric_id_with_database_id.tsv'))

# #adjust fasta files
# porthoMCL_prep(strepto_ids, os.path.join(p.parents[0], 'scripts', 'bash_scripts', '20012020_streptococcus_porthomcl_prep.sh'))
# porthoMCL_prep(flori_ids, os.path.join(p.parents[0], 'scripts', 'bash_scripts', '20012020_floricoccus_porthomcl_prep.sh'))
# porthoMCL_prep(lacto_ids, os.path.join(p.parents[0], 'scripts', 'bash_scripts', '20012020_lactococcus_porthomcl_prep.sh'))

'''
Get taxon list run following lines in the terminal (putput dir): 
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


3.3 Run Blast:
    Split the ids and generate several bash files to run blast on cluster
    Example bash line:
    blastp -query blastquery/<id>.fasta  -db blastdb/goodProteins.fasta  -seg yes  -dbsize 100000000  -evalue 1e-5  -outfmt 6 -num_threads 8 -out blastres/<id>.tab
'''

#get all ids in a single list and split it for blast run
# db_ids = get_taxon_list(os.path.join(p.parents[0], 'files', 'taxon_list'))
# splitted_ids = split_files(db_ids)

# for i, l_ids in enumerate(splitted_ids):
#     blast_run_bash(l_ids, os.path.join(p.parents[0], 'scripts', 'bash_scripts' , 'porthomcl', '200120_blastrun'+str(i)+'.sh'))

#Blastrun takes longer than my internship period, therefore, the streptococcus genmoes will be randomized
# for i in range(1, 16):
#     randomizer(os.path.join(p.parents[0], 'scripts', 'bash_scripts' , 'porthomcl', '200120_blastrun'+str(i)+'.sh'), os.path.join(p.parents[0], 'scripts', 'bash_scripts' , 'porthomcl', '210120_blastrun'+str(i)+'.sh'))

#%% runcell 3
'''
4. Parse Blast results:
$mkdir splitSimSeq

Make bash files that can be updated according to already blasted files
Example bash line:
porthomclBlastParser blastres/<id>.tab compliantFasta >> splitSimSeq/<id>.ss.tsv
'''

#get the date to keep track of the scripts (added to scriptname)
today = date.today().strftime("%d/%m/%Y")
today = today.split('/')
today = ''.join(today)

#Get ids_ of all genomes that are already blasted (located in blastres/)
dbs_done = []
for file in list(os.listdir('/home/meiker/orthomcl/blastres')):    
    id_ = file.strip().split('.')[0]
    dbs_done.append(id_)

print(len(dbs_done))
blast_Parser_bash(dbs_done, os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl', today+'_blastparser.sh'))

#%% runcell 4
'''
5. Finding best hits:
$mkdir paralogTemp
$mkdir besthit

Make bash script to finf best hits (-x <number>, index of taxon to work on)
Example bash line:
porthomclPairsBestHit.py -t taxon_list -s splitSimSeq -b besthit -q paralogTemp -x <1>
'''

finding_best_hits(dbs_done, os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl', today + '_find_best_hits.sh'))

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

find_orthologs(dbs_done, os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl', today + '_orthologs.sh'))

find_paralogs(dbs_done, os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl', today+ '_paralogs.sh'))

'''
8. Run MCL
$cat orthologs/*.tsv >> all.ort.tsv
$mcl all.ort.tsv  --abc -I 1.5 -t 4 -o .all.ort.group

$cat paralogs/*.tsv >> all.par.tsv
$mcl all.par.tsv  --abc -I 1.5 -t 4 -o all.par.group
'''

# test_ids = get_taxon_list(os.path.join(p.parents[0], "files", 'porthomcl', 'taxon_list'))
# ids_split = split_files(test_ids)

# for i, l_ids in enumerate(ids_split):
#     blast_run_bash(l_ids, os.path.join(p.parents[0], 'scripts', 'bash_scripts' , 'porthomcl', 'testset_blastrun'+str(i)+'.sh'))

# blast_Parser_bash(test_ids, os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl','testset_blastparser.sh'))
# finding_best_hits(test_ids, os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl','testset_besthits.sh'))
# find_orthologs(test_ids, os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl', 'testset_orthologs.sh'))
# find_paralogs(test_ids, os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl', 'testset_paralogs.sh'))    