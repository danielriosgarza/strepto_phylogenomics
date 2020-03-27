#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 12:10:37 2020

@author: meike
"""

'''
Finding orthologs and paralogs using the software Porthomcl.
'''
import random
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

def porthoMCL_prep(ids, savedir):
    '''
    Writes bash lines for PorthoMCL preparation. 
    '''
    with open (savedir, 'w') as f:
        for id_ in ids:
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


def blast_run_bash(ids, savedir):
    '''
    writes bash line to run blast for all ids in taxon_list:
    blastp -query blastquery/DB_ID.fasta  -db blastdb/goodProteins.fasta  -seg yes  -dbsize 100000000  
    -evalue 1e-5  -outfmt 6 -num_threads 8 -out blastres/DB_ID.tab
    '''
    homepath = '/home/meiker/orthomcl/'
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
            f.write("porthomclBlastParser /home/meiker/orthomcl/blastres/"+id_[0]+".tab /home/meiker/orthomcl/compliantFasta >> /home/meiker/orthomcl/splitSimSeq/"+id_[0]+".ss.tsv\n")
            
def finding_best_hits(indexes_ids, savedir):
    '''
     paralogs are found, and an unnormalized score is assigned to them. Step 5.3 will normalize 
     this score so that it be comparable among different genomes.
    '''
    with open(savedir, 'w') as f:
        for i in indexes_ids:
            f.write("porthomclPairsBestHit.py -t /home/meiker/orthomcl/taxon_list -s /home/meiker/orthomcl/splitSimSeq -b /home/meiker/orthomcl/besthit -q /home/meiker/orthomcl/paralogTemp -x "+str(i[0])+"\n")
  
def split_files(ids, nsplits = 16):
    '''
    Splits db_list depending on length and makes lists with savdir that can be used for bash file generation.
    Tuples w/ (db_id, original_index)
    '''    
    groupsize = round(int(len(ids)/nsplits))
    
    db_index = [(id_, i + 1) for i, id_ in enumerate(ids)]
    
    id_i = [db_index[i :i + groupsize] for i in range(0, len(db_index), groupsize)]
    
    return id_i

def find_orthologs(indexes_ids, savedir):
    '''
    Output of bash line is all the ortholog genes.
    '''
    with open (savedir, 'w') as f:
        for i in indexes_ids:
            f.write("porthomclPairsOrthologs.py -t /home/meiker/orthomcl/taxon_list -b /home/meiker/orthomcl/besthit -o /home/meiker/orthomcl/orthologs -x "+ str(i)+"\n")

def find_paralogs(indexes_ids, savedir):
    '''
    Bash lines for finding paralogs. Uses split lists (db_id, index).
    '''
    with open (savedir, 'w') as f:
        for i in indexes_ids:
            f.write ("porthomclPairsInParalogs.py -t /home/meiker/orthomcl/taxon_list -q /home/meiker/orthomcl/paralogTemp -o /home/meiker/orthomcl/ogenes -p /home/meiker/orthomcl/paralogs -x "+str(i)+"\n")

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
#For the Updater: Blast parser only parses results that are added to blastres folder (finished blasted), the rest of the orthologue search and MCL uses all indexes of the parsed blast results. If it gets updated, it needs to run the scripts also on the earlier files!(because of the newly generated tsv files that migth contain additional orthologues etc.)

#get the date to keep track of the scripts (added to scriptname)
today = date.today().strftime("%d/%m/%Y")
today = today.split('/')
today = ''.join(today)

#get all ids (genomes) that were already blasted and parsed
dbs_parsed = []
for file in list(os.listdir('/home/meiker/orthomcl/splitSimSeq/')):
    id_ = file.strip().split('.')[0]
    dbs_parsed.append(id_)

#Get ids_ of all genomes that are already blasted (located in blastres/) and not further analysed
dbs_ready2analyze = []
for file in list(os.listdir('/home/meiker/orthomcl/blastres')):    
    id_ = file.strip().split('.')[0]
    if id_ not in dbs_parsed:
        dbs_ready2analyze.append(id_)

#Get the index from the taxon list of the ready2analyse ids
taxon_list = []
with open(os.path.join(p.parents[0], 'files', 'taxon_list')) as f:
    for line in f:
        taxon_list.append(line.strip())

indexes = []        
for i, taxon in enumerate(taxon_list):
    if taxon in dbs_parsed:
        indexes.append(i+1)
        
#to give an indication how many are added
print(len(dbs_ready2analyze))

split_ids = split_files(dbs_ready2analyze, nsplits = 12)

for i, l_ids in enumerate(split_ids):
    i += 1
    blast_Parser_bash(l_ids, os.path.join(p.parents[0], 'scripts', 'bash_scripts' , 'porthomcl', 'blastparser', today + '_blastparser'+str(i)+'.sh'))
    


#%% runcell 4
'''
5. Finding best hits:
$mkdir paralogTemp
$mkdir besthit

Make bash script to find best hits (-x <number>, index of taxon to work on)
Example bash line:
porthomclPairsBestHit.py -t taxon_list -s splitSimSeq -b besthit -q paralogTemp -x <1>
'''
split_ids = split_files(indexes, nsplits = 12)

for i, l_inds in enumerate(split_ids):
    finding_best_hits(l_inds, os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl', 'besthit', + today + '_find_best_hits' + str(i) +'.sh'))

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
#Split the indexes into 10 sections (orthology search takes otherwise really long)

size = round(len(indexes)/10)

line = 'porthomclPairsOrthologs.py -t /home/meiker/orthomcl/taxon_list -b /home/meiker/orthomcl/besthit -o /home/meiker/orthomcl/orthologs -x '

for i in range(1, 11):
    with open(os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl', 'orthologs', today + '_orthologs' + str(i) + '.sh'), 'w') as f:
        if i == 1:
            numbers = indexes[:size]
            for n in numbers:
                f.write(line + str(n) + '\n')        
        if i == 2:
            numbers = indexes[size:size*2]
            for n in numbers:
                f.write(line + str(n) + '\n')
        if i == 3:
            numbers = indexes[size*2:size*3]
            for n in numbers:
                f.write(line + str(n) + '\n')
        if i == 4:
            numbers = indexes[size*3:size*4]
            for n in numbers:
                f.write(line + str(n) + '\n')
        if i == 5:
            numbers = indexes[size*4:size*5]
            for n in numbers:
                f.write(line + str(n) + '\n')
        if i == 6:
            numbers = indexes[size*5:size*6]
            for n in numbers:
                f.write(line + str(n) + '\n')
        if i == 7:
            numbers = indexes[size*6:size*7]
            for n in numbers:
                f.write(line + str(n) + '\n')   
        if i == 8:
            numbers = indexes[size*7:size*8]
            for n in numbers:
                f.write(line + str(n) + '\n')   
        if i == 9:
            numbers = indexes[size*8:size*9]
            for n in numbers:
                f.write(line + str(n) + '\n')   
        if i == 10:
            numbers = indexes[size*9::]
            for n in numbers:
                f.write(line + str(n) + '\n') 



line = 'porthomclPairsInParalogs.py -t /home/meiker/orthomcl/taxon_list -q /home/meiker/orthomcl/paralogTemp -o /home/meiker/orthomcl/ogenes -p /home/meiker/orthomcl/paralogs -x '

with open(os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl', 'paralogs', today + '_paralogs.sh'), 'w') as f:
    for i in indexes:
                f.write(line + str(i) + '\n')    


'''
8. Run MCL
$cat orthologs/*.tsv >> all.ort.tsv
$mcl all.ort.tsv  --abc -I 1.5 -t 4 -o all.ort.group

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