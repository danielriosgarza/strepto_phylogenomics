#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  2 13:05:10 2020

@author: meike
"""

'''
Getting information over gene lists using Enrichr. 
(used library = BioPlanet2019)

Runcell 1: Get information over associated paths to different Pan-genome parts

Runcell 2: get info over associated pathways to potential interesting genes per group (defined based on phylogenetic tree)
'''

import os
from pathlib import Path
import json
import requests
from datetime import date

def analyze_List(list_genes, description_genelist):
    '''
    Needs list of genes and description for Enrichr.
    Returns Dict with shortId and userListId (needed to download information)

    '''
    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/addList'
    genes_str = '\n'.join(list_genes)
    description = description_genelist
    payload = {
        'list': (None, genes_str),
        'description': (None, description)
    }
    
    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')
    data = json.loads(response.text)
    return data

def get_enrichment_results(ListID, Library):
    '''
    Give list id and the desired library. Returns enrichment results as Dict.
    '''
    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
    query_string = '?userListId=%s&backgroundType=%s'
    user_list_id = ListID
    gene_set_library = Library
    response = requests.get(
        ENRICHR_URL + query_string % (user_list_id, gene_set_library)
     )
    if not response.ok:
        raise Exception('Error fetching enrichment results')
    
    data = json.loads(response.text)
    return data

def download_enrichment_results(ListID, Library, file):
    '''
    Saves the enrichment results into given dir. Needs list id and chosen library that should be downloadded.
    '''
    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/export'
    query_string = '?userListId=%s&filename=%s&backgroundType=%s'
    user_list_id = ListID
    filename = file
    gene_set_library = Library
    
    url = ENRICHR_URL + query_string % (user_list_id, filename, gene_set_library)
    response = requests.get(url, stream=True)
    
    with open(filename + '_' + Library + '.tsv', 'wb') as f:
        for chunk in response.iter_content(chunk_size=1024): 
            if chunk:
                f.write(chunk)
            
                    
lib = 'BioPlanet_2019'

path = os.getcwd()
p = Path(path)

#%% runcell 1
#information over associated pathways in different pan-genome parts

bin_file = os.path.join(p.parents[0], 'files', 'binary_table', '14042020_binary_table_sorted.tsv')

saveDir = os.path.join(p.parents[0], 'files', 'pan-genome') + '/'

#Counts of genes belongin to different parts of the pan-genome
pan_count = {}
pan_genes = {}
with open(bin_file) as f:
    f.readline()
    for line in f:
        a = line.split('\t')
        genes = a[2].split(',')
        if 'N/A' in genes:
            new_genes = []
            for i in genes:
                if i != 'N/A':
                    new_genes.append(i)
            genes = new_genes
            
        if a[0] == 'Extended Core':
            a[0] = 'Extended_Core'
        
        if a[0] not in pan_count:
            pan_count[a[0]] = 1
        else:
            pan_count[a[0]] += 1

        if a[0] not in pan_genes:
            pan_genes[a[0]] = genes
        else:
            pan_genes[a[0]] = pan_genes[a[0]] + genes
print(pan_count)

for g in pan_genes:
    data = analyze_List(pan_genes[g], g)
    l_id = data['userListId']
    download_enrichment_results(l_id, lib, saveDir + g)
    

    
#%% runcell 2

#Look for associated pathways of potential interesting genes


#get the date to keep track of the scripts (added to scriptname)
today = date.today().strftime("%d/%m/%Y")
today = today.split('/')
today = ''.join(today)

bin_greater = os.path.join(p.parents[0], 'files', 'binary_table', '24042020_pvalues_genes_greater.tsv')

bin_less = os.path.join(p.parents[0], 'files', 'binary_table', '24042020_pvalues_genes_less.tsv')

group2grless2genes = {i : { 'greater': [], 'less' : []} for i in range(8)}

#determine number of cols and save the lines
with open(bin_greater) as f:
    f.readline()
    for line in f:
        a = line.strip().split('\t')
        gene_num = a[0].split('_')[1]
        for ind, value in enumerate(a[1::]):
            if float(value) <= 0.05:
                group2grless2genes[ind]['greater'].append(gene_num)
                
with open(bin_less) as f:
    f.readline()
    for line in f:
        a = line.strip().split('\t')
        gene_num = a[0].split('_')[1]
        for ind, value in enumerate(a[1::]):
            if float(value) <= 0.05:
                group2grless2genes[ind]['less'].append(gene_num)
                
#Dict: group --> less/greater --> gene names
g_names_needed = []
for gr in group2grless2genes:
    for di in group2grless2genes[gr]:
        for num in group2grless2genes[gr][di]:
            if num not in g_names_needed:
                g_names_needed.append(num)
g_names_needed = sorted(g_names_needed)


#get gene name (in binary table only number)
genes = {}
with open(os.path.join(p.parents[0], 'files', 'binary_table', '14042020_binary_table_sorted.tsv')) as f:
    f.readline()
    lines = [l for l in f]

for ind in g_names_needed:
    ge = lines[int(ind)].strip().split('\t')[2]
    ge = ge.split(',')
    genes[ind] = ge
        
#get gene lists per group (greater and less presence apart)
groups_genes = {}   

for gr in group2grless2genes:
    groups_genes[gr] = {}
    for Dir in group2grless2genes[gr]:
        groups_genes[gr][Dir] = []
        for num in group2grless2genes[gr][Dir]:
            name = genes[num]
            if 'N/A' not in name:
                groups_genes[gr][Dir] += name
 
#make folders to save the enrichr results in               
for i_gr, gr in enumerate(groups_genes):
    folder = os.path.join(p.parents[0], 'files', 'pathway_search','group' + str(i_gr + 1))
    if not os.path.isdir(folder):
        os.mkdir(folder)
    for grless in groups_genes[gr]:
        
        #where enrichr results should be saved to
        filename = folder + '/' + grless 
        #save genes in a list as input for enrichr
        gene_l = []
        for ge_name in groups_genes[gr][grless]:
            gene_l.append(ge_name)
            
        #search and save results of enrichr
        data = analyze_List(gene_l, lib)
        l_id = data['userListId']
        download_enrichment_results(l_id, lib, filename)
    