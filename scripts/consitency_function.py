#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 15:02:56 2019

@author: meike
"""
'''Getting out typos etc.'''
'''Makes dicitonary with all possible writing styles as keys asigning to same value to make it consistent.'''
from fuzzywuzzy import process

def compare(str2Match, strOptions):
    '''Give query string that is compared to other strings in a list. Returns list with all strings that had
    a score of at least 90% (Levenshtein Distance). No identical hits are in the list included.'''
    score_l = []
    names=[]
    scores = process.extract(str2Match, strOptions, limit = len(strOptions))
    for score in scores:
        if score[1] >= 95:
            score_l.append(score[0])
    for item in score_l:
        if item not in names:
            names.append(item)
    return names

def score_dict(comparison_list):
   '''Needs a list of strings that should be compared and uses compare function. Returns a dict with
   all strings as keys and all different strings with a high score as values'''
   all_scores ={} 
   for name in comparison_list: #make dict with field name as key and all scores higher than 95 as values
            if name not in all_scores:
                all_scores[name] = compare(name, comparison_list)
   return all_scores

def synonym_dict(comparison_list):
    '''Needs a list of strings that should be compared and uses the function score_dict. Returns a dict with
    all different writing styles of a string  and all the same value (first item of list)'''
    all_styles = score_dict(comparison_list)
    synonyms={}
    for k,v in all_styles.items(): #use all scores to make dict that contains all possible writing styles of a species 
        for synonym in v:
            if synonym not in synonyms:
                synonyms[synonym] = k  #gives all same value/writing style
    return synonyms

def get_synonyms(file):
    '''Give file with different columns. Goes through each coloum and stores all values in a dict with as key the
    column name (empty lines not included)'''
    with open (file) as f:
        headers =  f.readline().strip().split('\t')
        synomyms = {i:[] for i in headers} #make dict to store all lists with synonyms
        for line in f:
            line = line.strip().split('\t')
            if len(line) != len(headers):
                difference = len(headers)-len(line)
                for i in range(difference):
                    line.append('')
            for i, name in enumerate(headers): #go through headers and store lines in different lists to check for typos
                if line[i].isspace() == False and line[i] != '':
                        synomyms[name].append(line[i])
    return synomyms
#get typos out of the fields
    
lacto_synonyms = get_synonyms('/home/meike/strepto_phylogenomics/files/lactococcus_genomes_quality.tsv')
columns = ['genome.biovar', 'genome.cell_shape', 'genome.geographic_location', 'genome.habitat', 'genome.host_name',
           'genome.isolation_country']

for col in columns:
    col = synonym_dict(lacto_synonyms[col])
print(col)
lines =[]
with open ('/home/meike/strepto_phylogenomics/files/lactococcus_genomes_quality.tsv') as f:
    for line in f:
        line = line.strip().split('\t')
        lines.append(line)

with open('/home/meike/strepto_phylogenomics/files/lactococcus_genome_database.tsv', 'w') as f:
    for line in f:
        for i, column in enumerate(lines):
            if 'biovar' in column:
                if line[i] in lacto_synonyms:
                    line[i] = lacto_synonyms[column]
        f.write('\t'.join(line) + '\n')
            
#def compare_sets(sets_list):
#    #fields that can be merged --> token_sort_ratio?
#    for name in sets_list:
    
# set seperated by : and from other sets by ::
#
#with open ('/home/meike/strepto_phylogenomics/files/lactococcus_genomes_quality.tsv') as f:
#    headers = f.readline().strip().split('\t')
#    additional_metadata = []
#    for line in f:
#        line = line.strip().split('\t')
#        if len(line) != len(headers):
#            difference = len(headers)-len(line)
#            for i in range(difference):
#                line.append('')
#        for i, name in enumerate(headers):
#            if 'additional_metadata' in name:
#                additional_metadata.append(line[i])
#
#for i, item in enumerate(additional_metadata):
#    items = item.split('::')
#    additional_metadata[i] = items
#    #list with containg list with sets; compare sets of first list to rest of sets in other lists
#for item in additional_metadata:
#    for i in item:
#        
     
with open ('/home/meike/strepto_phylogenomics/files/lactococcus_genomes_quality.tsv') as f:
    headers = f.readline().strip().split('\t')
    data_dict = {i:[] for i in headers}
    #make comparison lists to remove typos
    additional_meta =[]
#    anti_resistance = []
#    assembly_method = []
#    body_sample_site = []
#    collection_date = []
#    collection_year
    for line in f:
        info = line.strip().split('\t')
        if len(info) != len(headers):
            difference = len(headers)-len(info)
            for i in range(difference):
                info.append('')
        for inde, name in enumerate(headers):
            data_dict[name] = info[inde]
            
    
#fields_ids = {}          
#field_l =[]
#temp =[]
#with open ('/home/meike/strepto_phylogenomics/files/lactococcus_genomes_quality.tsv') as inputfile:
#    with open ('/home/meike/strepto_phylogenomics/files/lactococcus_database.tsv', 'w') as outfile:
#        for line in inputfile:
#            line = line.strip().split('\t')
#            if line[0].startswith('genome.genome_id'):
#                for i, header in enumerate(line):
#                    fields_ids[header] = i
#                    if "antimicrobial" in fields_ids:
#                        temp.append(line[i])
#            outfile.writelines('\t'.join(line) + '\n')
#            
#
#test = synonym_dict(field_l)
#print(test)

#fields that can be merged --> token_sort_ratio?
#Str1 = "united states v. nixon"
#Str2 = "Nixon v. United States"
#Ratio = fuzz.ratio(Str1.lower(),Str2.lower())
#Partial_Ratio = fuzz.partial_ratio(Str1.lower(),Str2.lower())
#Token_Sort_Ratio = fuzz.token_sort_ratio(Str1,Str2)
#print(Ratio)
#print(Partial_Ratio)
#print(Token_Sort_Ratio)
#OUT:
#59
#74
#100