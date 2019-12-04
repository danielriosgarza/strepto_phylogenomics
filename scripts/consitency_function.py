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
        if score[1] >= 90:
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
    all different writing styles of a string (95 score) and all the same value (first item of list)'''
    all_styles = score_dict(comparison_list)
    synonyms={}
    for k,v in all_styles.items(): #use all scores to make dict that contains all possible writing styles of a species 
        for synonym in v:
            if synonym not in synonyms:
                synonyms[synonym] = k  #gives all same value/writing style
    return synonyms

#def synonym_dict(comparison_score_dict):
#    '''Needs a list with strings that are similar according to their similarity score.
#    Returns a dictionary with all similar strings as keys and all with the same value
#    (first string of the list)'''
#    species={}
#    for k,v in species_all_scores.items(): #use all scores to make dict that contains all possible writing styles of a species 
#        for synonym in v:
#            if synonym not in species:
#                species[synonym] = k  #gives all same value/writing style
#    return species
          
field_l =[]
with open ('/home/meike/strepto_phylogenomics/files/strepto_genomes_quality.tsv') as f:
    for line in f:
        line = line.strip().split('\t')
        if line[0].startswith('genome.genome_id'):
            field_i = line.index('genome.isolation_country')
        else:
            if line[field_i] != '':
                field_l.append(line[field_i])

test = synonym_dict(field_l)
print(test)

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