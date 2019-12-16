#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 13:24:53 2019

@author: meike
"""
import pandas as pd
from fuzzywuzzy import process

def compare(str2Match, strOptions, score_t):
    '''Give query string that is compared to other strings in a list. Returns list with all strings that had
    a score of at least the chosen score (Levenshtein Distance). No identical hits are in the list included.'''
    score_l = []
    names=[]
    scores = process.extract(str2Match, strOptions, limit = len(strOptions))
    for score in scores:
        if score[1] >= score_t:
            score_l.append(score[0])
    for item in score_l:
        if item not in names:
            names.append(item)
    return names

def score_dict(comparison_list, score_t):
   '''Needs a list of strings that should be compared and uses compare function. Returns a dict with
   all strings as keys and all different strings with a high score as values'''
   all_scores ={} 
   for name in comparison_list: #make dict with field name as key and all scores higher than 95 as values
            if name not in all_scores:
                all_scores[name] = compare(name, comparison_list, score_t)
   return all_scores

def synonym_dict(comparison_list, score_t):
    '''Needs a list of strings that should be compared and uses the function score_dict. Returns a dict with
    all different writing styles of a string  and all the same value (first item of list)'''
    all_styles = score_dict(comparison_list, score_t)
    synonyms={}
    for k,v in all_styles.items(): #use all scores to make dict that contains all possible writing styles of a species 
        for synonym in v:
            if synonym not in synonyms:
                synonyms[synonym] = k  #gives all same value/writing style
    return synonyms

def get_synonyms(file):
    '''Give file with different columns and scores that you wish for each column. Goes through each 
    column and stores all values in a dict with as key the column name (empty lines not included)'''
    table= pd.read_csv(file, sep = '\t')
    table = table.astype(str)
    headers =list(table.columns)
    synonyms = {i:[] for i in headers}
    for k in synonyms:
        values = table[k].tolist()
        temp =[]
        for item in values:
            if item != 'nan':
                temp.append(item)
        synonyms[k] = temp             
    return synonyms

def parse_sets(set_l):
    ''' makes sets, seperated by '::', returns list with lists of sets'''    
    name_sets = []
    for name in set_l:
        if len(name) >= 2:
            word_l = name.split('::')
            name_sets.append(word_l)
        else:
            name_sets.append(name)

    return name_sets

#in list per item: string mit jedem string in eigener liste vergleichen, dann mit jedem str von der nachsten etc.
def compare_sets(sets_l, score_t):
    '''Loops through list of sets and looks at each set and compares it to the other sets. Needs a list with in lists 
    all sets'''
    all_scores ={}
    sets =[]
    for name_l in sets_l: 
        for name in name_l:
            if name not in sets and name != '':
                sets.append(name)
    for name in sets: #make dict with set name as key and all scores higher than 95 as values
        if name not in all_scores:
            all_scores[name] = compare(name, sets, score_t)
    return all_scores
    

lacto_synonyms = get_synonyms('/home/meike/strepto_phylogenomics/files/lactococcus_genomes_quality.tsv')
lacto_input = '/home/meike/strepto_phylogenomics/files/lactococcus_genomes_quality.tsv'
lacto_output = '/home/meike/strepto_phylogenomics/files/lactococcus_genome_database.tsv'

names_sets = parse_sets(lacto_synonyms["genome.additional_metadata"])


      
                            
                            
def spelling(inputfile, outputfile, syn_dict):
    '''Replaces typos with consistent spelling. Changes columns ['genome.biovar',  'genome.geographic_location', 'genome.habitat', 'genome.host_name',
    'genome.isolation_country']. Removes 'genome.' from columns'''
    
    #list of fields that should be looked through for spelling
    fields_to_change = ['genome.additional_metadata', 'genome.biovar',  'genome.geographic_location', 'genome.habitat', 'genome.host_name',
    'genome.isolation_country']
    
    #make dictionary with all the synonym dictionaries inside. Keys = column name
    fc={}
    for name in fields_to_change:
        if name == 'genome.additional_metadata':
            name_sets = parse_sets(syn_dict[name])
            fc[name] = compare_sets(name_sets, 95)
        if name =='genome.biovar':
            fc[name] = synonym_dict(syn_dict[name],90)
        else:
            fc[name] = synonym_dict(syn_dict[name],95)
    
    
    with open (inputfile) as f:
        with open(outputfile, 'w') as f2:
            headers = f.readline().strip().split('\t')
            
            #look for indexes of the fields that need speeling check and add them to a list
            indexes = [i for i,name in enumerate(headers) if name in fields_to_change]
            
            #remove the 'genome.' from the headers
            for field in headers:
                if 'genome.' in field:
                    f2.write(field.replace('genome.',''))
                else:
                    f2.write(field)
                f2.write('\t')
            f2.write('\n')
            
            #go through lines and make consistent spelling in desireed columns
            for line in f:
                a = line.strip().split('\t')
                for i,name in enumerate(a):
                    #look if the current line index is in the list of indexes of interest
                    if i in indexes:
                        #look for the name of the field
                        field = fields_to_change[indexes.index(i)] 
                        
                        #if the name of the field is in the dictionary of synonym dictionaries
                        if name in fields_to_change:
                            #if column contains sets, change sets
                            if field == 'genome.additional_metadata':
                                b = a[i].split('::')
                                print(b)
                                for word_combi in b:
                                    f2.write(fc[field][word_combi] + '::')
                                
                            #if not sets: write what in synonym_dict is
                            else:
                                f2.write(fc[field][name]) #write what in the dicts of fields in the synonym dict is
                        #if name not in dict, write whats there
                        else:
                            f2.write(name) #if not just write whats there
                    #if not a field to change, just write the line        
                    else:
                        f2.write(name)
                    f2.write('\t') #seperate each column with tab
                f2.write('\n') #after each line
                    
lacto_synonyms = get_synonyms('/home/meike/strepto_phylogenomics/files/lactococcus_genomes_quality.tsv')
lacto_input = '/home/meike/strepto_phylogenomics/files/lactococcus_genomes_quality.tsv'
lacto_output = '/home/meike/strepto_phylogenomics/files/lactococcus_genome_database.tsv'
spelling(lacto_input,lacto_output, lacto_synonyms)
