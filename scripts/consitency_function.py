#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 15:02:56 2019

@author: meike
"""
'''Getting out typos etc.'''
'''Makes dicitonary with all possible writing styles as keys asigning to same value to make it consistent.'''
from fuzzywuzzy import process
import dateparser

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

def parse_geographic_location(gl_list):
    '''Parses column geographic location to seperate country from specific location.'''
    new_location_l = []
    for ele in gl_list:
        a = ele.split(':')
        
        if len(a)==2:
            
            while a[1][0]==' ':
                a[1] = a[1][1::]
            new_location_l.append(a[0]+ '('+a[1] + ')')
        else:
            new_location_l.append(a[0]+ '(unkownLocation)')
    return new_location_l


def get_bag_of_terms(field):
    ''' makes sets, seperated by '::', returns list with lists of sets'''    
    name_sets = []
    for line in field:
        if line=='':
            pass
        else:
            name_sets += line.split('::')
    return list(set(name_sets)) #trick to keep only one element of a list (sets are unordered)

def organize_dates(field):
    '''Makes list of strings with different formats of dates and returns a dict of all (complete) dates as strings'''
    dates = {}
    for number in field:
        if any(i.isdigit() for i in number) == True and len(number) > 4:
            if len(number) <= 8 and not'/' in number:
                org_date = str(dateparser.parse(number, date_formats= '%Y/%m'))[:-12]
                dates[number] = org_date
            else:
                org_date = str(dateparser.parse(number, date_formats= '%Y/%m/%d'))[:-9]
                dates[number] = org_date
    return dates


def get_typo_dicts(syn_dict):
    '''Uses synonym dict of all fields to make nested dictionary for spelling check. First key is field to change, second key
    synonyms of that field.'''
    #list of fields that should be looked through for spelling
    fields_to_change = ['genome.additional_metadata', 'genome.biovar',  'genome.collection_date', 
                        'genome.geographic_location', 'genome.habitat', 'genome.host_name', 'genome.isolation_country',
                        'genome.gram_stain', 'genome.isolation_source', 'genome.optimal_temperature', 
                        'genome.publication', 'genome.refseq_accessions']
    
    rv_min = {'-' : ''}
    
    #make dictionary with all the synonym dictionaries inside. Keys = column name
    fc={}
    for name in fields_to_change:
        if name == 'genome.additional_metadata': #if metadata than first get bag of terms
            name_sets = get_bag_of_terms(syn_dict[name])
            
            fc[name] = synonym_dict(name_sets, 95)
            
        elif name =='genome.biovar':
            fc[name] = synonym_dict(syn_dict[name],90)
        
        #make date format consistent with dateparser
        elif name == 'genome.collection_date': 
            fc[name] = organize_dates(syn_dict[name])
        
        #replace positive by '+'
        elif name == 'genome.gram_stain':
            dic = {}
            for item in syn_dict[name]:
                if item != '':
                    dic[item] = '+'
            fc[name] = dic
            
        #all different formats to just numbers
        elif name == 'genome.optimal_temperature':
            temperature = {}
            for item in syn_dict[name]:
                if item != '':
                    if any(i.isdigit() for i in item) == True:
                        tem = ''.join(filter(lambda i: i.isdigit(), item))
                        temperature [item] = str(tem)
                    else:
                        temperature[item] = ''
            fc[name] = temperature
            
        #remove '-' from these fields    
        elif name == 'genome.publication' or 'genome.refseq_accessions':
            fc[name] = rv_min
        
        else:
            fc[name] = synonym_dict(syn_dict[name],95)
    return fc
    
def spelling(inputfile, outputfile, fc):
    '''Replaces typos with consistent spelling. Removes 'genome.' from columns'''
    
    fields_to_change = ['genome.additional_metadata', 'genome.biovar', 'genome.collection_date', 
                        'genome.geographic_location', 'genome.gram_stain','genome.habitat', 'genome.host_name',
                        'genome.isolation_country', 'genome.isolation_source', 'genome.optimal_temperature',
                        'genome.publication', 'genome.refseq_accessions']
    fields_to_change.sort()
    
    with open (inputfile) as f:
        with open(outputfile, 'w') as f2:
            headers = f.readline().strip().split('\t')
            
            #look for indexes of the fields that need speeling check and add them to a list
            indexes = [i for i,name in enumerate(headers) if name in fields_to_change]
            print(indexes)
            #remove the 'genome.' from the headers
            for field in headers:
                if 'genome.' in field:
                    f2.write(field.replace('genome.',''))
                else:
                    f2.write(field)
                f2.write('\t')
            f2.write('\n')
            
            #go through lines and make consistent spelling in desired columns
            for line in f:
                a = line.strip().split('\t')
                for i,name in enumerate(a):
                    #look if the current line index is in the list of indexes of interest
                    if i in indexes:
                        #look for the name of the field
                        field = fields_to_change[indexes.index(i)] 
                    
                        
                        if "additional_metadata" in field:
                            
                            if a[i] == '':
                                pass
                            else:
                                items = []
                                items += a[i].split('::')
                                terms=[]
                                for item in items:
                                    if item in fc[field]:
                                        c = fc[field][item]
                                        terms.append(c)
                                    else:
                                        terms.append(item)
                                terms.sort()
                                f2.write("::".join(terms))
                                    
                                    
                        #if the name of the field is in the dictionary of synonym dictionaries
                        elif name in fc[field]:
                            f2.write(fc[field][name]) #write what in the dicts of fields in the synonym dict is
                        else:
                            f2.write(name) #if not just write whats there
                    else:
                        f2.write(name)
                    f2.write('\t') #seperate each column with tab
                f2.write('\n') #after each line
                    

#get typos out of the fields   
lacto_synonyms = get_synonyms('/home/meike/strepto_phylogenomics/files/lactococcus_genomes_quality.tsv')
fields2change = get_typo_dicts(lacto_synonyms)
lacto_input = '/home/meike/strepto_phylogenomics/files/lactococcus_genomes_quality.tsv'
lacto_output = '/home/meike/strepto_phylogenomics/files/lactococcus_genome_database.tsv'


spelling(lacto_input,lacto_output, fields2change)

meta_sets = get_bag_of_terms(lacto_synonyms['genome.additional_metadata'])
meta_dict = synonym_dict(meta_sets, 95)

