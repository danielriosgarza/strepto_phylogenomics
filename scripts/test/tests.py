#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 09:16:34 2019

@author: meike
"""

#--outdir /home/meiker/git/prokka_annotation/streptococcus_02451
# 'streptococcus_11895', 'streptococcus_11939']
# 'lactococcus_00188'
#wget ftp://ftp.patricbrc.org/genomes/1739284.3/1739284.3.fna -O /home/meiker/git/genomes/streptocuccus_11846.fna


from fuzzywuzzy import process
import dateparser
import re


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
                synonyms[synonym] = k.capitalize()  #gives all same value/writing style
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

def get_bag_of_terms(field, splitparameter):
    ''' makes sets, seperated by splitparamter, returns list with lists of sets'''    
    name_sets = []
    for line in field:
        if line=='':
            pass
        else:
            name_sets += line.split(splitparameter)
    return list(set(name_sets))

def organize_dates(field):
    '''Makes list of strings with different formats of dates and returns a dict of all (complete) dates as strings'''
    dates = {}
    for number in field:
        if any(i.isdigit() for i in number) == True and len(number) > 4:
            if len(number) <= 8 and not'/' in number:
                org_date = str(dateparser.parse(number, date_formats= '%Y/%m'))[:-12]
                dates[number] = org_date
            elif len(number) == 9 and '/' in number:
                pass
            else:
                org_date = str(dateparser.parse(number, date_formats= '%Y/%m/%d'))[:-9]
                dates[number] = org_date
    return dates

def organize_ages(field):
    
    ages = {}
    for item in field:
        item = item.lower()
        if any(i.isdigit() for i in item) == True:
            number = ''.join(re.findall(r'\d+', item))
            unit = ''.join(re.findall(r'[a-z]', item))
            number_count = sum(c.isdigit() for c in item)
            if 'm' in unit:
                time = 'M'
            else:
                time = 'Y'
            
            if item not in ages:
                if number_count >= 4:
                    ages[item] = item+time
                elif item.startswith('>') or item.startswith('<'):
                    ages[item] = item+time
                else:
                    ages[item] = number+time
    return ages
                
    
    
strepto_syn = get_synonyms('/home/meike/strepto_phylogenomics/files/streptococcus_genomes_quality.tsv')

strepto_input = '/home/meike/strepto_phylogenomics/files/streptococcus_genomes_quality.tsv.tsv'
strepto_output = '/home/meike/strepto_phylogenomics/files/streptococcus_genome_database.tsv'

hh_sets = get_bag_of_terms(strepto_syn['genome.host_health'], '; ')
hh = synonym_dict(hh_sets, 98)
