#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 16:10:12 2020

@author: meike
"""
'''
Making Pieplots of sources associated with certain hosts. 
'''


import matplotlib.pyplot as plt
import os
from pathlib import Path
from collections import Counter 
import pandas as pd
import numpy as np

def make_plot(data, labels, title, savedir):
    '''
    Makes a pie plot with legend at right sight a legend. Colors=Blues
    '''
    #Get color theme for pie
    theme = plt.get_cmap('Spectral')
    colors = [theme(1. * i / len(data)) for i in range(len(data))]
    colors[-1]='#0a5142'
    fig, ax = plt.subplots()
    #autopct='%1.1f%%' gives percentages to pie, pctdistance changes postition of percentages
    ax.pie(data, autopct='%1.1f%%', colors=colors, startangle=0, radius = 1, pctdistance=1.15) 
    
    #loc in combi w/ bbocx --> loc tells matplotlib which part of bounding box should be placed at the arguments of bbox_to_anchor
    ax.legend(labels, loc= "center left", bbox_to_anchor=(1, 0.5))

    ax.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    #pad defines location of title to plot
    plt.title(title, pad=20)
    
    #plt.savefig(os.path.join(savedir), bbox_inches='tight')
    
    
def find_iso_sources_hosts(infile, genus):
    '''
    Looks at hosts and the isolation sources that were used. Saves a pie chart of the 
    top 5 hosts (those with most divers sources) and indicates the rest as "others"
    '''
    with open(lactofile) as f:
        headers = f.readline().strip().split('\t')
        inds = {k: i for i, k in enumerate(headers)}
        hostname_i = inds['host_name']
        isolation_source_i = inds['isolation_source']
        
        #make lists of hosts and the related isolation source
        hosts = []
        source=[]
        for line in f:
            a = line.strip().split('\t')
            hosts.append(a[hostname_i])
            source.append(a[isolation_source_i])
    
    isolations = {k : [] for k in hosts if k != ''}
    
    for i in range(len(hosts)):
        if hosts[i] != '':
            isolations[hosts[i]].append(source[i])
    
    #count the different isolations sources that are associated with certain host
    host_count = {}       
    for k,v in isolations.items():
        host_count[k] = len(set(v))

    top5 = Counter(host_count).most_common(5)
    top5.sort()
     
    
    if len(host_count) > 5: 
        #check the other species and add the isolation sources to a temp list  
        others =[]
        for item in host_count:
            for i in range(len(top5)):
                if item != top5[i][0]:
                    others.append(isolations[item])
        #only select unique isolation sources          
        rest =[]
        for ls in others:
            for it in ls:
                if it not in rest and it != "":
                    rest.append(it)
        top5.append(('Others', len(rest)))
    
    #Defines data and labels    
    numbers = []
    labels =[]
    for item in top5:
        print(item)
        numbers.append(item[1])
        labels.append(item[0])
    
    make_plot(numbers, labels, 'Amount of isolation sources found in different hosts\n('+genus+ ')', os.path.join(p.parents[0], 'figures', '20200109_'+genus+'_isolation_sources_per_host.png'))
    
    
    

path = os.getcwd()
p = Path(path)



lactofile = os.path.join(p.parents[0], 'files', '06012020_lactococcus_database.tsv')
find_iso_sources_hosts(lactofile, "Lactococcus")

florifile = os.path.join(p.parents[0], 'files', '06012020_floricoccus_database.tsv')
find_iso_sources_hosts(florifile, "Floricoccus")

streptofile = os.path.join(p.parents[0], 'files', '06012020_streptococcus_database.tsv')
find_iso_sources_hosts(streptofile, "Streptococcus")