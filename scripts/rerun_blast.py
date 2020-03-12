#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 12:00:27 2020

@author: meike
"""

'''
Rerun blast on cluster.
'''
import os
from pathlib import Path
from datetime import date

def check_last_blasted_id(file, List_already_done):
    '''
    Checks which id was the last that was blasted (created) and the next that is missing in the file.
    '''    
    not_done = ''
    with open(file) as f:
        for line in f:
            a = line.strip().split('/')
            if a[-1] not in already_done:
                return not_done
                break
            not_done = a[-1]
               

#get the date to keep track of the scripts (added to scriptname)
today = date.today().strftime("%d/%m/%Y")
today = today.split('/')
today = ''.join(today)       
        
already_done = []
with open('/home/meike/tests/Files/already_blasted.txt') as f:
    for line in f:
        line = line.strip().split()
        already_done += line

    
path = os.getcwd()
p = Path(path)

plinc_folder = os.path.join(p, 'bash_scripts', 'porthomcl', 'plinc')
scripts = os.listdir(plinc_folder)  

plinc_paths = []
for scriptname in scripts:
    plinc_paths.append(plinc_folder + '/' + scriptname)

i = 0    
for script in plinc_paths:
    last_id = check_last_blasted_id(script, already_done)
    i += 1
    with open(script) as f:
        with open(plinc_folder + '/' + today + '_blastrun_' + str(i) + '.sh', 'w') as f2:
            passed_last_id = False
            for line in f:
                a = line.strip().split('/')
                if a[-1] == last_id:
                    passed_last_id = True
                    
                if passed_last_id:
                   f2.write(line)


script_folder = os.path.join(p, 'bash_scripts', 'porthomcl', 'narrativum2') 
s_scripts = os.listdir(script_folder) 

surprise_paths = []
for scriptname in s_scripts:
    surprise_paths.append(script_folder + '/' + scriptname)

i = 0    
for script in surprise_paths:
    last_id = check_last_blasted_id(script, already_done)
    i += 1
    with open(script) as f:
        with open(script_folder + '/' + today + '_blastrun_' + str(i) + '.sh', 'w') as f2:
            passed_last_id = False
            for line in f:
                a = line.strip().split('/')
                if a[-1] == last_id:
                    passed_last_id = True
                    
                if passed_last_id:
                   f2.write(line)          

