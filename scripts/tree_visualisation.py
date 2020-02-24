#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 14:21:00 2020

@author: meike
"""

'''
Visualisation and analysis of the tree.
'''

import ete3
import os
from pathlib import Path
from datetime import date

#get the date to keep track of the scripts (added to scriptname)
today = date.today().strftime("%d/%m/%Y")
today = today.split('/')
today = ''.join(today)

    
path = os.getcwd()
p = Path(path)


files_dir = os.path.join(p.parents[0], 'files', 'phylogenetic_tree/') 

hmm_contree = ete3.Tree(files_dir + '/concat_alignments.contree')

ts = ete3.TreeStyle()
ts.mode = 'c'
#ts.show_branch_support = True

hmm_contree.render(files_dir + today + '_firstTree.pdf', tree_style=ts)