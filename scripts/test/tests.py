#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 09:16:34 2019

@author: meike
"""
import random

def split_files(db_ids):
    '''
    Splits db_list depending on length and makes lists with savdir that can be used for bash file generation.
    Tuples w/ (db_id, original_index)
    '''    
    groupsize = int(len(db_ids)/6)
    dbwithi = [(id_, i) for i, id_ in enumerate(db_ids)]
    
    id_i = [dbwithi[i :i +groupsize] for i in range(0, len(dbwithi), groupsize)]
           
    for i in range(0, len(db_ids), groupsize):
        id_i.append([db_ids(db_ids[i],db_ids.index(db_ids[i])))   
    
    [lst[i:i + n] for i in range(0, len(lst), n)]
    
    return id_i


def histogram(item_l, title=None, xlabel=None, label=None, bins=10, color='#2166ac', rotation=0):
    '''
    Makes density histogram from given list.
    '''
    plt.hist(item_l, bins, density=1, alpha=0.75, color = color, label =label)
    plt.title(title, fontsize= 14)
    plt.xlabel(xlabel, fontsize= 12)
    plt.xticks(rotation=rotation)
    plt.ylabel("Density", fontsize= 12)
    plt.tight_layout() 
    
    #plt.savefig(os.path.join(p.parents[0], 'figures', )

test = 