#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 14:21:00 2020

@author: meike
"""

'''
Visualisation and analysis of the tree.
Tests with ete3.
'''

from ete3 import Tree, NodeStyle, TreeStyle
import ete3
import os
from pathlib import Path
from datetime import date
import matplotlib
import matplotlib.pyplot as plt
import itertools


def get_colors(data_list, color='Reds'):
    '''
    Needs a dat_list to determine how many different colours are needed. Cmap colour schemes can be provided, default are Reds.
    '''
    colors = []
    cmap = plt.get_cmap(color, len(data_list))
    
    for i in range(cmap.N):
        rgb = cmap(i)[:3] # will return rgba, we take only first 3 so we get rgb
        colors.append(matplotlib.colors.rgb2hex(rgb))
    return colors

def get_ids2species(speciesD):
    '''
    Get a dict with all ids as keys and the corresponding species as values.
    '''
    id_2_species = {}
    
    for genus in speciesD:
        for spec, ids in speciesD[genus].items():
            for _id in ids:
                id_2_species[_id] = spec
    return id_2_species
    
def species2ids(inputfile):
    '''
    Get all ids that belong to a certain species in a dict. Needs file containing species and id (first column of file).
    '''
    all_species = {}
    with open (inputfile) as f:
        headers = f.readline().strip().split('\t')
        species_ind = headers.index('species')
        for line in f:
            a = line.strip().split('\t')
            species = a[species_ind]
            if species not in all_species:
                all_species[species] = [a[0]] #assign id as value, type list
            else:
                all_species[species].append(a[0]) #if another id is the same species add the id
    return all_species


path = os.getcwd()
p = Path(path)
files_dir = os.path.join(p.parents[0], 'files')

#make dict with all ids belonging to a ceratin species sorted per genus
#Structure: genus --> species --> ids
species = {}
species['streptococcus'] = species2ids(files_dir+'/03032020_streptococcus_database_final.tsv')
species['lactococcus'] = species2ids(files_dir+'/06012020_lactococcus_database.tsv')
species['floricoccus'] = species2ids(files_dir+'/06012020_floricoccus_database.tsv')


#make dict to assign colors for leaves
leaf_colours = {}

streptos = list(species['streptococcus'].keys())
lactos = list(species['lactococcus'].keys())
floris =list(species['floricoccus'].keys())

strep_color = get_colors(streptos)
lacto_color = get_colors(lactos, color = 'Blues')
flori_color = get_colors(floris, color = 'Greens')

for i, spec in enumerate(streptos):
    for _id in species['streptococcus'][spec]:
        leaf_colours[_id] = strep_color[i]   

for i, spec in enumerate(lactos):
    for _id in species['lactococcus'][spec]:
        leaf_colours[_id] = lacto_color[i]  
     
for i, spec in enumerate(floris):
    for _id in species['floricoccus'][spec]:
        leaf_colours[_id] = flori_color[i]  


#Set style to tree
#t = Tree(files_dir + '/phylogenetic_tree/concat_alignments.contree')
t = Tree('/home/meike/tests/Files/09032020_reduced_concat_alignments.fa.contree')        
        
for node in t.traverse():
    node.img_style["size"] = 0 #removes dots at nodes  
    if node.is_leaf():
        color = leaf_colours.get(node.name, None)
        if color:
            node.img_style['bgcolor'] = color
            
ts = TreeStyle()
ts.mode = 'c'
#ts.show_leaf_name = False

t.show(tree_style=ts)

ids2species = get_ids2species(species)

#%%

#t = Tree(files_dir + '/phylogenetic_tree/concat_alignments.contree') 

ids2species = get_ids2species(species)

keepNodes = []
for node in t.traverse():
    children_nodes = node.children
    child2remove = []
    temp = []
    for child in children_nodes:
        #check if children are leafnodes and add the according species in a temp list
        if any(i.isdigit() for i in child.name):
                temp.append(ids2species[child.name])
                child2remove.append(child)
        #also add the node that is not a leaf
        else:
            temp.append(child.name)
    #check if the children of a given node are all leafnodes and belong to the same species
    if len(set(temp)) == 1 and any(c.isalpha() for c in temp[0]):
        keepNodes.append(child2remove[0].name)
        # for kid in child2remove:
        #     node.remove_child(kid)
        # node.name = child2remove[0].name
    elif child.is_leaf() and any(c.isalpha() for c in child.name):
            keepNodes.append(child.name)
        
for node in t.traverse():
    node.img_style["size"] = 0 #removes dots at nodes  
    if node.is_leaf():
        color = leaf_colours.get(node.name, None)
        node.name = ids2species[node.name]
        if color:
            node.img_style['bgcolor'] = color
            
t.prune(keepNodes)            
ts = TreeStyle()
ts.mode = 'c'            
t.show(tree_style = ts)

#t.prune([list]) list with all leafs that should be kept
            

#%%
#Set outgroups to root the tree

#t = Tree(files_dir + '/phylogenetic_tree/concat_alignments.contree')       

# #determine the two lactococcus nodes with greatest distance to search for the root node
lactoleaves = []
florileaves = []
for leaf in t.iter_leaves():
    if 'lacto' in leaf.name:
        lactoleaves.append(leaf)
    if 'flori' in leaf.name:
        florileaves.append(leaf)

# #look for all distances (takes some time (total 17578 combinations))
distances = []
order = []
for a, b in itertools.combinations(lactoleaves, 2):
    distances.append(a.get_distance(b))
    order.append((a,b))

# #get the pair (two nodes) that have greatest distance
greatest_dist = max(distances)
index = distances.index(greatest_dist)
pair_for_outgroup = order[index] 
print(pair_for_outgroup) # 'lactococcus_00009', 'lactococcus_00081'

# #search for common ancestor of lactococcus and floricoccus
# florinodes = [florileaves[0].name, florileaves[1].name]
# lactonodes = [pair_for_outgroup[0].name, pair_for_outgroup[1].name]

#root_node = t.get_common_ancestor('floricoccus_00001', 'floricoccus_00002', 'lactococcus_00009', 'lactococcus_00081')

root_node = t.get_common_ancestor('floricoccus_00001', 'lactococcus_00009', 'lactococcus_00173')

t.set_outgroup(root_node)

for node in t.traverse():
    node.img_style["size"] = 0 #removes dots at nodes  
    if node.is_leaf():
        #node.img_style['draw_descendants'] = False
        color = leaf_colours.get(node.name, None)
        node.name = ids2species[node.name]
        if color:
            node.img_style['bgcolor'] = color
        # if 'flori' in node.name:
        #     genus = ete3.RectFace(width = 50, height = 50, bgcolor = '#d1e5f0', fgcolor = '#d1e5f0')
        #     node.add_face(genus, column = 0, position="aligned")
        #     node.img_style['size'] = 50
        #     node.img_style['hz_line_width'] = 30
        #     node.img_style['bgcolor'] = 'Red'
        # if 'lacto' in node.name:
        #     node.img_style['size'] = 50
        #     node.img_style['hz_line_width'] = 30
        #     node.img_style['bgcolor'] = 'MediumBlue'

species_counter = {}

for i in t.iter_leaf_names():
    if i not in species_counter:
        species_counter[i] = 1
    else:
        species_counter[i] += 1

t.write(format = 1, outfile = '/home/meike/tests/Files/rooted_tree.nwk')
            
ts = TreeStyle()
ts.mode = 'c'
#ts.show_leaf_name = False

t.show(tree_style=ts)
        

#%%


t = Tree('/home/meike/tests/Files/09032020_reduced_concat_alignments.fa.contree') 

nodeId = 0

for n in t.traverse('levelorder'):
    n.add_features(ND = nodeId)
    nodeId += 1
    print(n.name, )

        
    

for node in t.traverse():
    node.img_style["size"] = 0 #removes dots at nodes  
    if node.is_leaf():
        color = leaf_colours.get(node.name, None)
        if color:
            node.img_style['bgcolor'] = color
            
         
ts = TreeStyle()
ts.mode = 'c'            
t.show(tree_style = ts)

#%%
from operator import itemgetter



t = Tree(files_dir + '/phylogenetic_tree/concat_alignments.contree') 

#counter to see how many branches the tree has before removing leaves
before = 0
for node in t.iter_leaves():
    before += 1
    
#dict that maps ids to its species
ids2species = get_ids2species(species) 

#counter to keep track how many leaves will be removed
check = 0


for node in t.traverse():
    t2 = Tree(t.write())
    
    children_nodes = node.children
    species_check = []
    
    #iterate over nodes and check the leaves is the same species is present
    for c in children_nodes:
        if c.is_leaf():
            species_check.append((c, ids2species[c.name]))
    #sorted(species_check, key = itemgetter(1)) #sort according to species, so all ids/nodes from one species are in a row
    
    #go through list and store all nodes/leaves where species is already seen in a list
    nodes2detach = []
    seen = set()
    for ele in species_check:
        if ele[1] not in seen:
            seen.add(ele[1])
        else:
            nodes2detach.append(ele[0])
    
    #if list of nodes to detach is not empty go through list and remove the nodes (duplicates on same branch)
    if len(nodes2detach) > 0:
        for item in nodes2detach:
            t2.item.detach()
            check += 1
                
#counter to keep track how many leaves the tree still has
after = 0
for node in t.iter_leaves():
    after += 1    
print(before, after, check)    

#set style for the tree
for node in t.traverse():
    node.img_style["size"] = 0 #removes dots at nodes  
    if node.is_leaf():
        color = leaf_colours.get(node.name, None)
        if color:
            node.img_style['bgcolor'] = color
            
         
ts = TreeStyle()
ts.mode = 'c'            
t.show(tree_style = ts)
        