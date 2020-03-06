#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 11:28:04 2020

@author: meike
"""
'''
Rename node/adding IDs to nodes
'''


import ete3
from ete3 import Tree, TreeStyle, NodeStyle
import matplotlib.pyplot as plt
import matplotlib
import os
from pathlib import Path
from skbio import TreeNode


path = os.getcwd()
p = Path(path)


def collapsed_leaf(node):
    if len(node2labels[node]) == 1:
       return True
    else:
       return False

t = Tree("((((a,a,a)a,a)aa, (b,b)b)ab, (c, (d,d)d)cd);", format=1)


i = 1
for node in t.traverse('levelorder'):
    if node.name is not None and node.name != '':
        node.name = '%s:N%d' % (node.name, i)
    else:
        node.name = 'N%d' % i
    i += 1

t.write(outfile = '/home/meike/tests/Files/tree/test.nwk')

tree = TreeNode.read('/home/meike/tests/Files/tree/test.nwk')

i = 1
for node in tree.levelorder():
    if not node.is_tip():
        if node.name is not None and node.name != '':
            node.name = '%s:N%d' % (node.name, i)
        else:
            node.name = 'N%d' % i
        i += 1

t2 = tree.write('/home/meike/tests/Files/tree/test_with_ids.nwk')


# We create a cache with every node content
node2labels = t.get_cached_content(store_attr="name")

collapsed_nodes = []
for node in t.traverse():
    if collapsed_leaf(node) == 1:
        collapsed_nodes.append(node.name)

with open ('/home/meike/tests/Files/tree/collapsing.txt', 'w') as f:
    f.write('COLLAPSE\nDATA\n')
    for name in collapsed_nodes:
        f.write(name + '\n')

#t.show()
#t.write(is_leaf_fn=collapsed_leaf)

# We can even load the collapsed version as a new tree
#t2 = Tree( t.write(is_leaf_fn=collapsed_leaf) )


#t2.write(format = 1, outfile = '/home/meike/tests/Files/tree/test_collapsed.nwk')
