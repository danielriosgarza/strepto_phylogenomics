# -*- coding: utf-8 -*-
"""
Ete3 test.
"""

import ete3
from ete3 import Tree, TreeStyle, NodeStyle
import matplotlib.pyplot as plt
import matplotlib

def get_colors(data_list, color='Reds'):
    '''
    Needs a dat_list to determine how many different colours are needed. Cmap colour schemes can be provided, default are Reds.
    '''
    colors = []
    cmap = plt.get_cmap(color, len(data_list))    # PiYG
    
    for i in range(cmap.N):
        rgb = cmap(i)[:3] # will return rgba, we take only first 3 so we get rgb
        colors.append(matplotlib.colors.rgb2hex(rgb))
    return colors

def layout(node):
    '''
    Set Nodestyle for tree. 

    '''
    #add face with color per genus (rectangle at outer ring of the tree)
    if node.is_leaf():
        if 'strepto' in node.name:
            genus = ete3.RectFace(width = 10, height = 10, bgcolor = '#d6604d', fgcolor = '#d6604d')
            ete3.faces.add_face_to_node(genus, node, column = 0, position="aligned")
        
        if 'lacto' in node.name:
            genus = ete3.RectFace(width = 10, height = 10, bgcolor = '#2166ac', fgcolor = '#2166ac')
            ete3.faces.add_face_to_node(genus, node, column = 0, position="aligned")
        
        if 'flori' in node.name:
            genus = ete3.RectFace(width = 10, height = 10, bgcolor = '#d1e5f0', fgcolor = '#d1e5f0')
            ete3.faces.add_face_to_node(genus, node, column = 0, position="aligned")


 
t = Tree("(flori_1:1,flori_2:1,(lacto_1:1,lacto_2:1,lacto_3:0.5,lacto_4:0.5,(strepto_1:1,strepto_2:1,strepto_3:0.5,strepto_4:0.25,strepto_5:0.75,strepto_6:0.1):0.5):0.6);")

#colour per species
species = {'s1' :  ['strepto_1', 'strepto_2'],
       's2' : ['strepto_3'],
       'l1' : ['lacto_1'],
       'l2' : ['lacto_2'],
       'f1' : ['flori_1']}     

#leaf_colours --> k = id_ and v = colour and the same species needs the same colour
leaf_colours = {}

streptos = []
lactos =[]
floris = []

#check values (ids) to sort the species per genus
for k,v  in species.items():
    if next((True for v in v if 'strepto' in v), False) == True:
        streptos.append(k)
    if next((True for v in v if 'lacto' in v), False) == True:
        lactos.append(k)
    if next((True for v in v if 'flori' in v), False) == True:
        floris.append(k)



strep_colour = get_colors(streptos)
lacto_colour = get_colors(lactos, color= 'Blues')
flori_colour = get_colors(floris, color = 'Greys')

for i, spec in enumerate(streptos):
    for _id in species[spec]:
        leaf_colours[_id] = strep_colour[i]   

for i, spec in enumerate(lactos):
    for _id in species[spec]:
        leaf_colours[_id] = lacto_colour[i]  
     
for i, spec in enumerate(floris):
    for _id in species[spec]:
        leaf_colours[_id] = flori_colour[i]  

for node in t.traverse():
    node.img_style["size"] = 0 #removes dots at nodes  
    if node.is_leaf():
        color = leaf_colours.get(node.name, None)
        if color:
            name_face = ete3.TextFace(node.name, fgcolor=color, fsize=10)
            node.add_face(name_face, column=0, position='branch-right') 

    
    
    

# for node in t.iter_search_nodes():
#     node.img_style["size"] = 0 #removes dots  nodes  
   
#     if 'strepto' in node.name:
#         node.img_style['fgcolor'] = '#d6604d'

#     if 'lacto' in node.name:
#         node.img_style['fgcolor'] = '#2166ac'

#     if 'flori' in node.name:
#         node.img_style['bgcolor'] = '#d1e5f0'

      


ts = TreeStyle()
ts.mode = 'c'
ts.show_leaf_name = False

ts.layout_fn = layout

t.show(tree_style=ts)

#%%

def get_colors(data_list, color='Reds'):
    '''
    Needs a dat_list to determine how many different colours are needed. Cmap colour schemes can be provided, default are Reds.
    '''
    colors = []
    cmap = plt.get_cmap(color, len(data_list))    # PiYG
    
    for i in range(cmap.N):
        rgb = cmap(i)[:3] # will return rgba, we take only first 3 so we get rgb
        colors.append(matplotlib.colors.rgb2hex(rgb))
    return colors

def layout(node):
    '''
    Set Nodestyle for tree. 

    '''
    #add face with color per genus (circle at outer ring of the tree)
    if node.is_leaf():
        if 'strepto' in node.name:
            genus = ete3.CircleFace(radius = 5, color = '#d6604d')
            ete3.faces.add_face_to_node(genus, node, column = 1, position = "aligned")
        
        if 'lacto' in node.name:
            genus = ete3.CircleFace(radius = 5, color = '#2166ac')
            ete3.faces.add_face_to_node(genus, node, column = 1, position = "aligned",)
        
        if 'flori' in node.name:
            genus = ete3.CircleFace(radius = 5, color = '#636363')
            ete3.faces.add_face_to_node(genus, node, column = 1, position = "aligned")


    
 
t = Tree("(flori_1:1,flori_2:1,(lacto_1:1,lacto_2:1,lacto_3:0.5,lacto_4:0.5,(strepto_1:1,strepto_2:1,strepto_3:0.5,strepto_4:0.25,strepto_5:0.75,strepto_6:0.1):0.5):0.6);")

#colour per species
species = {'s1' :  ['strepto_1', 'strepto_2'],#d1e5f0
       's2' : ['strepto_3'],
       's3' : ['strepto_4'],
       's4' : ['strepto_5', 'strepto_6'],
       'l1' : ['lacto_1'],
       'l2' : ['lacto_2'],
       'l3' : ['lacto_3', 'lacto_4'],
       'f1' : ['flori_1'],
       'f2' : ['flori_2']}     

#leaf_colours --> k = id_ and v = colour and the same species needs the same colour theme

leaf_colours = {}

streptos = []
lactos =[]
floris = []

#check values (ids) to sort the species per genus
for k,v  in species.items():
    if next((True for v in v if 'strepto' in v), False) == True:
        streptos.append(k)
    if next((True for v in v if 'lacto' in v), False) == True:
        lactos.append(k)
    if next((True for v in v if 'flori' in v), False) == True:
        floris.append(k)


#get as many colors (hex codes) as there are species
strep_colour = get_colors(streptos)
lacto_colour = get_colors(lactos, color= 'Blues')
flori_colour = get_colors(floris, color = 'Greys')


#make dict with leafcolors per species
for i, spec in enumerate(streptos):
    for _id in species[spec]:
        leaf_colours[_id] = strep_colour[i]   

for i, spec in enumerate(lactos):
    for _id in species[spec]:
        leaf_colours[_id] = lacto_colour[i]  
     
for i, spec in enumerate(floris):
    for _id in species[spec]:
        leaf_colours[_id] = flori_colour[i]  


check = 0
species_check = []
for node in t.traverse():
    children_nodes = node.children
    species_check = []
    
    #iterate over nodes and check the leaves is the same species is present
    for c in children_nodes:
        if c.is_leaf():
            for k,v in species.items():
                if c in v:
                    species_check.append((c, k))
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
            item.detach()
            check += 1
                


for node in t.traverse():
    node.img_style["size"] = 0 #removes dots at nodes  
    if node.is_leaf():
        color = leaf_colours.get(node.name, None)
        if color:
            node.img_style["bgcolor"] = color
            # strain = ete3.RectFace(width = 10, height = 10, fgcolor = color, bgcolor = color) 
            # node.add_face(strain, column = 0, position="branch-right")

    
    
    

# for node in t.iter_search_nodes():
#     node.img_style["size"] = 0 #removes dots  nodes  
   
#     if 'strepto' in node.name:
#         node.img_style['fgcolor'] = '#d6604d'

#     if 'lacto' in node.name:
#         node.img_style['fgcolor'] = '#2166ac'

#     if 'flori' in node.name:
#         node.img_style['bgcolor'] = '#d1e5f0'

      


ts = TreeStyle()
ts.mode = 'c'

#ts.show_leaf_name = False

# ts.legend.add_face(ete3.CircleFace(5, '#d6604d'), column=0)
# ts.legend.add_face(ete3.TextFace("Streptococcus"), column = 1)

# ts.legend.add_face(ete3.CircleFace(5, '#2166ac'), column=0)
# ts.legend.add_face(ete3.TextFace("Lactococcus"), column=1)

# ts.legend.add_face(ete3.CircleFace(5, '#636363'), column=0)
# ts.legend.add_face(ete3.TextFace("Floricoccus"), column=1)

# ts.legend_position = 1
# ts.title.add_face(ete3.TextFace("Genus"), column = 1)

#ts.layout_fn = layout

t.show(tree_style=ts)