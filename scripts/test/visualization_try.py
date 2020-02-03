#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 09:37:18 2020

@author: meike
"""

'''
Trying different methods to visualize data.
'''

import os
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import matplotlib.gridspec as gridspec
import heapq

def get_info(file, field):
    '''
    Give input file and column of interest. Returns list with all items of that field.
    '''
    item_l = []
    with open(file) as f:
        headers = f.readline().strip().split('\t')
        headers_inds = {name: i for i, name in enumerate(headers)}
            
        for line in f:
            a = line.strip().split('\t')
            for i, item in enumerate(a):
                if headers_inds[field] == i and item != '':
                    if item.isdigit()==1:
                        item_l.append(int(item))
                    else:
                        item_l.append(str(item))
    return item_l

def count(item_l):
    '''
    Counts items in a list, returns lists containg data and labels.
    '''    
    count = {}
    for item in item_l:
        if item not in count:
            count[item] = 1
        else:
            count[item] += 1
    return count
    
    
path = os.getcwd()
p = Path(path)


def get_df_count(l_for_counts, species):
    '''
    Giv elist with non-numerical data. counts them and give df with unique values and 
    corresponding count. Includes a column with the species.
    '''
    platforms =[]
    count =[]
    for item in l_for_counts:
        if item not in platforms:
            platforms.append(item)
            count.append(1)
        else:
            count[platforms.index(item)] += 1
    for i, number in enumerate(count):
        count[i] = round(number/sum(count)*100)
    df = pd.DataFrame()
    df['Platform'] = platforms
    df['Frequency'] = count
    df['Species'] = species
    return df

#Get sequencing platforms of all species
lac_splatform = get_info(os.path.join(p.parents[0], 'files', '06012020_lactococcus_database.tsv'), 'sequencing_platform')
lac_splatform_count = count(lac_splatform)

strepto_splatform = get_info(os.path.join(p.parents[0], 'files', '20012020_streptococcus_database.tsv'), 'sequencing_platform')
strepto_splatform_count = count(strepto_splatform)

flori_splatform = get_info(os.path.join(p.parents[0], 'files', '06012020_floricoccus_database.tsv'), 'sequencing_platform')
flori_splatform_count = count(flori_splatform)


fig = plt.figure(constrained_layout = 1)
grid = fig.add_gridspec(2,2)

ax1 = fig.add_subplot(grid[:,0])
ax2 = fig.add_subplot(grid[0,1])
ax3 = fig.add_subplot(grid[1,1])




#####Look how you can get counts for freq and col names for legend#####
df = pd.DataFrame({'Lactococcus': pd.Series(lac_splatform),
                   'Streptococcus': pd.Series(strepto_splatform),
                   'Floricoccus' : pd.Series(flori_splatform)})
lac_pcount = df['Lactococcus'].value_counts()
flori_pcount = df['Floricoccus'].value_counts()
strepto_pcount = df['Streptococcus'].value_counts()




lactodf = get_df_count(lac_splatform, 'Lactococcus')
floridf = get_df_count(flori_splatform, 'Floricoccus')
streptodf = get_df_count(strepto_splatform, 'Streptococcus')

df = streptodf.append(lactodf).append(floridf)

df2 = df.groupby(['Platform', 'Species'])['Platform'].count().unstack('Species').fillna(0)
df2[['Species']].plot(kind='bar', stacked=1)

df['Platform'].plot(kind='bar', stacked=1)


df1 = lactodf.sort_values('Frequency', ascending=0).head(10).reset_index(drop=1)
df2 = streptodf.sort_values('Frequency', ascending=0).head(10).reset_index(drop=1)

df = df1.append(df2).append(floridf)





fig, ax = plt.subplots()

ax = sns.barplot(floridf['Platform'], floridf['Frequency'], label = "Floriococcus", palette = 'Reds_r') 
ax = sns.barplot(df1['Platform'], df1['Frequency'], label = "Lactococcus", palette = 'Blues_r') 
ax = sns.barplot(df2['Platform'], df2['Frequency'], label = "Streptococcus", palette = 'Greens_r') 
plt.xticks(rotation='vertical')
plt.title("Frequency distribution of used sequencing platforms", fontsize= 16)
plt.legend(loc= 'upper right')


df1 = pd.DataFrame(lac_splatform, columns=['Platform'])
df2 = pd.DataFrame()
df2['Platform'] = df1['Platform'].unique()
df2['Count'] = df1['Platform'].value_counts()
df2['Species'] = "Lactococcus"

# df2 = pd.DataFrame(flori_splatform, columns=['Platform'])
# df2["Species"] = 'Floricoccus'
# df = df1.append(df2)

fig, ax = plt.subplots()

ax = sns.boxplot(x=df['Platform'].unique(), y=df['Platform'].value_counts())
plt.xticks(rotation='vertical')


#ax = sns.barplot(strepto_pcount.index, strepto_pcount.values, alpha=0.75, palette= 'Oranges')
ax = sns.barplot(flori_pcount.index, flori_pcount.values, alpha=0.75, palette= 'Reds')
ax= sns.barplot(lac_pcount.index, lac_pcount.values, alpha=0.75, palette= 'Blues')
plt.title("Frequency distribution of used sequencing platforms", fontsize= 16)
plt.ylabel('Frequency')
plt.xlabel("Sequencing Platform")
plt.xticks(rotation='vertical')
plt.legend(df.columns, loc= 'upper right')

###################################

theme = plt.get_cmap('RdYlBu')
colors = [theme(1. * i / len(df)) for i in range(len(df))]

fig = plt.plot(figsize =(8,8))
platform_count = df['Platform'].value_counts()
sns.set(style="whitegrid")
sns.barplot(platform_count.index, platform_count.values, alpha=0.9, palette= 'Spectral')
plt.title("Frequency distribution of used sequencing platforms", fontsize= 16)
plt.ylabel('Frequency')
plt.xlabel("Platform")
plt.xticks(rotation='vertical')
plt.tight_layout()



plt.table(df["Sequencing Platform"])





plt.hist(fd.values, histtype='bar', bins=10)


ax = sns.countplot(data=pd.melt(fd), x='value', hue='variable')
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)


sns.catplot(x="Platform", hue="Species",  data=lac_splatform)






#Make grid of plots (4 total) that can be filled
fig, axs = plt.subplots(2,2, figsize= (10,8)) #2x2 grid with determined figure size

#fill subplots with histograms
ax1 = sns.distplot(lac_splatform, color = '#91cf60', label = 'Lactococcus', ax=axs[0,0])
ax1 = sns.distplot(strepto_splatform, color = '#3288bd', label = 'Streptococcus', ax=axs[0,0])
# ax1 = sns.distplot(all_gs['Floricoccus'], color = '#b2182b', label = 'Floricoccus', ax=axs[0,0])
ax2 = sns.distplot(lac_splatform, color = '#91cf60', label = 'Lactococcus', ax=axs[0,1], rug=1)
ax3 = sns.distplot(strepto_splatform, color = '#3288bd', label = 'Streptococcus', ax=axs[1,0], rug=1)
ax4 = sns.distplot(flori_splatform,bins=20, color = '#b2182b', label = 'Floricoccus', ax=axs[1,1], kde=0, norm_hist=0, rug=1)
ax1.legend(loc='upper right')

#Titles
fig.suptitle('Sequencing Platforms', fontsize=16, y=0.95)
ax1.set_title('Lactococcus and Streptococcus',  weight = 'italic')
ax2.set_title('Lactococcus')
ax3.set_title('Streptococcus')
ax4.set_title('Floricoccus')

#Set x and y labels for all plots
# for ax in axs.flat:
    # ax.set(xlabel='Genome size (in Mb)', ylabel='Density')

ax1.set(xlabel='Genome size (in Mb)', ylabel='Density')
ax2.set(xlabel='Genome size (in Mb)', ylabel='Density')
ax3.set(xlabel='Genome size (in Mb)', ylabel='Density')
ax4.set(xlabel='Genome size (in Mb)', ylabel='Frequency')

