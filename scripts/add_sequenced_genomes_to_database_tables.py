#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 14:03:12 2020

@author: meike
"""

'''
Add 24 genomes to metatable and database_id table.
'''

import os
from pathlib import Path


path = os.getcwd()
p = Path(path)


#remove the two genomes that had no fna file (#11895, #11939)
with open (os.path.join(p.parents[0], 'files' , '03032020_streptococcus_database.tsv')) as f:
    with open (os.path.join(p.parents[0], 'files' , '03032020_streptococcus_database_corrected.tsv') , 'w') as f2:
           for line in f:
               if line.startswith('streptococcus_11895') or line.startswith('streptococcus_11939'):
                   pass
               else:
                   f2.write(line)

strepto_id_file = os.path.join(p.parents[0], 'files' , '03032020_streptococcus_patric_id_with_database_id.tsv')

#Add the 24 sequenced genomes to the end of the database to patric id table. Adds column species for the sequenced ones
with open(strepto_id_file) as f:
    with open (os.path.join(p.parents[0], 'files', 'sequenced_species_name_list')) as f2:
        with open(os.path.join(p.parents[0], 'files' , '03032020_streptococcus_patric_id_with_database_id_including_sequenced_genomes.tsv'), 'w') as f3:
            headers = f.readline().strip().split('\t')
            headers.append('species') #for mapping of the sequenced genomes with the db_id
            f3.write('\t'.join(headers)+'\n')
            for line in f:
                a = line.strip().split('\t')
                #remove the ids lacking the fna file from the databse mapping
                if a[0] == 'streptococcus_11895' or a[0] == 'streptococcus_11939':
                    pass
                else:
                    f3.write('\t'.join(a)+'\n')
            last_id = a[0]    
            last_id = last_id.split('_') #get last db_id to generate the following of the sequenced ones
            species = []
            for line in f2:
                line= line.strip()
                species.append(line)
            new_db_ids = []
            for i in range(len(species)):
                db_id = last_id[0]+'_'+ str(int(last_id[1])+i+1)
                new_db_ids.append(db_id)
                f3.write('\t'.join([db_id, '', species[i], '\n']))
                 
#Add genomes to database table
with open(os.path.join(p.parents[0], 'files' , '03032020_streptococcus_database_corrected.tsv')) as f:
    with open(os.path.join(p.parents[0], 'files' , '03032020_streptococcus_database_final.tsv'), 'w') as f2:
        headers = f.readline().strip().split('\t')
        headers_inds = {i : name for name, i in enumerate(headers)}
        
        db_index = headers_inds['database_id']
        species_index = headers_inds['species']
        
        f2.write('\t'.join(headers)+'\n')
        for line in f:
            a = line.strip().split('\t')
            number_columns = len(a)
            f2.write(line)
            
        #make new lines to add db_id and species to right columns
        for i in range(len(species)):
            newline = ['']*number_columns
            newline[db_index] = new_db_ids[i]
            newline[species_index] = species[i]
            f2.write('\t'.join(newline)+'\n')
            
