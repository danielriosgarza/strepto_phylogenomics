import os

path  = '/home/meiker/git/genomes/'
files = os.listdir(path)


with open('/home/meiker/git/strepto_phylogenomics/scripts/bash_scripts/rename_strepto.sh','w') as f:
    for i in files:
        if 'streptocuccus' in i:
            name = i.replace('streptocuccus', 'streptococcus')
            f.write('mv ' + path + i + ' ' + path + name + '\n\n')

