from pathlib import Path
import os

path = os.getcwd()
p = Path(path)
files_dir = os.path.join(p.parents[1], 'files')


with open (files_dir+'/20012020_streptococcus_database.tsv') as f:
    with open('/home/meike/tests/Files/same_ids.tsv', 'w') as f2:
        for line in f:
            if line.startswith('streptococcus_04796'):
                f2.write(line)