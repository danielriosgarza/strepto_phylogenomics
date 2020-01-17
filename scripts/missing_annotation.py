

import os


def missing_elements(L):
    start, end = L[0], L[-1]
    return sorted(set(range(start, end + 1)).difference(L))

processed_genomes = os.listdir('/home/meiker/git/prokka_annotation/')

missing = []
all_genomes=[]

with open('/home/meiker/git/strepto_phylogenomics/files/lactococcus_genomes_quality.tsv') as f:
    f.readline()
    for line in f:
        a=line.strip().split('\t')
        all_genomes.append(a[0])
        if a[0] not in processed_genomes:
            missing.append(a[0])
print(len(all_genomes))

print(len(set(all_genomes)))

print(missing)
