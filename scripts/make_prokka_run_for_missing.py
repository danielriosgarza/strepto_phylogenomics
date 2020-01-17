import os

processed_genomes = os.listdir('/home/meiker/git/prokka_annotation/')


with open('/home/meiker/git/strepto_phylogenomics/files/streptococcus_genomes_quality.tsv') as f:
    f.readline()
    with open('/home/meiker/git/strepto_phylogenomics/scripts/bash_scripts/prokka_annotate_missing_strepto.sh', 'w') as f2:
        for line in f:
            a = line.strip().split('\t')
            if a[0] not in processed_genomes:
                f2.write('prokka /home/meiker/git/genomes/'+a[0]+'.fna --outdir /home/meiker/git/prokka_annotation/' + a[0] + ' --prefix ' +a[0]+' --genus Streptococcus' + '\n')
