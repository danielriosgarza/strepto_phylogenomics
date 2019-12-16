# strepto_phylogenomics

Meike Report: https://docs.google.com/document/d/1oAlTzNKZa6KR_BxNZg9qsh-vp9L4EgcCgheUhfUHTdA/edit?usp=sharing

#Files:

- streptococcus_all_genomes.tsv 

all genomes in the Patric database, including plasmids and low quality genomes. (Total 18936 genome ids)

- strepto_genomes_quality.tsv

all the selected genomes for the project. The parameters for selection:
'''Parameters for genome quality: 
    1. genome quality = good
    2. Completness of more than or equal to 90%
    3. Status not plasmid
    4. If status is nothing ('') than cds of at least 700
    5. Consistencies above or equal to 95%'''
selection is done by the script: genome_quality_check.py. Further, some genomes were included to assure we had at least one per species (in the same script).


- strepto_genome_database.tsv

database identifier, mapping to all information that we have about the strain.


- floricoccus_all_genomes.tsv 

all genomes in the Patric database, including plasmids and low quality genomes. (Total 2 genome ids)


- floricoccus_genome_database.tsv

database identifier, mapping to all information that we have about the strain.

- lactococcus_all_genomes.tsv 

all genomes in the Patric database, including plasmids and low quality genomes. (Total 474 genome ids)

- lactococcus_genomes_quality.tsv

all the selected genomes for the project. The parameters for selection:
'''Parameters for genome quality: 
    1. genome quality = good
    2. Completness of more than or equal to 90%
    3. Status not plasmid
    4. If status is nothing ('') than cds of at least 700
    6. Consistencies above or equal to 95%'''
selection is done by the script: genome_quality_check.py. Further, some genomes were included to assure we had at least one per species (in the same script).


- lactococcus_genome_database.tsv

database identifier, mapping to all information that we have about the strain.
