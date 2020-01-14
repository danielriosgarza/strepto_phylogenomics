blastp -query blastquery/floricoccus_00001.fasta  -db blastdb/goodProteins.fasta  -seg yes  -dbsize 100000000  -evalue 1e-5  -outfmt 6 -num_threads 8 -out blastres/floricoccus_00001.tab
blastp -query blastquery/floricoccus_00002.fasta  -db blastdb/goodProteins.fasta  -seg yes  -dbsize 100000000  -evalue 1e-5  -outfmt 6 -num_threads 8 -out blastres/floricoccus_00002.tab
