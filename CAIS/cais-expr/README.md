# CAIS-expression relationships

All python and R scripts used to calculate CAIS and CAI of genes for mouse and human.

Steps:
0. Download the following from ensembl biomart
	Dataset Human or mouse genes
	Attributes	Gene stable ID, Gene stable ID version, Protein stable ID, Protein stable ID version, Coding sequence
	Results Export all results to fasta file

1. Download the following from PaxDB
	Organism - Whole organism (Integrated)

2. Use the first couple blocks of organism_cais_expr.Rmd to create organism_genes_100codons.fasta

OR UNZIP GENES.zip

3. Concatenate all protein-coding genes over 100 codons in length
awk '/^>/ {next} {seq = seq $0} END {print ">organism\n" seq}' organism_genes_100codons.fasta > organism_concat_100codons.fasta

4. Getet GC content
grep -v '^>' organism_concat_100codons.fasta | tr -d '\n' | awk '{
    seq = toupper($0);
    g_count = gsub(/G/, "", seq);
    c_count = gsub(/C/, "", seq);
    gc = g_count + c_count;
    total = length($0);
    print "GC content:", gc / total * 100 "%";
}'

which gives GC content:
Mouse all genes: 0.5190
Human all genes: 0.5002
Mouse >100 codons: 0.5102
Human > 100 codons: 0.5108

5. Input these values in CAIS_organismAAGC_ofgenes.py, CAI_fixed_organism_ofgenes.py, and CAI_fixed_organismAAGC_global.py.

6. Run:
	python CAIS_organismAAGC_ofgenes.py > CAIS_ofgenes_organism.txt
to calculate CAIS for each gene 

7. Run CAI_fixed_organismAAGC_global.py to calculate RSCUmax across all 100 genes

which gives RSCUmax: 
Mouse all genes: 1.2690
Human all genes: 1.3548
Mouse >100 codons: 0.3401
Human >100 codons: 0.3382

6. Input these values in CAI_fixed_organism_ofgenes.py and run:
	python CAI_organismAAGC_ofgenes.py > CAI_ofgenes_organism.txt










