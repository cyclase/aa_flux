import operator
import numpy as np
from Bio import SeqIO
from functools import reduce


# paths
genomes_dir = '/Users/karinazile/PhD/TRGs/data/2018_03_15_five_drosophilas_refseq'
results_dir = '/Users/karinazile/PhD/TRGs/results/2019_07_03_cai_analysis'
sequence_file = '/Users/karinazile/PhD/TRGs/results/2019_07_03_cai_analysis/trgf_relevant_elements.fna'
species_filename = dict([('DROME','GCF_000001215.4_Release_6_plus_ISO1_MT'), ('DROSI','GCF_000754195.2_ASM75419v2'), ('DROSE','GCF_000005215.3_dsec_caf1'), ('DROYA','GCF_000005975.2_dyak_caf1'), ('DROER','GCF_000005135.1_dere_caf1')])

list_of_species = ['DROME', 'DROSI', 'DROSE', 'DROYA', 'DROER']
list_of_codons = ['AAA', 'AAT', 'AAG', 'AAC', 'ATA', 'ATT', 'ATG', 'ATC', 'AGA', 'AGT', 'AGG', 'AGC', 'ACA', 'ACT', 'ACG', 'ACC', 
				  'TAA', 'TAT', 'TAG', 'TAC', 'TTA', 'TTT', 'TTG', 'TTC', 'TGA', 'TGT', 'TGG', 'TGC', 'TCA', 'TCT', 'TCG', 'TCC', 
				  'GAA', 'GAT', 'GAG', 'GAC', 'GTA', 'GTT', 'GTG', 'GTC', 'GGA', 'GGT', 'GGG', 'GGC', 'GCA', 'GCT', 'GCG', 'GCC',
				  'CAA', 'CAT', 'CAG', 'CAC', 'CTA', 'CTT', 'CTG', 'CTC', 'CGA', 'CGT', 'CGG', 'CGC', 'CCA', 'CCT', 'CCG', 'CCC']

aa_codon = dict([('A',['GCA', 'GCT', 'GCG', 'GCC']),
	             ('R',['AGG','AGA', 'CGA', 'CGT', 'CGG', 'CGC']),
	             ('N',['AAT', 'AAC']),
	             ('D',['GAT', 'GAC']),
	             ('C',['TGT', 'TGC']),
	             ('Q',['CAA', 'CAG']),
	             ('E',['GAA', 'GAG']),
	             ('G',['GGA', 'GGT', 'GGG', 'GGC']),
	             ('H',['CAT', 'CAC']),
	             ('I',['ATA', 'ATT', 'ATC']),
	             ('L',['TTA', 'TTG', 'CTA', 'CTT', 'CTG', 'CTC']),
	             ('K',['AAA', 'AAG']),
	             ('M',['ATG']),
	             ('F',['TTT', 'TTC']),
	             ('P',['CCA', 'CCT', 'CCG', 'CCC']),
	             ('S',['TCA', 'TCT', 'TCG', 'TCC', 'AGT', 'AGC']),
	             ('T',['ACA', 'ACT', 'ACG', 'ACC']),
	             ('W',['TGG']),
	             ('Y',['TAT', 'TAC']),
	             ('V',['GTA', 'GTT', 'GTG', 'GTC']),
	             ('*',['TAA', 'TAG', 'TGA'])])




### count how often each codon is used in each species (only using one longest isoform per gene)
### only count sequences whose length is multiple of 3 and that do not contain "N"
raw_codon_counts = np.zeros((len(list_of_species), len(list_of_codons)))

for index, species in enumerate(list_of_species):
	print(species)
	# collect all isoform sequences from fasta file
	fasta_dict = SeqIO.index('%s/%s_cds_from_genomic.fna'%(genomes_dir, species_filename[species]), 'fasta')
	locus_seqs = {}
	for k, v in fasta_dict.items():
		try:
			locus = v.description.split('gene=')[1].split(']')[0]
			sequence = str(v.seq)
			if len(sequence) % 3 == 0:
				if 'N' not in sequence:
					locus_seqs.setdefault(locus, []).append(sequence)
		except:
			pass
	print(len(locus_seqs))
	# for every locus select the longest isoform and count codons it contains
	for k, v in locus_seqs.items():
		list_of_len = [len(x) for x in v]
		isoform = list_of_len.index(max(list_of_len))
		sequence = v[isoform]
		for i in range(0, len(sequence), 3):
			raw_codon_counts[index, list_of_codons.index(sequence[i:i+3])] += 1
np.savetxt('%s/raw_codon_counts.tsv'%(results_dir), raw_codon_counts, delimiter='\t', fmt ='%.0f')




### calculate Relative Synonymous Codon Usage (RSCU) value and relative adaptedness value for every codon 
### RSCU = (number of occurances of this codon * number of codons for this aa) / number of occurances of this aa
### If all codons are used equally frequently, then all values will be 1.
### relative adaptedness = RSCU / maximum RSCU value for that amino acid
rscu = np.zeros((len(list_of_species), len(list_of_codons)))
rel_adapt = np.zeros((len(list_of_species), len(list_of_codons)))

for index, species in enumerate(list_of_species):
	for k, v in aa_codon.items():
		frequencies = [raw_codon_counts[index, list_of_codons.index(x)] for x in v]
		if sum(frequencies) == 0:
			print('%s does not occur in %s'%(k, species))
		else:
			multiplier = len(frequencies) / sum(frequencies)
			max_value = max(frequencies) * multiplier
			for codon in v:
				rscu[index, list_of_codons.index(codon)] = raw_codon_counts[index, list_of_codons.index(codon)] * multiplier
				rel_adapt[index, list_of_codons.index(codon)] = rscu[index, list_of_codons.index(codon)] / max_value

np.savetxt('%s/RSCUs.tsv'%(results_dir), rscu, delimiter='\t')
np.savetxt('%s/relative_adaptedness.tsv'%(results_dir), rel_adapt, delimiter='\t')

if np.count_nonzero(rel_adapt) != rel_adapt.shape[0]*rel_adapt.shape[1]:
	print('Not all codons have a non-zero relative adaptedness value!')



### CAI, or Codon Adaptation Index, measures the degree of to which a gene or gene segment 
### exhibits the same codon bias as the rest of the species
### CAI = the geometric mean of the relative adaptedness values of all codons in a sequence
### CAI can not be calculated for sequences containing only M (ATG) and W (TGG)
### Using CAI gene segments are either tentatively classified as ancestral or novel
### Higher CAI ==> ancestral (gene more closely resembles overall genome codon usage patterns)
### Lower CAI ==> novel (gene less closely resembles overall genome codon usage patterns)

interest_species = 'DROSI'
interest_sequence = 'ATGCCTCCAGCTGGGCTGGCCCCTCGTCCTTCGCTTCATATCCTGTCAGGAGCTGACTTTTTCCTGTTGCCGCTGCTTCCTTCCTTCCTTTCTTCTGGCTCTTATCTTCGTTTAACTTTCTACCATGTTATGGCCCGTTATCTGTCTGCCCGTGTACGTGTCCTTCGTGTCGGTGTGTGTTCGGTTTCTATGTGTGTGTGTGTGTATGTGTGTGCTTCCTCGTCTCTGTTTGCATTTTTGAGCAAATTAAGTGAAATTACATTTTTCTATGCGATGCTCTTAATGTTATTTGCATACGTATCCTGCAGGCTACAATATTTACAAATGAATGCAGTCCGGTGCCGTTCGCTTGGTCCTTGGTCCTCGGTCCTTGGGGGCTTGGTCCTTTCGAGAATCAGACAAGCGGCTGTTGAACCAACCACCTTAACTGCCCAACTTCCCAGCAGACCAACCGTCCAGCCGAACTAA'
# if len(interest_sequence) % 3 == 0:
# 	rel_adapt_values = [rel_adapt[list_of_species.index(interest_species), list_of_codons.index(interest_sequence[x:x+3])] for x in range(0, len(interest_sequence), 3)]
# 	product_of_rel_adapt = reduce(operator.mul, rel_adapt_values, 1)
# 	cai = product_of_rel_adapt**(3/len(interest_sequence))


if len(interest_sequence) % 3 == 0:
	product_of_rel_adapt = 1
	num_of_inf_codons = 0
	for i in range(0, len(interest_sequence), 3):
		codon = interest_sequence[i:i+3]
		if codon not in ['ATG', 'TGG']:
			product_of_rel_adapt = product_of_rel_adapt * rel_adapt[list_of_species.index(interest_species), list_of_codons.index(codon)]
			num_of_inf_codons += 1
	cai = product_of_rel_adapt**(1/num_of_inf_codons)