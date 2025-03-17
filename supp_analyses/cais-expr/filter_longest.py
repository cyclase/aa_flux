from Bio import SeqIO

def filter_longest_sequences(input_fasta, output_fasta):
    # Dictionary to store the longest sequence for each gene
    gene_to_longest_seq = {}

    for record in SeqIO.parse(input_fasta, "fasta"):
        gene_name = record.id.split('|')[0]  # Extract gene name
        seq_length = len(record.seq)

        # Update the record if it's the longest sequence for the gene
        if gene_name not in gene_to_longest_seq or seq_length > len(gene_to_longest_seq[gene_name].seq):
            gene_to_longest_seq[gene_name] = record

    # Write the filtered sequences to a new FASTA file
    with open(output_fasta, "w") as output_file:
        SeqIO.write(gene_to_longest_seq.values(), output_file, "fasta")

# Example usage
input_fasta = "mouse_genes.fasta"
output_fasta = "mouse_genes_onlylongest.fasta"
filter_longest_sequences(input_fasta, output_fasta)
