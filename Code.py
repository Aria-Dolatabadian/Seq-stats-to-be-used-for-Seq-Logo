filename = "seq.txt"

with open(filename, "r") as file:
    lines = file.readlines()

genes = []
current_gene = None

# Parse the genes from the file
for line in lines:
    if line.startswith(">"):
        if current_gene is not None:
            genes.append(current_gene)
        current_gene = {"name": line.strip()[1:], "sequence": ""}
    else:
        current_gene["sequence"] += line.strip()

if current_gene is not None:
    genes.append(current_gene)

# Calculate the nucleotide counts and percentages for each gene
for gene in genes:
    counts = {"A": 0, "T": 0, "C": 0, "G": 0}
    length = len(gene["sequence"])
    for nt in gene["sequence"]:
        if nt in counts:
            counts[nt] += 1
    percentages = {nt: count / length for nt, count in counts.items()}
    # percentages = {nt: count / length * 100 for nt, count in counts.items()}
    print(f"{gene['name']}: A={counts['A']} ({percentages['A']:.2f}), T={counts['T']} ({percentages['T']:.2f}), C={counts['C']} ({percentages['C']:.2f}), G={counts['G']} ({percentages['G']:.2f})")
    # print(f"{gene['name']}: A={counts['A']} ({percentages['A']:.2f}%), T={counts['T']} ({percentages['T']:.2f}%), C={counts['C']} ({percentages['C']:.2f}%), G={counts['G']} ({percentages['G']:.2f}%)")

import csv

# Define function to calculate nucleotide counts and percentages
def calc_nucleotide_percentages(seq):
    nucleotides = ['A', 'T', 'C', 'G']
    counts = {n: seq.count(n) for n in nucleotides}
    total_count = sum(counts.values())
    percentages = {n: count/total_count for n, count in counts.items()}
    # percentages = {n: count/total_count * 100 for n, count in counts.items()}
    return counts, percentages

# Open input file and read gene sequences
with open('seq.txt', 'r') as f:
    gene_sequences = {}
    current_gene = None
    for line in f:
        if line.startswith('>'):
            current_gene = line.strip()[1:]
            gene_sequences[current_gene] = ''
        else:
            gene_sequences[current_gene] += line.strip()

# Calculate nucleotide counts and percentages for each gene
results = {}
for gene, seq in gene_sequences.items():
    counts, percentages = calc_nucleotide_percentages(seq)
    results[gene] = percentages

# Export results to CSV file
with open('genes_stats.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Gene', 'A', 'T', 'C', 'G'])
    # writer.writerow(['Gene', 'A%', 'T%', 'C%', 'G%'])
    for gene, percentages in results.items():
        writer.writerow([gene, percentages['A'], percentages['T'], percentages['C'], percentages['G']])

# Transpose CSV table
with open('genes_stats.csv', 'r') as f:
    rows = zip(*csv.reader(f))
    with open('genes_stats_transposed.csv', 'w', newline='') as f_transposed:
        writer = csv.writer(f_transposed)
        writer.writerows(rows)
