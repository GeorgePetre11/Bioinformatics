import matplotlib.pyplot as plt
from collections import Counter

# Genetic code (same as before)
genetic_code = {
    'UUU':'Phe', 'UUC':'Phe', 'UUA':'Leu', 'UUG':'Leu',
    'CUU':'Leu', 'CUC':'Leu', 'CUA':'Leu', 'CUG':'Leu',
    'AUU':'Ile', 'AUC':'Ile', 'AUA':'Ile', 'AUG':'Met',
    'GUU':'Val', 'GUC':'Val', 'GUA':'Val', 'GUG':'Val',
    'UCU':'Ser', 'UCC':'Ser', 'UCA':'Ser', 'UCG':'Ser',
    'CCU':'Pro', 'CCC':'Pro', 'CCA':'Pro', 'CCG':'Pro',
    'ACU':'Thr', 'ACC':'Thr', 'ACA':'Thr', 'ACG':'Thr',
    'GCU':'Ala', 'GCC':'Ala', 'GCA':'Ala', 'GCG':'Ala',
    'UAU':'Tyr', 'UAC':'Tyr', 'UAA':'Stop', 'UAG':'Stop',
    'CAU':'His', 'CAC':'His', 'CAA':'Gln', 'CAG':'Gln',
    'AAU':'Asn', 'AAC':'Asn', 'AAA':'Lys', 'AAG':'Lys',
    'GAU':'Asp', 'GAC':'Asp', 'GAA':'Glu', 'GAG':'Glu',
    'UGU':'Cys', 'UGC':'Cys', 'UGA':'Stop', 'UGG':'Trp',
    'CGU':'Arg', 'CGC':'Arg', 'CGA':'Arg', 'CGG':'Arg',
    'AGU':'Ser', 'AGC':'Ser', 'AGA':'Arg', 'AGG':'Arg',
    'GGU':'Gly', 'GGC':'Gly', 'GGA':'Gly', 'GGG':'Gly'
}

def read_fasta(filename):
    """Reads a FASTA file and returns the sequence as a single string."""
    with open(filename, 'r') as f:
        lines = f.readlines()
    seq = ''.join(line.strip() for line in lines if not line.startswith('>'))
    seq = seq.upper().replace('T', 'U')  # convert DNA to RNA
    return seq

def count_codons(sequence):
    """Count codons in a given RNA sequence."""
    codons = [sequence[i:i+3] for i in range(0, len(sequence)-2, 3) if len(sequence[i:i+3]) == 3]
    return Counter(codons)

def codon_to_amino_counts(codon_counts):
    """Convert codon frequency counts into amino acid frequency counts."""
    aa_counts = Counter()
    for codon, count in codon_counts.items():
        aa = genetic_code.get(codon, None)
        if aa and aa != 'Stop':
            aa_counts[aa] += count
    return aa_counts

def plot_top_codons(codon_counts, title):
    """Plot top 10 codons."""
    top = codon_counts.most_common(10)
    codons, counts = zip(*top)
    plt.figure(figsize=(8,5))
    plt.bar(codons, counts, color='teal')
    plt.title(title)
    plt.xlabel('Codon')
    plt.ylabel('Frequency')
    plt.tight_layout()
    plt.show()

# ---------------- MAIN SCRIPT ----------------

covid_seq = read_fasta("covid.fasta")
flu_seq = read_fasta("influenza.fasta")

covid_counts = count_codons(covid_seq)
flu_counts = count_codons(flu_seq)

# (a) Chart for COVID-19
plot_top_codons(covid_counts, "Top 10 Most Frequent Codons - COVID-19")

# (b) Chart for Influenza
plot_top_codons(flu_counts, "Top 10 Most Frequent Codons - Influenza")

# (c) Compare top codons
common_codons = set([c for c, _ in covid_counts.most_common(10)]) & set([c for c, _ in flu_counts.most_common(10)])
print("Common top codons between COVID-19 and Influenza:", ", ".join(common_codons) if common_codons else "None")

# (d) Top 3 amino acids for each genome
covid_aa = codon_to_amino_counts(covid_counts)
flu_aa = codon_to_amino_counts(flu_counts)

print("\nTop 3 amino acids - COVID-19:")
for aa, count in covid_aa.most_common(3):
    print(f"{aa}: {count}")

print("\nTop 3 amino acids - Influenza:")
for aa, count in flu_aa.most_common(3):
    print(f"{aa}: {count}")