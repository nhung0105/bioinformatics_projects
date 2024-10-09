#read fasta file
def read_fasta(filename):
    with open(filename,"r") as file:
        lines = file.readlines()

        sequence = []
        for line in lines:
            if not line.startswith('>'):
                sequence.append(line.strip())

        full_sequence = ''.join(sequence)

    return full_sequence

filename = "rabbit.fna"
dna_sequence = read_fasta(filename)

#Create a codon dictionary
genetic_code = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}
def translate_dna_to_protein(dna_sequence):
    protein_sequence = []
    
    #group 3 nucleotides together
    for i in range(0, len(dna_sequence) - 2, 3):
        codon = dna_sequence[i:i+3]
        if codon in genetic_code:
            amino_acid = genetic_code[codon]
        if amino_acid == "_" : #stop codon
             break
        protein_sequence.append(amino_acid)
    return ''.join(protein_sequence)

protein_sequence = translate_dna_to_protein(dna_sequence)
print("Protein Sequence: ", protein_sequence)


#Plots

def count_nucleotides(dna_sequence):
    counts = {
        "A": dna_sequence.count('A'),
        "T": dna_sequence.count('T'),
        "C": dna_sequence.count('C'),
        "G": dna_sequence.count('G')
    }
    return counts

import matplotlib.pyplot as plt

def plot_nucleotide_counts(counts):
    nucleotides = list(counts.keys()) #retrieve all keys from dictionary
    values = list(counts.values())

    plt.figure(figsize = (8, 6))
    plt.bar(nucleotides, values, color=['green', 'blue', 'orange', 'red'])
    plt.title("Nucleotide Counts")
    plt. xlabel("Nucleotides")
    plt.ylabel("Count")
    plt.show()

# Plotting GC content
def plot_gc_content(counts):
    gc_count = counts['G'] + counts['C']
    at_count = counts['A'] + counts['T']
    
    plt.figure(figsize=(8, 6))
    plt.pie([gc_count, at_count], labels=["GC Content", "AT Content"], autopct='%1.1f%%', colors=['orange', 'blue'])
    plt.title("GC Content vs AT Content")
    plt.show()

# Example usage
filename = "rabbit.fna"
dna_sequence = read_fasta(filename)

# Count nucleotides
nucleotide_counts = count_nucleotides(dna_sequence)

# Plot nucleotide counts
plot_nucleotide_counts(nucleotide_counts)

# Plot GC content
plot_gc_content(nucleotide_counts)
