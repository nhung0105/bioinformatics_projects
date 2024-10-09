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
print("DNA Sequence: ", dna_sequence)

#count nucleotides

def count_nucleotides(sequence):
    nucleotide_count = {
        "A": sequence.count("A"),
        "T": sequence.count("T"),
        "C": sequence.count("C"),
        "G": sequence.count("G"),
    }
    return nucleotide_count

nucleotide_count = count_nucleotides(dna_sequence)
print("Nucleotide count: ", nucleotide_count)

#Calculate G-C content
def calculate_gc(sequence):
    gc_count = sequence.count("G") + sequence.count("C")
    gc_content = (gc_count / len(sequence)) * 100
    return gc_content

gc_content = calculate_gc(dna_sequence)
print("GC content: {:.2f}%".format(gc_content)) 

#Reverse complement

def get_reverse_comp(sequence):
    complement = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C"
    }
    
    reverse_comp_list = []
    
    for nucleotide in reversed(sequence):
        reverse_comp_list.append(complement[nucleotide])
        # appened the complement of the nucleotide to the list

    reverse_complement = ''.join(reverse_comp_list)
    return reverse_complement

sequence = dna_sequence
reverse_complement = get_reverse_comp(sequence)
print("Reverse complement: ", reverse_complement)

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
