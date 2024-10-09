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
