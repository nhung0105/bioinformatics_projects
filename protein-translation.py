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
