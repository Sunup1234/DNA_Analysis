
def clean_seq(seq):
    seq = seq.upper()
    valid = {'A','T','G','C'}
    return ''.join(base for base in seq if base in valid)

DNA = input('Enter a sequence: ')
DNA = clean_seq(DNA)

print('DNA:', DNA)

def reverse_complement(seq):
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return ''.join(complement[base] for base in reversed(seq))

rev_comp = reverse_complement(DNA)
print('Reverse complement:' , rev_comp)

def get_codons(seq, frame=0):
    return[seq[i:i+3] for i in range(frame, len(seq)-2,3)]

print("codon : ", get_codons(DNA,0)[:10])
print("reversed codon :", get_codons(rev_comp,0)[:10])

codon_table_full = {
    # Phenylalanine / Leucine
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    # Serine
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    # Tyrosine / Stop
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    # Cysteine / Stop / Tryptophan
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    # Leucine
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    # Proline
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    # Histidine / Glutamine
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    # Arginine
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    # Isoleucine / Methionine (start)
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    # Threonine
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    # Asparagine / Lysine
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    # Serine / Arginine
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    # Valine
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    # Alanine
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    # Aspartate / Glutamate
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    # Glycine
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

codon_name_table = {
    "F": "Phenylalanine",
    "L": "Leucine",
    "S": "Serine",
    "Y": "Tyrosine",
    "*": "Stop",
    "C": "Cysteine",
    "W": "Tryptophan",
    "P": "Proline",
    "H": "Histidine",
    "Q": "Glutamine",
    "R": "Arginine",
    "I": "Isoleucine",
    "M": "Methionine",
    "T": "Threonine",
    "N": "Asparagine",
    "K": "Lysine",
    "V": "Valine",
    "A": "Alanine",
    "D": "Aspartate",
    "E": "Glutamate",
    "G": "Glycine",
    "X": "Unknown"
}

def translate (seq, frame=0):
    codons = get_codons(seq, frame)
    protein = ''.join(codon_table_full.get(c,'X') for c in codons)
    return protein

protein = translate(DNA, frame=0)
protein_rev = translate(reverse_complement(DNA), frame=0)
print(f'Protein: {protein} \nReversed Protein: {protein_rev}')

def longest_codon(seq):
    best_orf = ""
    for strand, nuc in [("forward", seq), ("reverse", reverse_complement(seq))]:
        for frame in range(3):
            protein = translate(nuc, frame)
            orfs = protein.split('_')
            longest = max(orfs, key=len)
            if len(longest) > len(best_orf):
                best_orf = longest
                best = (strand, frame, len(longest))

    return best, best_orf


longest_prot_codon = longest_codon(DNA)[0]
sequence = longest_codon(DNA)[1]
longest_prot_codon_rev = longest_codon(reverse_complement(DNA))
print('The best orf is : ', longest_prot_codon)
print('protein sequence is : ', sequence)





