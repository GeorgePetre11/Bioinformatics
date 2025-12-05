GENETIC_CODE = {
    'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
    'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
    'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'Stop', 'UAG': 'Stop',
    'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'Stop', 'UGG': 'Trp',
    'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
    'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
    'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
    'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
    'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
    'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
    'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
    'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
    'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly',
}


def dna_to_rna(dna_sequence):
    return dna_sequence.replace('T', 'U')


def translate_to_amino_acids(rna_sequence):
    amino_acids = []
    
    for i in range(0, len(rna_sequence) - 2, 3):
        codon = rna_sequence[i:i+3]
        
        if len(codon) != 3 or not all(base in 'AUGC' for base in codon):
            continue
        
        if codon in GENETIC_CODE:
            amino_acid = GENETIC_CODE[codon]
            amino_acids.append(amino_acid)
            
            if amino_acid == 'Stop':
                break
    
    return amino_acids, ''.join(amino_acids)


def convert_gene_to_protein(dna_sequence):
    dna_sequence = dna_sequence.upper().strip()
    
    if not dna_sequence:
        return {"error": "Empty sequence"}
    
    if not all(base in 'ATGC' for base in dna_sequence):
        return {"error": "Invalid DNA sequence - contains non-ATGC characters"}
    
    if len(dna_sequence) % 3 != 0:
        return {"warning": f"Sequence length ({len(dna_sequence)}) is not a multiple of 3"}
    
    rna_sequence = dna_to_rna(dna_sequence)
    
    amino_acid_list, amino_acid_string = translate_to_amino_acids(rna_sequence)
    
    return {
        "dna_sequence": dna_sequence,
        "rna_sequence": rna_sequence,
        "amino_acids": amino_acid_list,
        "protein_sequence": amino_acid_string,
        "length": len(amino_acid_list)
    }


def main():
    print("=" * 60)
    print("DNA to Amino Acid Sequence Converter")
    print("=" * 60)
    
    while True:
        print("\nOptions:")
        print("1. Convert DNA sequence to amino acids")
        print("2. Exit")
        
        choice = input("\nEnter your choice (1 or 2): ").strip()
        
        if choice == '1':
            dna_input = input("\nEnter DNA sequence (use A, T, G, C): ").strip()
            
            result = convert_gene_to_protein(dna_input)
            
            if "error" in result:
                print(f"\nError: {result['error']}")
            else:
                if "warning" in result:
                    print(f"\nWarning: {result['warning']}")
                
                print(f"\nDNA Sequence:     {result['dna_sequence']}")
                print(f"RNA Sequence:     {result['rna_sequence']}")
                print(f"Amino Acids:      {' - '.join(result['amino_acids'])}")
                print(f"Protein Sequence: {result['protein_sequence']}")
                print(f"Length:           {result['length']} amino acids")
        
        elif choice == '2':
            print("\nGoodbye!")
            break
        else:
            print("\nInvalid choice. Please enter 1 or 2.")


if __name__ == "__main__":
    main()
