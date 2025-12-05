#!/usr/bin/env python3

from restriction_enzyme_analysis import RestrictionEnzyme, RestrictionAnalyzer

def get_custom_sequence():
    print("\n" + "=" * 70)
    print("CUSTOM DNA SEQUENCE INPUT")
    print("=" * 70)
    print("Enter your DNA sequence (A, T, C, G only)")
    print("Press Enter twice when done:\n")
    
    lines = []
    while True:
        line = input()
        if line == "":
            break
        lines.append(line)
    
    sequence = "".join(lines).upper().replace(" ", "").replace("\n", "")
    
    valid_bases = set("ATCG")
    if not all(base in valid_bases for base in sequence):
        print("\n[!] Error: Invalid characters in sequence. Only A, T, C, G allowed.")
        return None
    
    print(f"\n[✓] Valid DNA sequence entered: {len(sequence)} nucleotides")
    return sequence

def add_custom_enzyme(analyzer):
    print("\n" + "=" * 70)
    print("ADD CUSTOM RESTRICTION ENZYME")
    print("=" * 70)
    
    name = input("Enter enzyme name: ").strip()
    recognition_seq = input("Enter recognition sequence (e.g., GAATTC): ").strip().upper()
    
    try:
        cleavage_pos = int(input("Enter cleavage position (0-based index): "))
    except ValueError:
        print("\n[!] Error: Cleavage position must be an integer.")
        return False
    
    valid_bases = set("ATCG")
    if not all(base in valid_bases for base in recognition_seq):
        print("\n[!] Error: Invalid recognition sequence. Only A, T, C, G allowed.")
        return False
    
    if cleavage_pos < 0 or cleavage_pos > len(recognition_seq):
        print("\n[!] Error: Cleavage position out of range.")
        return False
    
    analyzer.enzymes[name] = RestrictionEnzyme(name, recognition_seq, cleavage_pos)
    print(f"\n[✓] Custom enzyme '{name}' added successfully!")
    return True

def analyze_custom_dna():
    print("\n" + "=" * 70)
    print("CUSTOM DNA RESTRICTION ENZYME ANALYSIS TOOL")
    print("=" * 70)
    
    sequence = get_custom_sequence()
    if not sequence:
        return
    
    analyzer = RestrictionAnalyzer()
    
    print("\n" + "=" * 70)
    print("ENZYME SELECTION")
    print("=" * 70)
    print("\nAvailable standard enzymes:")
    for i, (name, enzyme) in enumerate(analyzer.enzymes.items(), 1):
        print(f"  {i}. {name} - {enzyme.recognition_seq}")
    
    print("\nOptions:")
    print("  [1] Use all standard enzymes")
    print("  [2] Select specific enzymes")
    print("  [3] Add custom enzyme")
    
    choice = input("\nYour choice (1/2/3): ").strip()
    
    if choice == "2":
        selected = input("Enter enzyme names (comma-separated): ").strip().split(",")
        selected = [s.strip() for s in selected]
        analyzer.enzymes = {k: v for k, v in analyzer.enzymes.items() if k in selected}
        
        if not analyzer.enzymes:
            print("\n[!] No valid enzymes selected. Using all standard enzymes.")
            analyzer = RestrictionAnalyzer()
    
    elif choice == "3":
        if not add_custom_enzyme(analyzer):
            return
    
    print("\n" + "=" * 70)
    print("ANALYZING...")
    print("=" * 70)
    
    results = analyzer.analyze(sequence)
    
    print("\n" + "=" * 70)
    print("ANALYSIS RESULTS")
    print("=" * 70)
    print(f"\nDNA Sequence Length: {len(sequence)} nucleotides\n")
    
    for enzyme_name, result in results.items():
        print(f"\n{'─' * 70}")
        print(f"Enzyme: {enzyme_name}")
        print(f"Recognition Sequence: {result['enzyme'].recognition_seq}")
        print(f"Number of Cleavages: {result['num_cleavages']}")
        print(f"Number of Fragments: {result['num_fragments']}")
        
        if result['cleavage_sites']:
            print(f"Cleavage Positions: {result['cleavage_sites']}")
        
        if result['fragments']:
            print(f"Fragment Lengths (sorted): {result['fragments']}")
    
    print("\n" + "=" * 70)
    
    visualize = input("\nGenerate visualizations? (y/n): ").strip().lower()
    
    if visualize == 'y':
        print("\n[*] Generating gel electrophoresis visualization...")
        analyzer.visualize_gel('custom_gel.png')
        
        print("\n[*] Generating fragment distribution chart...")
        analyzer.visualize_fragment_distribution('custom_distribution.png')
        
        print("\n[✓] Visualizations saved!")
    
    print("\n" + "=" * 70)
    print("Analysis Complete!")
    print("=" * 70 + "\n")

if __name__ == "__main__":
    analyze_custom_dna()
