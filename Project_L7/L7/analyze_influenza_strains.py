import matplotlib.pyplot as plt
from dna_pattern_detector import DNAPatternDetector
import os
from collections import Counter

def analyze_strain(filename):
    with open(filename, 'r') as f:
        content = f.read()
    
    first_line = content.split('\n')[0]
    strain_name = first_line.replace('>', '').strip()
    if len(strain_name) > 80:
        strain_name = strain_name[:80] + "..."
    
    detector = DNAPatternDetector(content)
    
    results = detector.find_consecutive_repetitions(
        min_pattern_size=3,
        max_pattern_size=6,
        min_repetitions=2
    )
    
    pattern_counts = Counter()
    for result in results:
        pattern = result['pattern']
        pattern_counts[pattern] += 1
    
    return strain_name, pattern_counts, len(detector.sequence)

def plot_top_patterns(strain_name, pattern_counts, sequence_length, file_number):
    top_20 = pattern_counts.most_common(20)
    
    if not top_20:
        print(f"No patterns found for {strain_name}")
        return
    
    patterns = [item[0] for item in top_20]
    counts = [item[1] for item in top_20]
    
    plt.figure(figsize=(14, 8))
    bars = plt.bar(range(len(patterns)), counts, color='steelblue', edgecolor='navy', alpha=0.7)
    
    for i, (bar, count) in enumerate(zip(bars, counts)):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5, 
                str(count), ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    plt.xlabel('DNA Pattern', fontsize=12, fontweight='bold')
    plt.ylabel('Frequency (Number of Occurrences)', fontsize=12, fontweight='bold')
    plt.title(f'Top 20 Most Frequent DNA Pattern Repetitions\n{strain_name}\n(Sequence Length: {sequence_length:,} bp)', 
              fontsize=13, fontweight='bold', pad=20)
    
    plt.xticks(range(len(patterns)), patterns, rotation=45, ha='right', fontsize=10)
    plt.yticks(fontsize=10)
    
    plt.grid(axis='y', alpha=0.3, linestyle='--')
    plt.tight_layout()
    
    output_filename = f'influenza_strain_{file_number}_chart.png'
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"✓ Chart saved: {output_filename}")
    plt.close()

def main():
    print("="*80)
    print("INFLUENZA VIRUS STRAINS - PATTERN REPETITION ANALYSIS")
    print("="*80)
    print()
    
    fasta_files = [f'sequence_influenza_{i}.fasta' for i in range(1, 11)]
    
    for i, filename in enumerate(fasta_files, 1):
        if not os.path.exists(filename):
            print(f"⚠ File not found: {filename}")
            continue
        
        print(f"\n[{i}/10] Analyzing {filename}...")
        
        try:
            strain_name, pattern_counts, seq_length = analyze_strain(filename)
            
            total_patterns = sum(pattern_counts.values())
            unique_patterns = len(pattern_counts)
            
            print(f"  Strain: {strain_name[:60]}...")
            print(f"  Sequence length: {seq_length:,} bp")
            print(f"  Total pattern occurrences: {total_patterns}")
            print(f"  Unique patterns found: {unique_patterns}")
            
            plot_top_patterns(strain_name, pattern_counts, seq_length, i)
            
        except Exception as e:
            print(f"  ✗ Error analyzing {filename}: {str(e)}")
    
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE!")
    print(f"Generated 10 charts showing top 20 pattern repetitions for each strain.")
    print("="*80)

if __name__ == "__main__":
    main()
