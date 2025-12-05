import matplotlib.pyplot as plt
from dna_pattern_detector import DNAPatternDetector
import os
from collections import Counter

def analyze_strain(filename):
    with open(filename, 'r') as f:
        content = f.read()
    
    first_line = content.split('\n')[0]
    strain_name = first_line.replace('>', '').strip()
    
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
    
    return strain_name, pattern_counts, len(detector.sequence), len(results)

def extract_short_name(full_name):
    if "California" in full_name and "H1N1" in full_name:
        return "H1N1 (2009)"
    elif "Beijing" in full_name and "H3N2" in full_name:
        return "H3N2 (1992)"
    elif "Brisbane" in full_name:
        return "H3 (2006)"
    elif "Hangzhou" in full_name and "H7N9" in full_name:
        return "H7N9 (2013)"
    elif "Hong Kong" in full_name and "Influenza B" in full_name:
        return "Influenza B (HK)"
    elif "Influenza B" in full_name:
        return "Influenza B"
    elif "WSN" in full_name:
        return "H1N1 (1933)"
    elif "Oesophagostomum" in full_name:
        return "Sample 9"
    elif "parvovirus" in full_name:
        return "Sample 10"
    else:
        parts = full_name.split()
        return " ".join(parts[:3]) if len(parts) >= 3 else full_name[:20]

def create_overview_chart():
    print("="*80)
    print("CREATING OVERVIEW CHART FOR ALL INFLUENZA STRAINS")
    print("="*80)
    print()
    
    fasta_files = [f'sequence_influenza_{i}.fasta' for i in range(1, 11)]
    
    strain_names = []
    total_patterns_list = []
    unique_patterns_list = []
    sequence_lengths = []
    
    for i, filename in enumerate(fasta_files, 1):
        if not os.path.exists(filename):
            print(f"⚠ File not found: {filename}")
            continue
        
        print(f"[{i}/10] Processing {filename}...")
        
        try:
            full_name, pattern_counts, seq_length, total_results = analyze_strain(filename)
            short_name = extract_short_name(full_name)
            
            total_patterns = sum(pattern_counts.values())
            unique_patterns = len(pattern_counts)
            
            strain_names.append(f"Strain {i}\n{short_name}")
            total_patterns_list.append(total_patterns)
            unique_patterns_list.append(unique_patterns)
            sequence_lengths.append(seq_length)
            
            print(f"  ✓ {short_name}: {total_patterns} total, {unique_patterns} unique")
            
        except Exception as e:
            print(f"  ✗ Error: {str(e)}")
    
    fig, axes = plt.subplots(2, 2, figsize=(18, 12))
    fig.suptitle('Influenza Virus Strains - Pattern Repetition Overview', 
                 fontsize=16, fontweight='bold', y=0.995)
    
    x_pos = range(len(strain_names))
    
    ax1 = axes[0, 0]
    bars1 = ax1.bar(x_pos, total_patterns_list, color='steelblue', edgecolor='navy', alpha=0.7)
    for bar, count in zip(bars1, total_patterns_list):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1, 
                str(count), ha='center', va='bottom', fontsize=9, fontweight='bold')
    ax1.set_xlabel('Strain', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Total Pattern Occurrences', fontsize=11, fontweight='bold')
    ax1.set_title('Total Pattern Repetitions Found', fontsize=12, fontweight='bold', pad=10)
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(strain_names, rotation=45, ha='right', fontsize=8)
    ax1.grid(axis='y', alpha=0.3, linestyle='--')
    
    ax2 = axes[0, 1]
    bars2 = ax2.bar(x_pos, unique_patterns_list, color='coral', edgecolor='darkred', alpha=0.7)
    for bar, count in zip(bars2, unique_patterns_list):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5, 
                str(count), ha='center', va='bottom', fontsize=9, fontweight='bold')
    ax2.set_xlabel('Strain', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Unique Patterns', fontsize=11, fontweight='bold')
    ax2.set_title('Number of Unique Patterns', fontsize=12, fontweight='bold', pad=10)
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(strain_names, rotation=45, ha='right', fontsize=8)
    ax2.grid(axis='y', alpha=0.3, linestyle='--')
    
    ax3 = axes[1, 0]
    bars3 = ax3.bar(x_pos, sequence_lengths, color='mediumseagreen', edgecolor='darkgreen', alpha=0.7)
    for bar, length in zip(bars3, sequence_lengths):
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 50, 
                f'{length:,}', ha='center', va='bottom', fontsize=8, fontweight='bold')
    ax3.set_xlabel('Strain', fontsize=11, fontweight='bold')
    ax3.set_ylabel('Sequence Length (bp)', fontsize=11, fontweight='bold')
    ax3.set_title('Genome Sequence Length', fontsize=12, fontweight='bold', pad=10)
    ax3.set_xticks(x_pos)
    ax3.set_xticklabels(strain_names, rotation=45, ha='right', fontsize=8)
    ax3.grid(axis='y', alpha=0.3, linestyle='--')
    
    ax4 = axes[1, 1]
    pattern_density = [total/length*1000 if length > 0 else 0 
                       for total, length in zip(total_patterns_list, sequence_lengths)]
    bars4 = ax4.bar(x_pos, pattern_density, color='mediumpurple', edgecolor='indigo', alpha=0.7)
    for bar, density in zip(bars4, pattern_density):
        ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5, 
                f'{density:.1f}', ha='center', va='bottom', fontsize=9, fontweight='bold')
    ax4.set_xlabel('Strain', fontsize=11, fontweight='bold')
    ax4.set_ylabel('Pattern Density (per 1000 bp)', fontsize=11, fontweight='bold')
    ax4.set_title('Pattern Density Analysis', fontsize=12, fontweight='bold', pad=10)
    ax4.set_xticks(x_pos)
    ax4.set_xticklabels(strain_names, rotation=45, ha='right', fontsize=8)
    ax4.grid(axis='y', alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    
    output_filename = 'influenza_strains_overview.png'
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print("\n" + "="*80)
    print(f"✓ Overview chart saved: {output_filename}")
    print("="*80)
    plt.close()

if __name__ == "__main__":
    create_overview_chart()
