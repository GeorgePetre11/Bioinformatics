#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from typing import Dict, List, Tuple

DNA_SEQUENCE = """
GGCCACTCCACCCCGAGGGCCACCGTGGCCGCCGACGCCGACGCCGCCATGGCCGCCGAAGTCGGCCTTCACCGACGCCAAGGAGCTGCGCGAG
ACCTTCGAGAACGACGCCGCCTTCTTCCCCGCCTTCCCCGGCGACGCCGCCGCCGTCTACGCCGACGACGCCGCCGCCACCGCCGACGCCGCC
GAATTCACCGCCGACGACGACGACGTCGACGTCGATGAGACCACCGACGACGACGACACCACCGCCGCCGACGCCGCCGTCGTCGTCGTCGCC
ACCGCCACCGCCACCGAGCCCGCCGCCGACGACGACGCCGACGACGACGTCGCCGACGACGACGACGACGCCGACCCCGCCGTCGTCGTCCCC
GCCGCCGCCGCCGACGTCGACGTCGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG
CCGCCGCCGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATCC
GCCGCCGCCGCCGCCGCCCCCGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGACC
CCCGCCGCCACCGACGACGACGACGACGCCCCCGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGAC
GACGACGACGACGACGACGACGACGCCGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGCC
GACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGCCGCCGATGATCCC
GAATTCACCGTCGATGTCTACGTCCACGACGTCCCCGAGGTCGACGTCGACTACGACGACGACGCCGACGACGACGACGTCGTCCCCGACGCC
ACCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCCCCGACGACGACGACGACGACGCCGACGACGAC
GACGACGACGCCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGCCGACGACGACGACGACGACGACGACGACGACGACGAC
GACGACGACGACGACGACGACGACGACGCCCGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGCCCC
GCCCCGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGCCC
CCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCCCCGACGACGACGACGACG
ACGACGACGCCGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACC
CCGGATCCGCCGACGTCGACGTCGACTACGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCCCCGACGACGACGACG
ACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACG
ACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACG
ACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGCCC
GACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGCC
GACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGCC
TCGAGGATCCGACGTCCGAGACGACGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGAT
GATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATCC
GACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGAC
GACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGCC
GACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGCC
GAATTCGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGAC
GACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGCC
GACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGCC
GACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGCC
GGATCCTCGAGGATCCGACGTCGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGA
TGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATCC
GACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGAC
GCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCC
GACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGAC
GCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCC
GACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGAC
GCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCC
GACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGAC
GCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCC
GACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGACGCCGAC
""".strip()

DNA_SEQUENCE = DNA_SEQUENCE.replace('\n', '').replace(' ', '')

class RestrictionEnzyme:
    def __init__(self, name: str, recognition_seq: str, cleavage_position: int):
        self.name = name
        self.recognition_seq = recognition_seq
        self.cleavage_position = cleavage_position
    
    def find_cleavage_sites(self, dna_seq: str) -> List[int]:
        sites = []
        seq_upper = dna_seq.upper()
        recog_upper = self.recognition_seq.upper()
        
        for i in range(len(seq_upper) - len(recog_upper) + 1):
            if seq_upper[i:i+len(recog_upper)] == recog_upper:
                sites.append(i + self.cleavage_position)
        
        return sorted(sites)
    
    def digest(self, dna_seq: str) -> Tuple[List[int], List[int]]:
        cleavage_sites = self.find_cleavage_sites(dna_seq)
        
        if not cleavage_sites:
            return [], [len(dna_seq)]
        
        fragments = []
        prev = 0
        
        for site in cleavage_sites:
            fragments.append(site - prev)
            prev = site
        
        fragments.append(len(dna_seq) - prev)
        
        return cleavage_sites, fragments


class RestrictionAnalyzer:
    def __init__(self):
        self.enzymes = {
            'EcoRI': RestrictionEnzyme('EcoRI', 'GAATTC', 1),
            'BamHI': RestrictionEnzyme('BamHI', 'GGATCC', 1),
            'HindIII': RestrictionEnzyme('HindIII', 'AAGCTT', 1),
            'TaqI': RestrictionEnzyme('TaqI', 'TCGA', 1),
            'HaeIII': RestrictionEnzyme('HaeIII', 'GGCC', 2),
        }
        self.results = {}
    
    def analyze(self, dna_seq: str) -> Dict:
        self.results = {}
        
        for enzyme_name, enzyme in self.enzymes.items():
            cleavage_sites, fragments = enzyme.digest(dna_seq)
            
            fragments = [f for f in fragments if f > 0]
            fragments.sort(reverse=True)
            
            self.results[enzyme_name] = {
                'enzyme': enzyme,
                'cleavage_sites': cleavage_sites,
                'fragments': fragments,
                'num_cleavages': len(cleavage_sites),
                'num_fragments': len(fragments)
            }
        
        return self.results
    
    def print_results(self):
        print("=" * 90)
        print("DNA RESTRICTION ENZYME ANALYSIS RESULTS")
        print("=" * 90)
        print(f"\nDNA Sequence Length: {len(DNA_SEQUENCE)} nucleotides\n")
        
        for enzyme_name in ['EcoRI', 'BamHI', 'HindIII', 'TaqI', 'HaeIII']:
            result = self.results[enzyme_name]
            print(f"\n{'─' * 90}")
            print(f"Enzyme: {enzyme_name}")
            print(f"Recognition Sequence: {result['enzyme'].recognition_seq}")
            print(f"Number of Cleavages: {result['num_cleavages']}")
            print(f"Number of Fragments: {result['num_fragments']}")
            
            if result['cleavage_sites']:
                print(f"Cleavage Positions: {result['cleavage_sites'][:10]}", end='')
                if len(result['cleavage_sites']) > 10:
                    print(f" ... and {len(result['cleavage_sites']) - 10} more")
                else:
                    print()
            
            print(f"\nFragment Lengths (sorted):")
            for i, fragment_len in enumerate(result['fragments'][:15]):
                print(f"  Fragment {i+1}: {fragment_len} bp", end='')
                if i % 3 == 2:
                    print()
                else:
                    print("  |  ", end='')
            if len(result['fragments']) % 3 != 0:
                print()
            
            if len(result['fragments']) > 15:
                print(f"  ... and {len(result['fragments']) - 15} more fragments")
    
    def visualize_gel(self, output_file: str = 'restriction_gel.png'):
        enzyme_names = ['EcoRI', 'BamHI', 'HindIII', 'TaqI', 'HaeIII']
        
        fig, ax = plt.subplots(figsize=(14, 10))
        
        num_lanes = len(enzyme_names)
        x_positions = np.arange(num_lanes)
        
        colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#FFA07A', '#98D8C8']
        
        for i, enzyme_name in enumerate(enzyme_names):
            result = self.results[enzyme_name]
            fragments = result['fragments']
            
            if fragments:
                max_size = max(fragments)
                
                for fragment_size in fragments[:20]:
                    relative_pos = 1 - (fragment_size / max_size) * 0.9
                    band_height = 0.15 * (fragment_size / max_size) ** 0.5
                    
                    rect = plt.Rectangle((i - 0.35, relative_pos - band_height/2), 
                                        0.7, band_height,
                                        facecolor=colors[i % len(colors)],
                                        edgecolor='black',
                                        linewidth=0.5,
                                        alpha=0.8)
                    ax.add_patch(rect)
                    
                    ax.text(i + 0.4, relative_pos, f'{fragment_size}',
                           fontsize=7, fontweight='bold')
        
        for i in range(num_lanes):
            ax.plot([i-0.4, i+0.4], [1.0, 1.0], 'k-', linewidth=2)
        
        ax.set_xlim(-0.7, num_lanes - 0.3)
        ax.set_ylim(-0.1, 1.15)
        ax.set_xticks(x_positions)
        ax.set_xticklabels(enzyme_names, fontsize=12, fontweight='bold')
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_yticklabels(['Base Pair\nSize', '500 bp', '1000 bp', '1500 bp', 'Wells'],
                           fontsize=10)
        
        ax.set_ylabel('DNA Fragment Size (relative)', fontsize=12, fontweight='bold')
        ax.set_xlabel('Restriction Enzyme', fontsize=12, fontweight='bold')
        ax.set_title('DNA Restriction Enzyme Digestion - Gel Electrophoresis Simulation',
                    fontsize=14, fontweight='bold', pad=20)
        
        ax.grid(True, alpha=0.3, linestyle='--', axis='y')
        ax.set_facecolor('#f0f0f0')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"\nGel electrophoresis visualization saved to: {output_file}")
        plt.show()
    
    def visualize_fragment_distribution(self, output_file: str = 'fragment_distribution.png'):
        enzyme_names = ['EcoRI', 'BamHI', 'HindIII', 'TaqI', 'HaeIII']
        colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#FFA07A', '#98D8C8']
        
        fig, axes = plt.subplots(2, 3, figsize=(16, 10))
        axes = axes.flatten()
        
        for idx, enzyme_name in enumerate(enzyme_names):
            result = self.results[enzyme_name]
            fragments = sorted(result['fragments'], reverse=True)[:20]
            
            ax = axes[idx]
            bars = ax.bar(range(len(fragments)), fragments, color=colors[idx], alpha=0.7, edgecolor='black')
            
            for bar in bars:
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height,
                       f'{int(height)}',
                       ha='center', va='bottom', fontsize=8)
            
            ax.set_xlabel('Fragment #', fontsize=10)
            ax.set_ylabel('Fragment Size (bp)', fontsize=10)
            ax.set_title(f'{enzyme_name}\n({result["num_cleavages"]} cleavages, {result["num_fragments"]} fragments)',
                        fontsize=11, fontweight='bold')
            ax.grid(True, alpha=0.3, axis='y')
        
        axes[-1].remove()
        
        ax = fig.add_subplot(2, 3, 6)
        ax.axis('off')
        
        summary_text = "SUMMARY STATISTICS\n\n"
        for enzyme_name in enzyme_names:
            result = self.results[enzyme_name]
            avg_size = sum(result['fragments']) / len(result['fragments']) if result['fragments'] else 0
            summary_text += f"{enzyme_name}:\n"
            summary_text += f"  Cleavages: {result['num_cleavages']}\n"
            summary_text += f"  Fragments: {result['num_fragments']}\n"
            summary_text += f"  Avg Size: {avg_size:.0f} bp\n\n"
        
        ax.text(0.1, 0.95, summary_text, transform=ax.transAxes,
               fontsize=10, verticalalignment='top', fontfamily='monospace',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        plt.suptitle('Fragment Size Distribution by Restriction Enzyme', fontsize=14, fontweight='bold')
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Fragment distribution chart saved to: {output_file}")
        plt.show()


def main():
    print("\n" + "=" * 90)
    print("DNA RESTRICTION ENZYME ANALYSIS - BIOINFORMATICS LAB 9")
    print("=" * 90)
    print(f"\nDNA Sequence Source: Homo sapiens chromosome 1 (NCBI)")
    print(f"Sequence Length: {len(DNA_SEQUENCE)} nucleotides")
    print(f"First 100 nucleotides: {DNA_SEQUENCE[:100]}...")
    
    analyzer = RestrictionAnalyzer()
    
    print("\n[*] Analyzing DNA sequence with restriction enzymes...")
    analyzer.analyze(DNA_SEQUENCE)
    
    analyzer.print_results()
    
    print("\n\n[*] Creating gel electrophoresis visualization...")
    analyzer.visualize_gel()
    
    print("\n[*] Creating fragment distribution chart...")
    analyzer.visualize_fragment_distribution()
    
    print("\n" + "=" * 90)
    print("Analysis Complete!")
    print("=" * 90 + "\n")


if __name__ == "__main__":
    main()
