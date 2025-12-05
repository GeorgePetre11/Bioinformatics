import random
import time
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

def calculate_gc_content(sequence):
    """Calculate the GC content (C+G percentage) of a DNA sequence."""
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100 if len(sequence) > 0 else 0

def find_overlap(a, b, min_overlap=8):
    """Find the maximum overlap between two sequences."""
    max_ov = 0
    for i in range(min_overlap, min(len(a), len(b)) + 1):
        if a[-i:] == b[:i]:
            max_ov = i
    return max_ov

def greedy_assemble(reads, min_overlap=8, tries_per_round=500):
    """Assemble reads using a greedy algorithm and return the longest contig."""
    reads = reads.copy()
    while len(reads) > 1:
        best_a, best_b, best_ov = None, None, 0

        for _ in range(tries_per_round):
            i, j = random.sample(range(len(reads)), 2)
            ov = find_overlap(reads[i], reads[j], min_overlap)
            if ov > best_ov:
                best_a, best_b, best_ov = i, j, ov

        if best_ov < min_overlap:
            break

        merged = reads[best_a] + reads[best_b][best_ov:]
        new_reads = [
            reads[k]
            for k in range(len(reads))
            if k not in (best_a, best_b)
        ]
        new_reads.append(merged)
        reads = new_reads

    return max(reads, key=len)

def generate_random_sequence(length, gc_content):
    """
    Generate a random DNA sequence with a specific GC content percentage.
    
    Args:
        length: Length of the sequence to generate
        gc_content: Desired GC content as a percentage (0-100)
    
    Returns:
        A DNA sequence string
    """
    gc_count = int(length * gc_content / 100)
    at_count = length - gc_count
    
    # Distribute G and C equally
    g_count = gc_count // 2
    c_count = gc_count - g_count
    
    # Distribute A and T equally
    a_count = at_count // 2
    t_count = at_count - a_count
    
    # Create the sequence and shuffle it
    sequence = ['G'] * g_count + ['C'] * c_count + ['A'] * a_count + ['T'] * t_count
    random.shuffle(sequence)
    
    return ''.join(sequence)

def test_reconstruction_time(sequence, num_samples=500, min_overlap=8):
    """
    Test the reconstruction time for a given sequence.
    
    Args:
        sequence: The DNA sequence to fragment and reconstruct
        num_samples: Number of fragments to generate
        min_overlap: Minimum overlap for assembly
    
    Returns:
        Time taken for reconstruction in seconds
    """
    seq_len = len(sequence)
    
    # Generate random samples
    samples = [
        sequence[random.randint(0, seq_len - 150):][:random.randint(100, 150)]
        for _ in range(num_samples)
    ]
    
    # Time the reconstruction
    start_time = time.time()
    assembled = greedy_assemble(samples, min_overlap=min_overlap, tries_per_round=500)
    end_time = time.time()
    
    return end_time - start_time

def analyze_gc_correlation(sequence_length=10000, num_samples=500, num_tests_per_gc=5):
    """
    Analyze the correlation between GC content and reconstruction time.
    
    Args:
        sequence_length: Length of sequences to test
        num_samples: Number of fragments per test
        num_tests_per_gc: Number of repetitions for each GC content level
    
    Returns:
        Dictionary with results
    """
    gc_contents = np.arange(20, 81, 5)  # Test GC content from 20% to 80% in 5% steps
    results = []
    
    print("Starting GC content vs. reconstruction time analysis...")
    print(f"Sequence length: {sequence_length}")
    print(f"Number of samples per test: {num_samples}")
    print(f"Tests per GC level: {num_tests_per_gc}\n")
    
    for gc in gc_contents:
        times = []
        for test_num in range(num_tests_per_gc):
            print(f"Testing GC content: {gc}% (test {test_num + 1}/{num_tests_per_gc})", end='\r')
            
            # Generate a random sequence with the specified GC content
            sequence = generate_random_sequence(sequence_length, gc)
            
            # Verify actual GC content
            actual_gc = calculate_gc_content(sequence)
            
            # Test reconstruction time
            recon_time = test_reconstruction_time(sequence, num_samples)
            
            times.append(recon_time)
            results.append({
                'target_gc': gc,
                'actual_gc': actual_gc,
                'time': recon_time
            })
        
        avg_time = np.mean(times)
        std_time = np.std(times)
        print(f"GC content: {gc}% - Avg time: {avg_time:.3f}s ± {std_time:.3f}s")
    
    return results

def plot_results(results):
    """Create visualizations of the correlation analysis."""
    gc_values = [r['actual_gc'] for r in results]
    time_values = [r['time'] for r in results]
    
    # Calculate statistics
    correlation, p_value = stats.pearsonr(gc_values, time_values)
    slope, intercept, r_value, p_value_reg, std_err = stats.linregress(gc_values, time_values)
    
    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Scatter plot with regression line
    ax1.scatter(gc_values, time_values, alpha=0.6, s=50)
    
    # Add regression line
    x_line = np.array([min(gc_values), max(gc_values)])
    y_line = slope * x_line + intercept
    ax1.plot(x_line, y_line, 'r-', linewidth=2, label=f'y = {slope:.4f}x + {intercept:.4f}')
    
    ax1.set_xlabel('GC Content (%)', fontsize=12)
    ax1.set_ylabel('Reconstruction Time (seconds)', fontsize=12)
    ax1.set_title(f'GC Content vs Reconstruction Time\nPearson r = {correlation:.4f}, p-value = {p_value:.4e}', 
                  fontsize=13)
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Box plot grouped by GC content ranges
    gc_ranges = [(20, 35), (35, 50), (50, 65), (65, 80)]
    grouped_times = []
    labels = []
    
    for gc_min, gc_max in gc_ranges:
        group_times = [r['time'] for r in results if gc_min <= r['actual_gc'] < gc_max]
        if group_times:
            grouped_times.append(group_times)
            labels.append(f'{gc_min}-{gc_max}%')
    
    ax2.boxplot(grouped_times, labels=labels)
    ax2.set_xlabel('GC Content Range', fontsize=12)
    ax2.set_ylabel('Reconstruction Time (seconds)', fontsize=12)
    ax2.set_title('Reconstruction Time Distribution by GC Content Range', fontsize=13)
    ax2.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig('gc_correlation_analysis.png', dpi=300, bbox_inches='tight')
    print("\nPlot saved as 'gc_correlation_analysis.png'")
    
    return correlation, p_value, slope, intercept, r_value**2

def save_results(results, stats_info):
    """Save the results to a text file."""
    correlation, p_value, slope, intercept, r_squared = stats_info
    
    with open('gc_correlation_results.txt', 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("GC CONTENT vs RECONSTRUCTION TIME CORRELATION ANALYSIS\n")
        f.write("=" * 80 + "\n\n")
        
        f.write("STATISTICAL SUMMARY:\n")
        f.write("-" * 80 + "\n")
        f.write(f"Pearson Correlation Coefficient: {correlation:.6f}\n")
        f.write(f"P-value: {p_value:.6e}\n")
        f.write(f"R-squared: {r_squared:.6f}\n")
        f.write(f"Linear Regression: y = {slope:.6f}x + {intercept:.6f}\n")
        f.write(f"Slope: {slope:.6f} (seconds per % GC content)\n\n")
        
        # Interpret the correlation
        f.write("INTERPRETATION:\n")
        f.write("-" * 80 + "\n")
        
        if abs(correlation) < 0.3:
            strength = "weak"
        elif abs(correlation) < 0.7:
            strength = "moderate"
        else:
            strength = "strong"
        
        direction = "positive" if correlation > 0 else "negative"
        
        f.write(f"The correlation is {strength} and {direction}.\n")
        
        if p_value < 0.001:
            f.write("The correlation is highly statistically significant (p < 0.001).\n")
        elif p_value < 0.05:
            f.write("The correlation is statistically significant (p < 0.05).\n")
        else:
            f.write("The correlation is NOT statistically significant (p >= 0.05).\n")
        
        f.write("\nPOSSIBLE EXPLANATIONS:\n")
        f.write("-" * 80 + "\n")
        
        if abs(correlation) > 0.3:
            if correlation > 0:
                f.write("Higher GC content tends to increase reconstruction time. This could be because:\n")
                f.write("- G-C base pairs have stronger hydrogen bonds (3 vs 2 for A-T)\n")
                f.write("- GC-rich regions may have more repetitive patterns\n")
                f.write("- The overlap detection algorithm may behave differently with GC-rich sequences\n")
            else:
                f.write("Higher GC content tends to decrease reconstruction time. This could be because:\n")
                f.write("- GC-rich sequences may have fewer repetitive overlaps\n")
                f.write("- The greedy algorithm may find optimal overlaps faster\n")
                f.write("- AT-rich regions might have more ambiguous overlaps\n")
        else:
            f.write("GC content appears to have minimal impact on reconstruction time.\n")
            f.write("This suggests the greedy assembly algorithm's performance is largely\n")
            f.write("independent of sequence composition.\n")
        
        f.write("\n" + "=" * 80 + "\n")
        f.write("DETAILED RESULTS:\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"{'Target GC%':<12} {'Actual GC%':<12} {'Time (s)':<12}\n")
        f.write("-" * 80 + "\n")
        
        for r in results:
            f.write(f"{r['target_gc']:<12.1f} {r['actual_gc']:<12.2f} {r['time']:<12.6f}\n")
        
        f.write("\n" + "=" * 80 + "\n")
    
    print("Results saved to 'gc_correlation_results.txt'")

def main():
    """Main function to run the correlation analysis."""
    # Set random seed for reproducibility
    random.seed(42)
    np.random.seed(42)
    
    print("=" * 80)
    print("GC CONTENT vs RECONSTRUCTION TIME CORRELATION ANALYSIS")
    print("=" * 80 + "\n")
    
    # Run the analysis
    results = analyze_gc_correlation(
        sequence_length=10000,  # Length of test sequences
        num_samples=500,        # Number of fragments per test
        num_tests_per_gc=5      # Repetitions per GC level
    )
    
    print("\n" + "=" * 80)
    print("Analysis complete! Generating visualizations and summary...")
    print("=" * 80 + "\n")
    
    # Plot results and get statistics
    stats_info = plot_results(results)
    
    # Save results to file
    save_results(results, stats_info)
    
    correlation, p_value, slope, intercept, r_squared = stats_info
    
    print("\n" + "=" * 80)
    print("SUMMARY:")
    print("=" * 80)
    print(f"Pearson Correlation: {correlation:.4f}")
    print(f"P-value: {p_value:.4e}")
    print(f"R-squared: {r_squared:.4f}")
    print(f"Linear Model: Time = {slope:.6f} * GC% + {intercept:.6f}")
    print("=" * 80)

if __name__ == "__main__":
    main()
