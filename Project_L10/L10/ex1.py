import matplotlib.pyplot as plt
import numpy as np


def calculate_cg_content(sequence: str) -> float:
    """
    (C+G)% as defined in Gagniuc et al.
    CG% = 100 * (C + G) / (A + C + G + T)
    """
    sequence = sequence.upper()
    length = len(sequence)
    if length == 0:
        return 0.0

    c = sequence.count("C")
    g = sequence.count("G")
    a = sequence.count("A")
    t = sequence.count("T")

    total = a + c + g + t
    if total == 0:
        return 0.0

    return (c + g) / total * 100.0


def calculate_kappa_ic(sequence: str) -> float:
    """
    Kappa Index of Coincidence for a *single* sequence,
    implemented exactly as in Gagniuc's VB code (DNAKAPPA):

        Function IC(ByVal s1 As String) As Variant
            max = Len(s1) - 1
            For u = 1 To max
                s2 = Mid(s1, u + 1)
                For i = 1 To Len(s2)
                    If Mid(s1, i, 1) = Mid(s2, i, 1) Then
                        count = count + 1
                    End If
                Next i
                total = total + (count / Len(s2) * 100)
                count = 0
            Next u
            IC = Round((total / max), 2)
        End Function

    Here we implement the same logic in Python.
    """
    sequence = sequence.upper()
    n = len(sequence)
    if n <= 1:
        return 0.0

    max_shift = n - 1
    total = 0.0

    # u corresponds to the shift (1..max_shift in VB, 0-based here)
    for u in range(1, max_shift + 1):
        s2 = sequence[u:]  # suffix starting at position u
        count = 0

        # Compare prefix of original sequence to this suffix
        # (Mid(s1, i, 1) vs Mid(s2, i, 1))
        for i in range(len(s2)):
            if sequence[i] == s2[i]:
                count += 1

        total += (count / len(s2)) * 100.0

    # VB uses Round(total / max, 2). We return the float
    ic_value = total / max_shift
    return ic_value


def sliding_window_analysis(sequence: str, window_size: int = 30):
    """
    Perform sliding-window analysis (step = 1 nt).
    Returns:
        positions: start index of each window
        cg_values: (C+G)% for each window
        kappa_values: Kappa IC for each window
    """
    sequence = sequence.upper()
    n = len(sequence)
    if n < window_size:
        return [], [], []

    positions = []
    cg_values = []
    kappa_values = []

    for start in range(n - window_size + 1):
        window = sequence[start:start + window_size]
        cg = calculate_cg_content(window)
        kappa = calculate_kappa_ic(window)

        positions.append(start)
        cg_values.append(cg)
        kappa_values.append(kappa)

    return positions, cg_values, kappa_values


def calculate_center_of_weight(x_values, y_values):
    """
    Center of weight (mean point) of the pattern.
    """
    if not x_values or not y_values:
        return 0.0, 0.0

    center_x = float(np.mean(x_values))
    center_y = float(np.mean(y_values))
    return center_x, center_y


def plot_pattern(cg_values, kappa_values, title="DNA Pattern: (C+G)% vs Kappa IC"):
    """
    Scatter plot of (C+G)% vs Kappa IC for all sliding windows.
    """
    plt.figure(figsize=(10, 8))
    plt.scatter(cg_values, kappa_values, alpha=0.6, s=50,
                edgecolors='black')

    plt.xlabel('(C+G) Content (%)', fontsize=12)
    plt.ylabel('Kappa Index of Coincidence (%)', fontsize=12)
    plt.title(title, fontsize=14, weight='bold')
    plt.grid(True, alpha=0.3)

    # Center of weight
    center_x, center_y = calculate_center_of_weight(cg_values, kappa_values)
    plt.scatter([center_x], [center_y],
                s=200, marker='X', color='red',
                edgecolors='black', linewidths=2,
                label=f'Center ({center_x:.2f}, {center_y:.2f})')
    plt.legend()

    plt.tight_layout()
    plt.show()


def plot_centers(centers, labels=None):
    """
    Plot the centers of multiple patterns on one chart.
    centers: list of (cg_center, kappa_center)
    labels: list of labels for each pattern
    """
    plt.figure(figsize=(10, 8))

    for i, (cg, kappa) in enumerate(centers):
        label = labels[i] if labels and i < len(labels) else f'Pattern {i+1}'
        plt.scatter([cg], [kappa],
                    s=200, marker='X',
                    edgecolors='black', linewidths=2,
                    label=label, alpha=0.8)

    plt.xlabel('(C+G) Content (%)', fontsize=12)
    plt.ylabel('Kappa Index of Coincidence (%)', fontsize=12)
    plt.title('Pattern Centers Comparison', fontsize=14, weight='bold')
    plt.grid(True, alpha=0.3)
    plt.legend()

    plt.tight_layout()
    plt.show()


def main():
    # Test sequence from the assignment
    S = "CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGATCTATCTAGCGATCTATCTACTACG"
    window_size = 30

    print("=" * 70)
    print("DNA PATTERN ANALYSIS - (C+G)% and Kappa Index of Coincidence")
    print("=" * 70)
    print()
    print(f"Sequence: {S}")
    print(f"Sequence length: {len(S)} bp")
    print(f"Window size: {window_size} bp")
    print()

    # Global (C+G)% and KappaIC for the entire sequence
    global_cg = calculate_cg_content(S)
    global_kappa = calculate_kappa_ic(S)
    print(f"Global (C+G)% for full sequence: {global_cg:.2f}%  (expected ≈ 29.27)")
    print(f"Global Kappa IC for full sequence: {global_kappa:.2f}%")
    print()

    # Sliding-window analysis
    positions, cg_values, kappa_values = sliding_window_analysis(S, window_size)
    print(f"Number of windows: {len(cg_values)}")
    print()

    avg_cg = float(np.mean(cg_values)) if cg_values else 0.0
    avg_kappa = float(np.mean(kappa_values)) if kappa_values else 0.0

    print(f"Average (C+G)% across windows: {avg_cg:.2f}%")
    print(f"Average Kappa IC across windows: {avg_kappa:.2f}%")
    print()

    # Test on the first window
    first_window = S[0:window_size]
    test_cg = calculate_cg_content(first_window)
    test_kappa = calculate_kappa_ic(first_window)

    print("First window test:")
    print(f"  Window: {first_window}")
    print(f"  (C+G)%: {test_cg:.2f}%")
    print(f"  Kappa IC: {test_kappa:.2f}%")
    print()

    # Center of weight of the whole pattern
    center_cg, center_kappa = calculate_center_of_weight(cg_values, kappa_values)
    print(f"Center of weight of pattern: ({center_cg:.2f}, {center_kappa:.2f})")
    print()

    # Plot pattern and centers
    print("Generating pattern plot...")
    plot_pattern(cg_values, kappa_values,
                 "DNA Pattern: (C+G)% vs Kappa IC")

    centers = [(center_cg, center_kappa)]
    labels = ["Test Sequence"]

    print("Generating centers plot...")
    plot_centers(centers, labels)

    print()
    print("=" * 70)
    print("Analysis complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
