import random

with open("gene.fna") as f:
    lines = f.readlines()

if lines[0].startswith(">"):
    lines = lines[1:]
sequence = "".join(line.strip() for line in lines)
seq_len = len(sequence)

NUM_SAMPLES = 2000
samples = [
    sequence[random.randint(0, seq_len - 150):][:random.randint(100, 150)]
    for _ in range(NUM_SAMPLES)
]

def find_overlap(a, b, min_overlap=8):
    max_ov = 0
    for i in range(min_overlap, min(len(a), len(b)) + 1):
        if a[-i:] == b[:i]:
            max_ov = i
    return max_ov

def greedy_assemble(reads, min_overlap=8, tries_per_round=500):
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

        print(f"Merged with overlap {best_ov}, remaining reads: {len(reads)}")

    return max(reads, key=len)

assembled = greedy_assemble(samples, min_overlap=8)

with open("reconstructed_gene.fna", "w") as f:
    f.write(">reconstructed_gene\n")
    for i in range(0, len(assembled), 70):
        f.write(assembled[i:i + 70] + "\n")

explanation = """\
Main Issue with This Algorithmic Approach:

Even though this script merges fragments that share overlaps ≥10 bases, it still fails to
reliably rebuild the original DNA sequence because:

1. The algorithm is greedy, so it can pick wrong overlaps when multiple fragments share
   similar subsequences.
2. Repetitive DNA regions make overlap decisions ambiguous.
3. Random sampling loses order and orientation information.
4. Comparing fragments is computationally expensive (O(n²)).
5. Sequencing errors or mutations would break overlaps completely.

Real genome assemblers fix these problems using graph-based methods (de Bruijn graphs or
overlap-layout-consensus) and probabilistic alignment scoring.
"""
with open("answer.txt", "w", encoding="utf-8") as f:
    f.write(explanation)

print("Done")
