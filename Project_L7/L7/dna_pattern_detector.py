class DNAPatternDetector:
    def __init__(self, sequence):
        if sequence.startswith('>'):
            lines = sequence.split('\n')
            sequence = ''.join(lines[1:])
        
        self.sequence = sequence.upper().replace('\n', '').replace('\r', '').replace(' ', '').replace('\t', '')
        self.sequence = ''.join(c for c in self.sequence if c in 'ATCG')
    
    def find_consecutive_repetitions(self, min_pattern_size=3, max_pattern_size=6, min_repetitions=2):
        results = []
        sequence_len = len(self.sequence)
        
        used_positions = set()
        
        for pattern_size in range(max_pattern_size, min_pattern_size - 1, -1):
            position = 0
            
            while position <= sequence_len - pattern_size:
                if position in used_positions:
                    position += 1
                    continue
                
                pattern = self.sequence[position:position + pattern_size]
                
                repetition_count = 1
                current_pos = position + pattern_size
                
                while current_pos + pattern_size <= sequence_len:
                    next_segment = self.sequence[current_pos:current_pos + pattern_size]
                    if next_segment == pattern:
                        repetition_count += 1
                        current_pos += pattern_size
                    else:
                        break
                
                if repetition_count >= min_repetitions:
                    end_position = position + (pattern_size * repetition_count)
                    results.append({
                        'pattern': pattern,
                        'start': position,
                        'end': end_position,
                        'count': repetition_count,
                        'pattern_size': pattern_size
                    })
                    
                    for p in range(position, end_position):
                        used_positions.add(p)
                    
                    position = end_position
                else:
                    position += 1
        
        results.sort(key=lambda x: x['start'])
        return results


def main():
    print("="*70)
    print("DNA PATTERN REPETITION DETECTOR")
    print("="*70)
    print()
    
    try:
        with open('sequence.fasta', 'r') as f:
            content = f.read()
            detector = DNAPatternDetector(content)
            
        print(f"Loaded DNA sequence: {len(detector.sequence)} nucleotides")
        print()
        print("DNA Sequence Preview (first 200 bp):")
        print("-" * 70)
        preview = detector.sequence[:200]
        for i in range(0, len(preview), 60):
            print(preview[i:i+60])
        print("-" * 70)
        print()
        
        print("Analysis Parameters:")
        print(f"  - Pattern size range: 3-6 base pairs")
        print(f"  - Minimum consecutive repetitions: 3")
        print()
        
        print("Analyzing sequence for consecutive pattern repetitions...")
        print()
        
        results = detector.find_consecutive_repetitions(
            min_pattern_size=3,
            max_pattern_size=6,
            min_repetitions=3
        )
        
        if not results:
            print("No consecutive repetitions found with the specified criteria.")
            return
        
        print(f"RESULTS: Found {len(results)} repetition pattern(s)")
        print("="*70)
        print()
        
        for i, result in enumerate(results, 1):
            pattern = result['pattern']
            start = result['start']
            end = result['end']
            count = result['count']
            pattern_size = result['pattern_size']
            
            print(f"Repetition #{i}")
            print("-" * 70)
            print(f"  Pattern:     {pattern} (size: {pattern_size} bp)")
            print(f"  Position:    {start} - {end}")
            print(f"  Repetitions: {count}x")
            print(f"  Visual:      {' '.join([pattern] * count)}")
            print(f"  Full seq:    {detector.sequence[start:end]}")
            print()
        
        print("="*70)
        print("Analysis complete!")
        print()
        
        print("="*70)
        print("EXAMPLES BY PATTERN SIZE:")
        print("="*70)
        print()
        
        by_size = {}
        for result in results:
            size = result['pattern_size']
            if size not in by_size:
                by_size[size] = []
            by_size[size].append(result)
        
        for size in sorted(by_size.keys(), reverse=True):
            patterns = by_size[size]
            print(f"--- {size} BASE PAIR PATTERNS ({size}bp) ---")
            print(f"Total found: {len(patterns)}")
            print()
            
            for i, result in enumerate(patterns[:3], 1):
                pattern = result['pattern']
                start = result['start']
                end = result['end']
                count = result['count']
                
                print(f"Example {i}:")
                print(f"  Pattern:     {pattern}")
                print(f"  Position:    {start} - {end}")
                print(f"  Repetitions: {count}x")
                print(f"  Visual:      {' '.join([pattern] * count)}")
                print(f"  Full seq:    {detector.sequence[start:end]}")
                print()
            
            if len(patterns) > 3:
                print(f"  ... and {len(patterns) - 3} more {size}bp patterns")
                print()
        
        print("="*70)
        
    except FileNotFoundError:
        print("ERROR: sequence.fasta file not found!")
    except Exception as e:
        print(f"ERROR: {str(e)}")


if __name__ == "__main__":
    main()
