"""
Advanced Scoring Methods for Sequence Alignment Similarity Assessment

This module implements three different scoring equations to quantify
the level of similarity between aligned sequences:

1. Normalized Alignment Score (NAS) - Score normalized by optimal score
2. Weighted Similarity Index (WSI) - Considers match quality and gap patterns
3. Conservation Score (CS) - Biological conservation metric
"""

import numpy as np


class SequenceSimilarityScorer:
    """
    Implements multiple scoring methods to assess sequence similarity.
    Each method provides a different perspective on alignment quality.
    """
    
    def __init__(self, match_score=2, mismatch_score=-1, gap_score=-2):
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_score = gap_score
    
    def calculate_all_scores(self, aligned_seq1, aligned_seq2, raw_score=None):
        """
        Calculate all three similarity scores for an alignment.
        
        Args:
            aligned_seq1: First aligned sequence (with gaps as '-')
            aligned_seq2: Second aligned sequence (with gaps as '-')
            raw_score: Raw alignment score (optional, will be calculated if not provided)
        
        Returns:
            Dictionary with all three scores and their interpretations
        """
        if len(aligned_seq1) != len(aligned_seq2):
            raise ValueError("Aligned sequences must have the same length")
        
        # Calculate raw score if not provided
        if raw_score is None:
            raw_score = self._calculate_raw_score(aligned_seq1, aligned_seq2)
        
        # Calculate each scoring method
        nas = self.normalized_alignment_score(aligned_seq1, aligned_seq2, raw_score)
        wsi = self.weighted_similarity_index(aligned_seq1, aligned_seq2)
        cs = self.conservation_score(aligned_seq1, aligned_seq2)
        
        # Get basic statistics
        stats = self._get_alignment_statistics(aligned_seq1, aligned_seq2)
        
        return {
            'NAS': nas,
            'WSI': wsi,
            'CS': cs,
            'raw_score': raw_score,
            'statistics': stats,
            'interpretation': self._interpret_scores(nas, wsi, cs)
        }
    
    def normalized_alignment_score(self, aligned_seq1, aligned_seq2, raw_score=None):
        """
        Scoring Method 1: Normalized Alignment Score (NAS)
        
        Formula: NAS = (Actual_Score - Min_Score) / (Max_Score - Min_Score)
        
        Where:
        - Actual_Score: The raw alignment score obtained
        - Max_Score: Maximum possible score (all matches, no gaps)
        - Min_Score: Minimum possible score (all mismatches/gaps)
        
        Range: [0, 1] where 1 indicates perfect similarity
        
        This score normalizes the raw alignment score to account for sequence
        length and scoring parameters, making it comparable across different
        alignments.
        """
        if raw_score is None:
            raw_score = self._calculate_raw_score(aligned_seq1, aligned_seq2)
        
        alignment_length = len(aligned_seq1)
        
        # Calculate maximum possible score (all matches)
        max_score = alignment_length * self.match_score
        
        # Calculate minimum possible score (all mismatches or gaps)
        min_score = alignment_length * min(self.mismatch_score, self.gap_score)
        
        # Normalize the score
        if max_score == min_score:
            return 1.0 if raw_score >= 0 else 0.0
        
        nas = (raw_score - min_score) / (max_score - min_score)
        
        # Clamp to [0, 1] range
        return max(0.0, min(1.0, nas))
    
    def weighted_similarity_index(self, aligned_seq1, aligned_seq2):
        """
        Scoring Method 2: Weighted Similarity Index (WSI)
        
        Formula: WSI = (w_m * M + w_mm * MM + w_g * G) / L
        
        Where:
        - M: Number of exact matches
        - MM: Number of mismatches
        - G: Number of gaps
        - L: Alignment length
        - w_m, w_mm, w_g: Weights derived from scoring parameters
        
        Normalized weights:
        - w_m = 1.0 (matches contribute positively)
        - w_mm = mismatch_score / match_score (mismatches reduce score)
        - w_g = gap_score / match_score (gaps reduce score)
        
        Range: [-1, 1] where 1 indicates perfect similarity
        
        This score provides a weighted assessment that considers the relative
        importance of matches, mismatches, and gaps based on the scoring scheme.
        """
        matches = 0
        mismatches = 0
        gaps = 0
        
        for a, b in zip(aligned_seq1, aligned_seq2):
            if a == '-' or b == '-':
                gaps += 1
            elif a == b:
                matches += 1
            else:
                mismatches += 1
        
        alignment_length = len(aligned_seq1)
        
        # Calculate normalized weights
        w_match = 1.0
        w_mismatch = self.mismatch_score / self.match_score if self.match_score != 0 else -0.5
        w_gap = self.gap_score / self.match_score if self.match_score != 0 else -1.0
        
        # Calculate weighted sum
        weighted_sum = (w_match * matches + 
                       w_mismatch * mismatches + 
                       w_gap * gaps)
        
        # Normalize by length
        wsi = weighted_sum / alignment_length if alignment_length > 0 else 0.0
        
        return wsi
    
    def conservation_score(self, aligned_seq1, aligned_seq2):
        """
        Scoring Method 3: Conservation Score (CS)
        
        Formula: CS = (M - α*MM - β*G_o - γ*G_e) / L
        
        Where:
        - M: Number of matches
        - MM: Number of mismatches
        - G_o: Number of gap openings (consecutive gaps counted once)
        - G_e: Number of gap extensions (additional gaps in runs)
        - L: Alignment length
        - α, β, γ: Penalty factors
        
        Penalty factors:
        - α = 0.5 (mismatches are moderately penalized)
        - β = 1.0 (gap openings heavily penalized)
        - γ = 0.3 (gap extensions lightly penalized)
        
        Range: [0, 1] where 1 indicates perfect conservation
        
        This score reflects biological conservation, recognizing that gap
        openings are more disruptive than gap extensions, and that runs of
        matches indicate conserved functional regions.
        """
        matches = 0
        mismatches = 0
        gap_openings = 0
        gap_extensions = 0
        
        in_gap_seq1 = False
        in_gap_seq2 = False
        
        for i, (a, b) in enumerate(zip(aligned_seq1, aligned_seq2)):
            if a == '-' or b == '-':
                # Check if this is a gap opening or extension
                if a == '-':
                    if not in_gap_seq1:
                        gap_openings += 1
                        in_gap_seq1 = True
                    else:
                        gap_extensions += 1
                else:
                    in_gap_seq1 = False
                
                if b == '-':
                    if not in_gap_seq2:
                        gap_openings += 1
                        in_gap_seq2 = True
                    else:
                        gap_extensions += 1
                else:
                    in_gap_seq2 = False
            else:
                in_gap_seq1 = False
                in_gap_seq2 = False
                
                if a == b:
                    matches += 1
                else:
                    mismatches += 1
        
        alignment_length = len(aligned_seq1)
        
        # Penalty factors
        alpha = 0.5   # Mismatch penalty
        beta = 1.0    # Gap opening penalty
        gamma = 0.3   # Gap extension penalty
        
        # Calculate conservation score
        numerator = matches - (alpha * mismatches) - (beta * gap_openings) - (gamma * gap_extensions)
        cs = numerator / alignment_length if alignment_length > 0 else 0.0
        
        # Normalize to [0, 1] range
        cs = max(0.0, min(1.0, cs))
        
        return cs
    
    def _calculate_raw_score(self, aligned_seq1, aligned_seq2):
        """Calculate the raw alignment score."""
        score = 0
        for a, b in zip(aligned_seq1, aligned_seq2):
            if a == '-' or b == '-':
                score += self.gap_score
            elif a == b:
                score += self.match_score
            else:
                score += self.mismatch_score
        return score
    
    def _get_alignment_statistics(self, aligned_seq1, aligned_seq2):
        """Get detailed alignment statistics."""
        matches = 0
        mismatches = 0
        gaps = 0
        gap_openings = 0
        
        in_gap = False
        
        for a, b in zip(aligned_seq1, aligned_seq2):
            if a == '-' or b == '-':
                gaps += 1
                if not in_gap:
                    gap_openings += 1
                    in_gap = True
            else:
                in_gap = False
                if a == b:
                    matches += 1
                else:
                    mismatches += 1
        
        alignment_length = len(aligned_seq1)
        percent_identity = (matches / alignment_length * 100) if alignment_length > 0 else 0
        
        return {
            'length': alignment_length,
            'matches': matches,
            'mismatches': mismatches,
            'gaps': gaps,
            'gap_openings': gap_openings,
            'percent_identity': percent_identity
        }
    
    def _interpret_scores(self, nas, wsi, cs):
        """Provide interpretation of the scores."""
        interpretations = []
        
        # Interpret NAS
        if nas >= 0.9:
            interpretations.append("NAS: Excellent similarity (≥90%)")
        elif nas >= 0.7:
            interpretations.append("NAS: Good similarity (70-90%)")
        elif nas >= 0.5:
            interpretations.append("NAS: Moderate similarity (50-70%)")
        else:
            interpretations.append("NAS: Low similarity (<50%)")
        
        # Interpret WSI
        if wsi >= 0.7:
            interpretations.append("WSI: Strong weighted similarity")
        elif wsi >= 0.4:
            interpretations.append("WSI: Moderate weighted similarity")
        elif wsi >= 0.0:
            interpretations.append("WSI: Weak weighted similarity")
        else:
            interpretations.append("WSI: Poor weighted similarity")
        
        # Interpret CS
        if cs >= 0.8:
            interpretations.append("CS: Highly conserved sequences")
        elif cs >= 0.6:
            interpretations.append("CS: Well conserved sequences")
        elif cs >= 0.4:
            interpretations.append("CS: Moderately conserved sequences")
        else:
            interpretations.append("CS: Poorly conserved sequences")
        
        return interpretations
    
    def compare_scoring_methods(self, alignments):
        """
        Compare multiple alignments using all three scoring methods.
        
        Args:
            alignments: List of tuples (seq1, seq2, score, description)
        
        Returns:
            Comparison report
        """
        results = []
        
        for seq1, seq2, score, desc in alignments:
            scores = self.calculate_all_scores(seq1, seq2, score)
            results.append({
                'description': desc,
                'scores': scores
            })
        
        return results
    
    def print_score_explanation(self):
        """Print detailed explanation of all three scoring methods."""
        explanation = """
        ╔══════════════════════════════════════════════════════════════════════╗
        ║           SEQUENCE SIMILARITY SCORING METHODS                        ║
        ╚══════════════════════════════════════════════════════════════════════╝
        
        Method 1: NORMALIZED ALIGNMENT SCORE (NAS)
        ────────────────────────────────────────────────────────────────────
        Formula: NAS = (Actual_Score - Min_Score) / (Max_Score - Min_Score)
        
        • Normalizes raw alignment score to [0, 1] range
        • Accounts for sequence length and scoring parameters
        • 1.0 = perfect match, 0.0 = worst possible alignment
        • Best for: Comparing alignments with different lengths
        
        
        Method 2: WEIGHTED SIMILARITY INDEX (WSI)
        ────────────────────────────────────────────────────────────────────
        Formula: WSI = (w_m*Matches + w_mm*Mismatches + w_g*Gaps) / Length
        
        • Weights based on scoring scheme (match, mismatch, gap)
        • Range: [-1, 1] where 1 = perfect similarity
        • Reflects relative importance of different alignment features
        • Best for: Understanding alignment composition
        
        
        Method 3: CONSERVATION SCORE (CS)
        ────────────────────────────────────────────────────────────────────
        Formula: CS = (M - 0.5*MM - 1.0*G_open - 0.3*G_ext) / Length
        
        • Biological conservation metric
        • Distinguishes gap openings from extensions
        • Emphasizes consecutive matches (conserved regions)
        • Range: [0, 1] where 1 = perfectly conserved
        • Best for: Identifying functionally conserved regions
        
        
        INTERPRETATION GUIDELINES:
        ────────────────────────────────────────────────────────────────────
        All three scores together provide comprehensive similarity assessment:
        
        • High NAS + High WSI + High CS → Strong overall similarity
        • High NAS + Low CS → Similar but with many gaps
        • High CS + Lower NAS → Well conserved functional regions
        • Low all three → Sequences are dissimilar
        
        ═══════════════════════════════════════════════════════════════════════
        """
        print(explanation)


def format_score_report(score_dict):
    """Format a comprehensive score report."""
    report = []
    report.append("\n" + "="*70)
    report.append("SIMILARITY SCORING ANALYSIS")
    report.append("="*70)
    
    # Display scores
    report.append("\n▶ SIMILARITY SCORES:")
    report.append(f"  • Normalized Alignment Score (NAS): {score_dict['NAS']:.4f}")
    report.append(f"  • Weighted Similarity Index (WSI):  {score_dict['WSI']:.4f}")
    report.append(f"  • Conservation Score (CS):          {score_dict['CS']:.4f}")
    report.append(f"  • Raw Alignment Score:              {score_dict['raw_score']:.2f}")
    
    # Display statistics
    stats = score_dict['statistics']
    report.append("\n▶ ALIGNMENT STATISTICS:")
    report.append(f"  • Length:          {stats['length']} bp")
    report.append(f"  • Matches:         {stats['matches']} ({stats['percent_identity']:.1f}%)")
    report.append(f"  • Mismatches:      {stats['mismatches']}")
    report.append(f"  • Gaps:            {stats['gaps']}")
    report.append(f"  • Gap Openings:    {stats['gap_openings']}")
    
    # Display interpretation
    report.append("\n▶ INTERPRETATION:")
    for interp in score_dict['interpretation']:
        report.append(f"  • {interp}")
    
    report.append("\n" + "="*70)
    
    return "\n".join(report)


# Example usage and testing
if __name__ == "__main__":
    # Create scorer
    scorer = SequenceSimilarityScorer(match_score=2, mismatch_score=-1, gap_score=-2)
    
    # Print explanation
    scorer.print_score_explanation()
    
    # Test alignments
    print("\n\n" + "="*70)
    print("EXAMPLE ALIGNMENTS")
    print("="*70)
    
    test_cases = [
        ("ACCGTGAAGCCAATAC", "AGCGTGCAGCCAATAC", "Perfect match regions"),
        ("ACCGT-GAAGCCAATAC", "AGCGTGC-AGCCAATAC", "Alignment with gaps"),
        ("AAAAAAAAAA", "TTTTTTTTTT", "Complete mismatch"),
        ("ACGTACGTACGT", "ACGTACGTACGT", "Perfect alignment"),
        ("ACGT--ACGT", "ACGTACACGT", "Gap opening and extension"),
    ]
    
    for seq1, seq2, desc in test_cases:
        print(f"\n▼ Test Case: {desc}")
        print(f"   Seq1: {seq1}")
        print(f"   Seq2: {seq2}")
        
        scores = scorer.calculate_all_scores(seq1, seq2)
        print(format_score_report(scores))
        print()
