"""
Sequence alignment functionality for comparing chromatogram sequences 
to reference genomes, specifically optimized for mtDNA analysis.
"""

from typing import Dict, List, Tuple, Optional
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np


class SequenceAligner:
    """
    Aligner for comparing Sanger sequencing results to reference genomes.
    Optimized for mitochondrial DNA analysis.
    """
    
    def __init__(self, reference_sequence: str, match_score: float = 2, 
                 mismatch_score: float = -1, gap_open: float = -2, gap_extend: float = -0.5):
        """
        Initialize the sequence aligner.
        
        Args:
            reference_sequence: The reference genome sequence
            match_score: Score for matching bases
            mismatch_score: Penalty for mismatched bases
            gap_open: Penalty for opening a gap
            gap_extend: Penalty for extending a gap
        """
        self.reference = reference_sequence.upper()
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_open = gap_open
        self.gap_extend = gap_extend
        
    def align(self, query_sequence: str, local: bool = False) -> Dict:
        """
        Perform sequence alignment.
        
        Args:
            query_sequence: The sequence to align to the reference
            local: If True, perform local alignment; otherwise global
            
        Returns:
            Dictionary containing alignment results
        """
        query = query_sequence.upper().strip()
        
        if local:
            alignments = pairwise2.align.localms(
                self.reference, query,
                self.match_score, self.mismatch_score,
                self.gap_open, self.gap_extend,
                one_alignment_only=True
            )
        else:
            alignments = pairwise2.align.globalms(
                self.reference, query,
                self.match_score, self.mismatch_score,
                self.gap_open, self.gap_extend,
                one_alignment_only=True
            )
        
        if not alignments:
            return {
                'aligned_reference': '',
                'aligned_query': '',
                'score': 0,
                'identity': 0,
                'gaps': 0,
                'mismatches': [],
                'start_pos': 0,
                'end_pos': 0
            }
        
        alignment = alignments[0]
        ref_aligned, query_aligned, score, start, end = alignment
        
        # Calculate alignment statistics
        stats = self._calculate_alignment_stats(ref_aligned, query_aligned, start, end)
        
        return {
            'aligned_reference': ref_aligned,
            'aligned_query': query_aligned,
            'score': score,
            'start_pos': start,
            'end_pos': end,
            **stats
        }
    
    def _calculate_alignment_stats(self, ref_aligned: str, query_aligned: str, 
                                 start: int, end: int) -> Dict:
        """Calculate alignment statistics."""
        matches = 0
        mismatches = []
        gaps = 0
        
        for i, (r, q) in enumerate(zip(ref_aligned, query_aligned)):
            if r == '-' or q == '-':
                gaps += 1
            elif r == q:
                matches += 1
            else:
                mismatches.append({
                    'position': start + i,
                    'reference': r,
                    'query': q
                })
        
        total_length = len([c for c in ref_aligned if c != '-'])
        identity = (matches / total_length * 100) if total_length > 0 else 0
        
        return {
            'identity': identity,
            'matches': matches,
            'gaps': gaps,
            'mismatches': mismatches,
            'alignment_length': len(ref_aligned)
        }
    
    def find_best_region(self, query_sequence: str, window_size: int = 100) -> Dict:
        """
        Find the best matching region in the reference for the query sequence.
        
        Args:
            query_sequence: Query sequence to find best match for
            window_size: Size of sliding window for comparison
            
        Returns:
            Dictionary with best matching region information
        """
        query = query_sequence.upper().strip()
        best_score = float('-inf')
        best_position = 0
        best_alignment = None
        
        # Slide window across reference sequence
        for i in range(0, len(self.reference) - len(query) + 1, window_size // 2):
            ref_window = self.reference[i:i + len(query) + window_size]
            
            # Perform local alignment within this window
            alignments = pairwise2.align.localms(
                ref_window, query,
                self.match_score, self.mismatch_score,
                self.gap_open, self.gap_extend,
                one_alignment_only=True
            )
            
            if alignments and alignments[0][2] > best_score:
                best_score = alignments[0][2]
                best_position = i
                best_alignment = alignments[0]
        
        if best_alignment:
            ref_aligned, query_aligned, score, start, end = best_alignment
            stats = self._calculate_alignment_stats(ref_aligned, query_aligned, 
                                                  best_position + start, best_position + end)
            
            return {
                'reference_start': best_position + start,
                'reference_end': best_position + end,
                'aligned_reference': ref_aligned,
                'aligned_query': query_aligned,
                'score': score,
                **stats
            }
        
        return {
            'reference_start': 0,
            'reference_end': 0,
            'aligned_reference': '',
            'aligned_query': '',
            'score': 0,
            'identity': 0,
            'gaps': 0,
            'mismatches': []
        }