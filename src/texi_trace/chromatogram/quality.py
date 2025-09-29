"""
Quality assessment tools for chromatogram data.
"""

from typing import Dict, List, Tuple
import numpy as np


class QualityAssessment:
    """
    Tools for assessing chromatogram and sequence quality.
    """
    
    @staticmethod
    def assess_overall_quality(quality_scores: List[int]) -> Dict[str, float]:
        """
        Assess overall sequence quality metrics.
        
        Args:
            quality_scores: List of Phred quality scores
            
        Returns:
            Dictionary of quality metrics
        """
        if not quality_scores:
            return {
                'mean_quality': 0.0,
                'median_quality': 0.0,
                'q20_percentage': 0.0,
                'q30_percentage': 0.0,
                'low_quality_bases': 0
            }
        
        scores = np.array(quality_scores)
        
        # Basic statistics
        mean_quality = np.mean(scores)
        median_quality = np.median(scores)
        
        # Quality thresholds
        q20_count = np.sum(scores >= 20)
        q30_count = np.sum(scores >= 30)
        low_quality_count = np.sum(scores < 20)
        
        total_bases = len(scores)
        
        return {
            'mean_quality': float(mean_quality),
            'median_quality': float(median_quality),
            'q20_percentage': (q20_count / total_bases) * 100,
            'q30_percentage': (q30_count / total_bases) * 100,
            'low_quality_bases': int(low_quality_count),
            'total_bases': total_bases
        }
    
    @staticmethod
    def find_quality_regions(quality_scores: List[int], 
                           min_quality: int = 20, 
                           min_length: int = 10) -> List[Tuple[int, int]]:
        """
        Find continuous regions of good quality sequence.
        
        Args:
            quality_scores: List of Phred quality scores
            min_quality: Minimum quality threshold
            min_length: Minimum length of good quality region
            
        Returns:
            List of tuples (start, end) for good quality regions
        """
        if not quality_scores:
            return []
        
        good_regions = []
        current_start = None
        
        for i, score in enumerate(quality_scores):
            if score >= min_quality:
                if current_start is None:
                    current_start = i
            else:
                if current_start is not None:
                    region_length = i - current_start
                    if region_length >= min_length:
                        good_regions.append((current_start, i))
                    current_start = None
        
        # Handle region that extends to end of sequence
        if current_start is not None:
            region_length = len(quality_scores) - current_start
            if region_length >= min_length:
                good_regions.append((current_start, len(quality_scores)))
        
        return good_regions
    
    @staticmethod
    def trim_low_quality_ends(sequence: str, quality_scores: List[int], 
                            min_quality: int = 20) -> Tuple[str, List[int], int, int]:
        """
        Trim low quality bases from sequence ends.
        
        Args:
            sequence: DNA sequence string
            quality_scores: Corresponding quality scores
            min_quality: Minimum quality threshold
            
        Returns:
            Tuple of (trimmed_sequence, trimmed_scores, start_pos, end_pos)
        """
        if not sequence or not quality_scores:
            return sequence, quality_scores, 0, len(sequence)
        
        # Find first good quality base
        start_pos = 0
        for i, score in enumerate(quality_scores):
            if score >= min_quality:
                start_pos = i
                break
        
        # Find last good quality base
        end_pos = len(quality_scores)
        for i in range(len(quality_scores) - 1, -1, -1):
            if quality_scores[i] >= min_quality:
                end_pos = i + 1
                break
        
        # Trim sequence and scores
        trimmed_sequence = sequence[start_pos:end_pos]
        trimmed_scores = quality_scores[start_pos:end_pos]
        
        return trimmed_sequence, trimmed_scores, start_pos, end_pos
    
    @staticmethod
    def identify_problematic_regions(quality_scores: List[int], 
                                   traces: Dict[str, List] = None) -> List[Dict]:
        """
        Identify regions with potential sequencing problems.
        
        Args:
            quality_scores: List of quality scores
            traces: Optional trace data for additional analysis
            
        Returns:
            List of problematic regions with descriptions
        """
        problems = []
        
        if not quality_scores:
            return problems
        
        # Find consecutive low quality regions
        low_qual_start = None
        for i, score in enumerate(quality_scores):
            if score < 15:  # Very low quality
                if low_qual_start is None:
                    low_qual_start = i
            else:
                if low_qual_start is not None:
                    if i - low_qual_start >= 5:  # At least 5 consecutive low quality bases
                        problems.append({
                            'type': 'low_quality_region',
                            'start': low_qual_start,
                            'end': i,
                            'description': f'Low quality region ({i - low_qual_start} bases)',
                            'severity': 'high' if i - low_qual_start >= 10 else 'medium'
                        })
                    low_qual_start = None
        
        # Handle low quality region at end
        if low_qual_start is not None:
            length = len(quality_scores) - low_qual_start
            if length >= 5:
                problems.append({
                    'type': 'low_quality_region',
                    'start': low_qual_start,
                    'end': len(quality_scores),
                    'description': f'Low quality region at end ({length} bases)',
                    'severity': 'high' if length >= 10 else 'medium'
                })
        
        # Additional analysis with trace data
        if traces:
            # Look for regions with very weak signal
            for base, trace in traces.items():
                if trace:
                    trace_array = np.array(trace)
                    weak_signal_threshold = np.max(trace_array) * 0.1
                    
                    # Find regions where all traces are weak
                    weak_regions = []
                    for i in range(0, len(trace_array), 10):  # Check every 10 positions
                        window = trace_array[i:i+10]
                        if np.max(window) < weak_signal_threshold:
                            weak_regions.append(i)
                    
                    if weak_regions:
                        problems.append({
                            'type': 'weak_signal',
                            'positions': weak_regions,
                            'description': f'Weak signal detected in {base} trace',
                            'severity': 'medium'
                        })
        
        return problems