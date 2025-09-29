"""
Visualization tools for chromatogram data including trace plots,
quality scores, and sequence alignment visualization.
"""

import matplotlib.pyplot as plt
import numpy as np
from typing import Dict, List, Optional, Tuple
import seaborn as sns


class ChromatogramVisualizer:
    """
    Visualizer for Sanger chromatogram data.
    """
    
    def __init__(self, style: str = 'seaborn-v0_8'):
        """
        Initialize the visualizer.
        
        Args:
            style: Matplotlib style to use
        """
        try:
            plt.style.use(style)
        except:
            plt.style.use('default')
        
        # Color scheme for bases
        self.base_colors = {
            'A': '#FF6B6B',  # Red
            'T': '#4ECDC4',  # Teal  
            'G': '#45B7D1',  # Blue
            'C': '#96CEB4'   # Green
        }
        
    def plot_traces(self, traces: Dict[str, List], sequence: str = '', 
                   peak_positions: List[int] = None, 
                   start_pos: int = 0, end_pos: int = None,
                   figsize: Tuple[int, int] = (15, 8)) -> plt.Figure:
        """
        Plot chromatogram traces for all four bases.
        
        Args:
            traces: Dictionary with trace data for each base (A, T, G, C)
            sequence: Called sequence to overlay on traces
            peak_positions: List of peak positions for base calls
            start_pos: Start position for plotting window
            end_pos: End position for plotting window
            figsize: Figure size tuple
            
        Returns:
            Matplotlib figure object
        """
        fig, ax = plt.subplots(figsize=figsize)
        
        # Determine plotting range
        max_length = max(len(trace) for trace in traces.values() if trace)
        if end_pos is None:
            end_pos = max_length
        
        end_pos = min(end_pos, max_length)
        start_pos = max(0, start_pos)
        
        # Plot traces for each base
        x_positions = np.arange(start_pos, end_pos)
        
        for base, trace in traces.items():
            if trace and len(trace) > start_pos:
                trace_segment = trace[start_pos:end_pos]
                ax.plot(x_positions[:len(trace_segment)], trace_segment, 
                       color=self.base_colors.get(base, 'black'), 
                       label=f'Base {base}', linewidth=1.5, alpha=0.8)
        
        # Add base calls if provided
        if sequence and peak_positions:
            self._add_base_calls(ax, sequence, peak_positions, start_pos, end_pos)
        
        ax.set_xlabel('Position')
        ax.set_ylabel('Signal Intensity')
        ax.set_title('Chromatogram Traces')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        return fig
    
    def _add_base_calls(self, ax, sequence: str, peak_positions: List[int], 
                       start_pos: int, end_pos: int):
        """Add base call annotations to the trace plot."""
        # Filter peak positions within the plotting range
        for i, pos in enumerate(peak_positions):
            if start_pos <= pos <= end_pos and i < len(sequence):
                base = sequence[i]
                color = self.base_colors.get(base, 'black')
                
                # Add vertical line at peak position
                ax.axvline(x=pos, color=color, linestyle='--', alpha=0.5)
                
                # Add base label
                ax.text(pos, ax.get_ylim()[1] * 0.9, base, 
                       ha='center', va='bottom', fontsize=10, 
                       color=color, fontweight='bold')
    
    def plot_quality_scores(self, quality_scores: List[int], sequence: str = '',
                          figsize: Tuple[int, int] = (15, 6)) -> plt.Figure:
        """
        Plot quality scores along the sequence.
        
        Args:
            quality_scores: List of Phred quality scores
            sequence: Corresponding sequence (optional)
            figsize: Figure size tuple
            
        Returns:
            Matplotlib figure object
        """
        fig, ax = plt.subplots(figsize=figsize)
        
        positions = range(len(quality_scores))
        
        # Color-code quality scores
        colors = []
        for score in quality_scores:
            if score >= 30:
                colors.append('#2ECC71')  # Good quality - green
            elif score >= 20:
                colors.append('#F39C12')  # Medium quality - orange  
            else:
                colors.append('#E74C3C')  # Poor quality - red
        
        bars = ax.bar(positions, quality_scores, color=colors, alpha=0.7)
        
        # Add quality thresholds
        ax.axhline(y=20, color='orange', linestyle='--', alpha=0.7, 
                  label='Q20 (99% accuracy)')
        ax.axhline(y=30, color='green', linestyle='--', alpha=0.7,
                  label='Q30 (99.9% accuracy)')
        
        ax.set_xlabel('Position')
        ax.set_ylabel('Phred Quality Score')
        ax.set_title('Sequence Quality Scores')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Add sequence bases as x-axis labels if sequence is short enough
        if sequence and len(sequence) <= 50:
            ax.set_xticks(positions[::max(1, len(sequence)//20)])
            ax.set_xticklabels([sequence[i] for i in positions[::max(1, len(sequence)//20)]])
        
        plt.tight_layout()
        return fig
    
    def plot_sequence_overview(self, traces: Dict[str, List], sequence: str,
                             quality_scores: List[int],
                             figsize: Tuple[int, int] = (15, 10)) -> plt.Figure:
        """
        Create a comprehensive overview plot with traces and quality scores.
        
        Args:
            traces: Trace data for each base
            sequence: Called sequence
            quality_scores: Quality scores for each base
            figsize: Figure size tuple
            
        Returns:
            Matplotlib figure object
        """
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, height_ratios=[3, 1])
        
        # Plot traces in upper panel
        max_length = max(len(trace) for trace in traces.values() if trace)
        x_positions = np.arange(max_length)
        
        for base, trace in traces.items():
            if trace:
                ax1.plot(x_positions[:len(trace)], trace, 
                        color=self.base_colors.get(base, 'black'),
                        label=f'Base {base}', linewidth=1.5, alpha=0.8)
        
        ax1.set_ylabel('Signal Intensity')
        ax1.set_title('Chromatogram Overview')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Plot quality scores in lower panel
        if quality_scores:
            positions = range(len(quality_scores))
            colors = ['#2ECC71' if q >= 30 else '#F39C12' if q >= 20 else '#E74C3C' 
                     for q in quality_scores]
            
            ax2.bar(positions, quality_scores, color=colors, alpha=0.7)
            ax2.axhline(y=20, color='orange', linestyle='--', alpha=0.7)
            ax2.axhline(y=30, color='green', linestyle='--', alpha=0.7)
            ax2.set_xlabel('Position')
            ax2.set_ylabel('Quality')
            ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        return fig