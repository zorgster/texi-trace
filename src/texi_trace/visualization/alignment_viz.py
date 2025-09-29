"""
Visualization tools for sequence alignment results.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from typing import Dict, List, Tuple, Optional
import seaborn as sns


class AlignmentVisualizer:
    """
    Visualizer for sequence alignment results.
    """
    
    def __init__(self, style: str = 'seaborn-v0_8'):
        """
        Initialize the alignment visualizer.
        
        Args:
            style: Matplotlib style to use
        """
        try:
            plt.style.use(style)
        except:
            plt.style.use('default')
    
    def plot_alignment_overview(self, alignment_result: Dict, 
                              figsize: Tuple[int, int] = (15, 8)) -> plt.Figure:
        """
        Create an overview plot of the alignment.
        
        Args:
            alignment_result: Alignment result dictionary
            figsize: Figure size tuple
            
        Returns:
            Matplotlib figure object
        """
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, height_ratios=[1, 3])
        
        # Top panel: Alignment statistics
        self._plot_alignment_stats(ax1, alignment_result)
        
        # Bottom panel: Alignment visualization
        self._plot_alignment_detail(ax2, alignment_result)
        
        plt.tight_layout()
        return fig
    
    def _plot_alignment_stats(self, ax, alignment_result: Dict):
        """Plot alignment statistics in the top panel."""
        stats = [
            f"Score: {alignment_result.get('score', 0):.2f}",
            f"Identity: {alignment_result.get('identity', 0):.1f}%",
            f"Matches: {alignment_result.get('matches', 0)}",
            f"Mismatches: {len(alignment_result.get('mismatches', []))}",
            f"Gaps: {alignment_result.get('gaps', 0)}",
            f"Length: {alignment_result.get('alignment_length', 0)}"
        ]
        
        # Create text display
        ax.text(0.05, 0.5, "  |  ".join(stats), 
               fontsize=12, va='center', ha='left',
               bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue", alpha=0.7))
        
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        ax.set_title('Alignment Statistics', fontsize=14, fontweight='bold')
    
    def _plot_alignment_detail(self, ax, alignment_result: Dict):
        """Plot detailed alignment visualization."""
        ref_seq = alignment_result.get('aligned_reference', '')
        query_seq = alignment_result.get('aligned_query', '')
        
        if not ref_seq or not query_seq:
            ax.text(0.5, 0.5, 'No alignment data available', 
                   ha='center', va='center', fontsize=12)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            return
        
        # Determine display window (show first 100 positions if sequence is long)
        max_display = 100
        if len(ref_seq) > max_display:
            ref_display = ref_seq[:max_display]
            query_display = query_seq[:max_display]
            truncated = True
        else:
            ref_display = ref_seq
            query_display = query_seq
            truncated = False
        
        # Create color map for alignment
        colors = []
        for i, (r, q) in enumerate(zip(ref_display, query_display)):
            if r == '-' or q == '-':
                colors.append('lightcoral')  # Gap
            elif r == q:
                colors.append('lightgreen')  # Match
            else:
                colors.append('orange')      # Mismatch
        
        # Plot alignment bars
        positions = np.arange(len(ref_display))
        
        # Reference sequence bar
        for i, (pos, color) in enumerate(zip(positions, colors)):
            ax.barh(1, 1, left=pos, height=0.3, color=color, alpha=0.7)
            
        # Query sequence bar  
        for i, (pos, color) in enumerate(zip(positions, colors)):
            ax.barh(0, 1, left=pos, height=0.3, color=color, alpha=0.7)
        
        # Add sequence text
        for i, (pos, r, q) in enumerate(zip(positions, ref_display, query_display)):
            if i % 5 == 0 or len(ref_display) <= 50:  # Show every 5th base or all if short
                ax.text(pos + 0.5, 1.15, r, ha='center', va='center', 
                       fontsize=8, fontfamily='monospace')
                ax.text(pos + 0.5, -0.15, q, ha='center', va='center', 
                       fontsize=8, fontfamily='monospace')
        
        # Add labels and formatting
        ax.set_yticks([0.15, 1.15])
        ax.set_yticklabels(['Query', 'Reference'])
        ax.set_xlim(-0.5, len(ref_display) - 0.5)
        ax.set_ylim(-0.5, 1.5)
        ax.set_xlabel('Position')
        
        if truncated:
            ax.set_title(f'Alignment Detail (first {max_display} positions)', 
                        fontsize=12, fontweight='bold')
        else:
            ax.set_title('Alignment Detail', fontsize=12, fontweight='bold')
        
        # Add legend
        match_patch = mpatches.Patch(color='lightgreen', label='Match')
        mismatch_patch = mpatches.Patch(color='orange', label='Mismatch')
        gap_patch = mpatches.Patch(color='lightcoral', label='Gap')
        ax.legend(handles=[match_patch, mismatch_patch, gap_patch], 
                 loc='upper right', bbox_to_anchor=(1, 1))
    
    def plot_mismatch_distribution(self, alignment_result: Dict,
                                 figsize: Tuple[int, int] = (12, 6)) -> plt.Figure:
        """
        Plot distribution of mismatches along the alignment.
        
        Args:
            alignment_result: Alignment result dictionary
            figsize: Figure size tuple
            
        Returns:
            Matplotlib figure object
        """
        fig, ax = plt.subplots(figsize=figsize)
        
        mismatches = alignment_result.get('mismatches', [])
        
        if not mismatches:
            ax.text(0.5, 0.5, 'No mismatches found', 
                   ha='center', va='center', fontsize=12)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.set_title('Mismatch Distribution')
            return fig
        
        # Extract positions and create histogram
        positions = [mm['position'] for mm in mismatches]
        
        # Create bins for histogram
        alignment_length = alignment_result.get('alignment_length', max(positions) + 1)
        bins = np.linspace(0, alignment_length, min(50, alignment_length // 10 + 1))
        
        ax.hist(positions, bins=bins, alpha=0.7, color='orange', edgecolor='black')
        
        # Add individual mismatch markers
        ax.scatter(positions, [0.1] * len(positions), 
                  alpha=0.6, c='red', s=20, marker='|')
        
        ax.set_xlabel('Position')
        ax.set_ylabel('Number of Mismatches')
        ax.set_title('Distribution of Mismatches Along Alignment')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        return fig
    
    def plot_identity_sliding_window(self, alignment_result: Dict, 
                                   window_size: int = 50,
                                   figsize: Tuple[int, int] = (12, 6)) -> plt.Figure:
        """
        Plot sequence identity in sliding windows.
        
        Args:
            alignment_result: Alignment result dictionary
            window_size: Size of sliding window
            figsize: Figure size tuple
            
        Returns:
            Matplotlib figure object
        """
        fig, ax = plt.subplots(figsize=figsize)
        
        ref_seq = alignment_result.get('aligned_reference', '')
        query_seq = alignment_result.get('aligned_query', '')
        
        if not ref_seq or not query_seq:
            ax.text(0.5, 0.5, 'No alignment data available', 
                   ha='center', va='center', fontsize=12)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            return fig
        
        # Calculate identity in sliding windows
        window_centers = []
        identities = []
        
        for i in range(0, len(ref_seq) - window_size + 1, window_size // 4):
            ref_window = ref_seq[i:i + window_size]
            query_window = query_seq[i:i + window_size]
            
            matches = sum(1 for r, q in zip(ref_window, query_window) 
                         if r == q and r != '-' and q != '-')
            non_gap_positions = sum(1 for r, q in zip(ref_window, query_window)
                                  if r != '-' and q != '-')
            
            if non_gap_positions > 0:
                identity = (matches / non_gap_positions) * 100
                window_centers.append(i + window_size // 2)
                identities.append(identity)
        
        if identities:
            ax.plot(window_centers, identities, 'b-', linewidth=2, alpha=0.8)
            ax.fill_between(window_centers, identities, alpha=0.3, color='blue')
            
            # Add horizontal line for overall identity
            overall_identity = alignment_result.get('identity', 0)
            ax.axhline(y=overall_identity, color='red', linestyle='--', 
                      label=f'Overall Identity: {overall_identity:.1f}%')
        
        ax.set_xlabel('Position')
        ax.set_ylabel('Identity (%)')
        ax.set_title(f'Sequence Identity (sliding window size: {window_size})')
        ax.set_ylim(0, 100)
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        plt.tight_layout()
        return fig