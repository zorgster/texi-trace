"""
Additional chromatogram parsing utilities.
"""

from typing import Dict, List, Tuple
import numpy as np


class ChromatogramParser:
    """
    Additional parsing utilities for chromatogram data processing.
    """
    
    @staticmethod
    def smooth_traces(traces: Dict[str, List], window_size: int = 3) -> Dict[str, List]:
        """
        Apply smoothing to trace data to reduce noise.
        
        Args:
            traces: Dictionary of trace data for each base
            window_size: Size of smoothing window
            
        Returns:
            Dictionary of smoothed trace data
        """
        smoothed = {}
        
        for base, trace in traces.items():
            if not trace:
                smoothed[base] = trace
                continue
                
            # Apply moving average smoothing
            trace_array = np.array(trace)
            if len(trace_array) < window_size:
                smoothed[base] = trace
                continue
                
            # Pad the array for edge handling
            padded = np.pad(trace_array, (window_size//2, window_size//2), mode='edge')
            
            # Apply moving average
            kernel = np.ones(window_size) / window_size
            smoothed_trace = np.convolve(padded, kernel, mode='valid')
            
            smoothed[base] = smoothed_trace.tolist()
        
        return smoothed
    
    @staticmethod
    def find_peaks(trace: List[float], height: float = None, distance: int = 10) -> List[int]:
        """
        Find peaks in a single trace.
        
        Args:
            trace: Single trace data
            height: Minimum peak height
            distance: Minimum distance between peaks
            
        Returns:
            List of peak positions
        """
        if not trace or len(trace) < 3:
            return []
        
        trace_array = np.array(trace)
        
        # Set default height to 20% of max signal
        if height is None:
            height = np.max(trace_array) * 0.2
        
        peaks = []
        
        for i in range(1, len(trace_array) - 1):
            # Check if this is a local maximum above threshold
            if (trace_array[i] > trace_array[i-1] and 
                trace_array[i] > trace_array[i+1] and 
                trace_array[i] >= height):
                
                # Check distance constraint
                if not peaks or i - peaks[-1] >= distance:
                    peaks.append(i)
        
        return peaks
    
    @staticmethod
    def calculate_signal_to_noise(traces: Dict[str, List]) -> Dict[str, float]:
        """
        Calculate signal-to-noise ratio for each base trace.
        
        Args:
            traces: Dictionary of trace data
            
        Returns:
            Dictionary of SNR values for each base
        """
        snr_values = {}
        
        for base, trace in traces.items():
            if not trace:
                snr_values[base] = 0.0
                continue
            
            trace_array = np.array(trace)
            
            # Calculate signal (mean of top 10% values)
            sorted_values = np.sort(trace_array)
            top_10_percent = sorted_values[int(0.9 * len(sorted_values)):]
            signal = np.mean(top_10_percent)
            
            # Calculate noise (standard deviation of bottom 50% values)
            bottom_50_percent = sorted_values[:int(0.5 * len(sorted_values))]
            noise = np.std(bottom_50_percent)
            
            # Calculate SNR
            if noise > 0:
                snr_values[base] = signal / noise
            else:
                snr_values[base] = float('inf')
        
        return snr_values