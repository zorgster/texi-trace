"""
texi-trace: A tool for scrutinizing Sanger chromatograms for mtDNA analysis.

This package provides functionality to:
- Parse and analyze Sanger chromatogram files
- Align sequences to reference genomes
- Visualize chromatogram data and alignment results
"""

__version__ = "0.1.0"
__author__ = "zorgster"

# Import core modules with graceful handling of missing dependencies
try:
    from .chromatogram import ChromatogramReader
    _CHROMATOGRAM_AVAILABLE = True
except ImportError:
    ChromatogramReader = None
    _CHROMATOGRAM_AVAILABLE = False

try:
    from .alignment import SequenceAligner
    _ALIGNMENT_AVAILABLE = True
except ImportError:
    SequenceAligner = None
    _ALIGNMENT_AVAILABLE = False

try:
    from .visualization import ChromatogramVisualizer
    _VISUALIZATION_AVAILABLE = True
except ImportError:
    ChromatogramVisualizer = None
    _VISUALIZATION_AVAILABLE = False

__all__ = [
    "ChromatogramReader",
    "SequenceAligner", 
    "ChromatogramVisualizer"
]