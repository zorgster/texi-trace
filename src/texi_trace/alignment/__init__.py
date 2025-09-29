"""Sequence alignment module for mtDNA analysis."""

try:
    from .aligner import SequenceAligner
    _ALIGNER_AVAILABLE = True
except ImportError:
    SequenceAligner = None
    _ALIGNER_AVAILABLE = False

from .mtdna import MtDNAReference

__all__ = ["SequenceAligner", "MtDNAReference"]