"""
Chromatogram file reader for various formats including ABI, SCF, and others.
"""

import os
from typing import Dict, Optional, Any
from Bio import SeqIO
import numpy as np


class ChromatogramReader:
    """
    Reader class for Sanger chromatogram files.
    
    Supports common formats like ABI (.ab1), SCF, and others.
    """
    
    def __init__(self, file_path: str):
        """
        Initialize the chromatogram reader.
        
        Args:
            file_path: Path to the chromatogram file
        """
        self.file_path = file_path
        self.format = self._detect_format()
        self._data = None
        
    def _detect_format(self) -> str:
        """Detect the file format based on extension."""
        ext = os.path.splitext(self.file_path)[1].lower()
        if ext in ['.ab1', '.abi']:
            return 'abi'
        elif ext in ['.scf']:
            return 'scf'
        else:
            raise ValueError(f"Unsupported file format: {ext}")
    
    def read(self) -> Dict[str, Any]:
        """
        Read and parse the chromatogram file.
        
        Returns:
            Dictionary containing chromatogram data including traces, sequence, and quality scores
        """
        if self._data is None:
            if self.format == 'abi':
                self._data = self._read_abi()
            elif self.format == 'scf':
                self._data = self._read_scf()
        
        return self._data
    
    def _read_abi(self) -> Dict[str, Any]:
        """Read ABI format chromatogram."""
        try:
            with open(self.file_path, 'rb') as handle:
                record = SeqIO.read(handle, 'abi')
                
            # Extract trace data
            traces = {
                'A': record.annotations.get('abif_raw', {}).get('DATA9', []),
                'T': record.annotations.get('abif_raw', {}).get('DATA10', []),
                'G': record.annotations.get('abif_raw', {}).get('DATA11', []),
                'C': record.annotations.get('abif_raw', {}).get('DATA12', [])
            }
            
            return {
                'sequence': str(record.seq),
                'quality_scores': record.letter_annotations.get('phred_quality', []),
                'traces': traces,
                'peak_positions': record.annotations.get('abif_raw', {}).get('PLOC1', []),
                'annotations': record.annotations
            }
            
        except Exception as e:
            raise ValueError(f"Error reading ABI file: {e}")
    
    def _read_scf(self) -> Dict[str, Any]:
        """Read SCF format chromatogram."""
        try:
            with open(self.file_path, 'rb') as handle:
                record = SeqIO.read(handle, 'scf')
                
            return {
                'sequence': str(record.seq),
                'quality_scores': record.letter_annotations.get('phred_quality', []),
                'traces': record.annotations.get('traces', {}),
                'peak_positions': record.annotations.get('peak_positions', []),
                'annotations': record.annotations
            }
            
        except Exception as e:
            raise ValueError(f"Error reading SCF file: {e}")
    
    @property
    def sequence(self) -> str:
        """Get the called sequence."""
        data = self.read()
        return data.get('sequence', '')
    
    @property
    def quality_scores(self) -> list:
        """Get quality scores for each base."""
        data = self.read()
        return data.get('quality_scores', [])
    
    @property
    def traces(self) -> Dict[str, list]:
        """Get the trace data for each base (A, T, G, C)."""
        data = self.read()
        return data.get('traces', {})