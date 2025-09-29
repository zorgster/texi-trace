"""Tests for chromatogram reading functionality."""

import pytest
import tempfile
import os
from unittest.mock import patch, MagicMock

from texi_trace.chromatogram.reader import ChromatogramReader


class TestChromatogramReader:
    """Test cases for ChromatogramReader class."""
    
    def test_format_detection_abi(self):
        """Test ABI format detection."""
        with tempfile.NamedTemporaryFile(suffix='.ab1', delete=False) as tmp:
            tmp_path = tmp.name
        
        try:
            reader = ChromatogramReader(tmp_path)
            assert reader.format == 'abi'
        finally:
            os.unlink(tmp_path)
    
    def test_format_detection_scf(self):
        """Test SCF format detection."""
        with tempfile.NamedTemporaryFile(suffix='.scf', delete=False) as tmp:
            tmp_path = tmp.name
        
        try:
            reader = ChromatogramReader(tmp_path)
            assert reader.format == 'scf'
        finally:
            os.unlink(tmp_path)
    
    def test_unsupported_format(self):
        """Test unsupported format raises exception."""
        with tempfile.NamedTemporaryFile(suffix='.txt', delete=False) as tmp:
            tmp_path = tmp.name
        
        try:
            with pytest.raises(ValueError, match="Unsupported file format"):
                ChromatogramReader(tmp_path)
        finally:
            os.unlink(tmp_path)
    
    @patch('texi_trace.chromatogram.reader.SeqIO')
    def test_read_abi_success(self, mock_seqio):
        """Test successful ABI file reading."""
        # Mock SeqIO.read to return a mock record
        mock_record = MagicMock()
        mock_record.seq = "ATCG"
        mock_record.letter_annotations = {'phred_quality': [30, 25, 35, 20]}
        mock_record.annotations = {
            'abif_raw': {
                'DATA9': [100, 120, 110, 90],   # A trace
                'DATA10': [80, 90, 95, 100],    # T trace  
                'DATA11': [90, 85, 200, 95],    # G trace
                'DATA12': [200, 180, 90, 85],   # C trace
                'PLOC1': [10, 20, 30, 40]       # Peak positions
            }
        }
        mock_seqio.read.return_value = mock_record
        
        with tempfile.NamedTemporaryFile(suffix='.ab1', delete=False) as tmp:
            tmp_path = tmp.name
        
        try:
            reader = ChromatogramReader(tmp_path)
            data = reader.read()
            
            assert data['sequence'] == "ATCG"
            assert data['quality_scores'] == [30, 25, 35, 20]
            assert 'traces' in data
            assert 'A' in data['traces']
            assert data['peak_positions'] == [10, 20, 30, 40]
            
        finally:
            os.unlink(tmp_path)
    
    @patch('texi_trace.chromatogram.reader.SeqIO')
    def test_read_abi_error(self, mock_seqio):
        """Test ABI file reading error handling."""
        mock_seqio.read.side_effect = Exception("File corrupted")
        
        with tempfile.NamedTemporaryFile(suffix='.ab1', delete=False) as tmp:
            tmp_path = tmp.name
        
        try:
            reader = ChromatogramReader(tmp_path)
            with pytest.raises(ValueError, match="Error reading ABI file"):
                reader.read()
        finally:
            os.unlink(tmp_path)
    
    @patch('texi_trace.chromatogram.reader.SeqIO')
    def test_sequence_property(self, mock_seqio):
        """Test sequence property access."""
        mock_record = MagicMock()
        mock_record.seq = "ATCG"
        mock_record.letter_annotations = {'phred_quality': [30, 25, 35, 20]}
        mock_record.annotations = {'abif_raw': {}}
        mock_seqio.read.return_value = mock_record
        
        with tempfile.NamedTemporaryFile(suffix='.ab1', delete=False) as tmp:
            tmp_path = tmp.name
        
        try:
            reader = ChromatogramReader(tmp_path)
            assert reader.sequence == "ATCG"
        finally:
            os.unlink(tmp_path)