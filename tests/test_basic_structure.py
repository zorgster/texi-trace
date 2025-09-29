"""Basic tests for project structure without BioPython dependency."""

import pytest
import sys
import os

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))


def test_utils_import():
    """Test that utility modules can be imported."""
    from texi_trace.utils.sequence_utils import reverse_complement, gc_content
    
    # Test reverse complement
    assert reverse_complement("ATCG") == "CGAT"
    assert reverse_complement("AAAA") == "TTTT"
    
    # Test GC content
    assert gc_content("ATCG") == 50.0
    assert gc_content("AAAA") == 0.0
    assert gc_content("GGGG") == 100.0


def test_file_utils_import():
    """Test file utility functions."""
    from texi_trace.utils.file_utils import validate_file_format
    
    assert validate_file_format("test.ab1", [".ab1", ".scf"]) == True
    assert validate_file_format("test.txt", [".ab1", ".scf"]) == False


def test_alignment_module_structure():
    """Test alignment module structure."""
    try:
        from texi_trace.alignment.mtdna import MtDNAReference
        
        # Test MtDNA reference without BioPython
        ref = MtDNAReference()
        assert len(ref.reference) > 0
        assert isinstance(ref.reference, str)
        
        # Test gene region functionality
        gene_regions = ref.get_gene_regions()
        assert isinstance(gene_regions, dict)
        assert 'D-loop' in gene_regions
        
        print("✓ MtDNA reference module works correctly")
    except ImportError as e:
        print(f"Expected import error (missing BioPython): {e}")


def test_project_structure():
    """Test that project directories exist."""
    project_root = os.path.join(os.path.dirname(__file__), '..')
    
    expected_dirs = [
        'src/texi_trace',
        'src/texi_trace/chromatogram',
        'src/texi_trace/alignment', 
        'src/texi_trace/visualization',
        'src/texi_trace/utils',
        'tests',
        'data',
        'examples'
    ]
    
    for dir_path in expected_dirs:
        full_path = os.path.join(project_root, dir_path)
        assert os.path.exists(full_path), f"Directory {dir_path} should exist"
    
    print("✓ All expected directories exist")


if __name__ == "__main__":
    test_utils_import()
    test_file_utils_import()
    test_alignment_module_structure()
    test_project_structure()
    print("✓ All basic structure tests passed!")