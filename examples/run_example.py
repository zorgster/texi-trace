#!/usr/bin/env python3
"""
Simple example to demonstrate texi-trace functionality without full dependencies.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

def test_basic_functionality():
    """Test basic functionality that works without BioPython."""
    print("texi-trace Basic Functionality Demo")
    print("=" * 40)
    
    # Test sequence utilities
    print("\n1. Testing sequence utilities:")
    from texi_trace.utils.sequence_utils import reverse_complement, gc_content, translate_dna
    
    test_seq = "ATCGATCGATCG"
    rc_seq = reverse_complement(test_seq)
    gc_pct = gc_content(test_seq)
    
    print(f"   Original sequence: {test_seq}")
    print(f"   Reverse complement: {rc_seq}")
    print(f"   GC content: {gc_pct:.1f}%")
    
    # Test translation
    protein = translate_dna("ATGAAATAG")  # Met-Lys-Stop
    print(f"   Translation of ATG AAA TAG: {protein}")
    
    # Test mtDNA reference
    print("\n2. Testing mtDNA reference:")
    from texi_trace.alignment.mtdna import MtDNAReference
    
    ref = MtDNAReference()
    print(f"   Reference length: {len(ref.reference)} bases")
    print(f"   First 50 bases: {ref.reference[:50]}")
    
    # Test gene regions
    gene_regions = ref.get_gene_regions()
    print(f"   Number of gene regions: {len(gene_regions)}")
    for gene_name, (start, end) in list(gene_regions.items())[:3]:
        print(f"   {gene_name}: {start}-{end}")
    
    # Test variant detection
    print("\n3. Testing variant detection:")
    # Create a sequence with some differences
    modified_seq = ref.reference[:30] + "T" + ref.reference[31:60]
    variants = ref.find_variants(modified_seq, 0)
    print(f"   Found {len(variants)} variant(s)")
    if variants:
        for var in variants:
            print(f"   Position {var['position']}: {var['reference']} -> {var['alternate']}")
    
    # Test file utilities
    print("\n4. Testing file utilities:")
    from texi_trace.utils.file_utils import validate_file_format, create_output_directory
    
    print(f"   'test.ab1' is valid ABI file: {validate_file_format('test.ab1', ['.ab1', '.scf'])}")
    print(f"   'test.txt' is valid ABI file: {validate_file_format('test.txt', ['.ab1', '.scf'])}")
    
    # Create a test output directory
    output_dir = create_output_directory("example_output")
    print(f"   Created output directory: {output_dir}")
    
    print("\n✓ All basic functionality tests completed successfully!")
    
    return True

def demonstrate_analysis_workflow():
    """Demonstrate the analysis workflow with mock data."""
    print("\n" + "=" * 50)
    print("Mock Analysis Workflow Demo")
    print("=" * 50)
    
    from texi_trace.alignment.mtdna import MtDNAReference
    from texi_trace.utils.sequence_utils import gc_content
    
    # Simulate reading a chromatogram (normally from .ab1 file)
    print("\n1. Simulating chromatogram data:")
    ref = MtDNAReference()
    mock_sequence = ref.reference[:100]  # Use part of reference as mock sequence
    mock_quality_scores = [30, 35, 28, 32, 29, 31, 33, 27, 36, 30] * 10  # Mock quality scores
    
    print(f"   Mock sequence length: {len(mock_sequence)} bases")
    print(f"   Average quality score: {sum(mock_quality_scores)/len(mock_quality_scores):.1f}")
    print(f"   GC content: {gc_content(mock_sequence):.1f}%")
    
    # Quality assessment
    print("\n2. Quality assessment:")
    high_quality_bases = sum(1 for q in mock_quality_scores if q >= 30)
    print(f"   High quality bases (Q≥30): {high_quality_bases}/{len(mock_quality_scores)} ({high_quality_bases/len(mock_quality_scores)*100:.1f}%)")
    
    # Alignment simulation
    print("\n3. Alignment analysis:")
    # Compare sequence to reference (should be identical in this mock case)
    alignment_identity = 100.0  # Perfect match since we used reference sequence
    print(f"   Alignment identity: {alignment_identity:.1f}%")
    
    # Gene region identification
    print("\n4. Gene region analysis:")
    gene = ref.identify_gene_region(50)  # Check what gene is at position 50
    print(f"   Position 50 is in gene region: {gene if gene else 'intergenic'}")
    
    # Variant calling simulation
    print("\n5. Variant analysis:")
    # Create a mock variant by changing one base
    variant_seq = mock_sequence[:30] + "A" + mock_sequence[31:]
    variants = ref.find_variants(variant_seq, 0)
    print(f"   Detected {len(variants)} variant(s)")
    
    print("\n✓ Mock analysis workflow completed!")

def show_installation_status():
    """Show which components are available based on installed dependencies."""
    print("\n" + "=" * 50)
    print("Installation Status")
    print("=" * 50)
    
    # Check core dependencies
    modules_to_check = [
        ("NumPy", "numpy"),
        ("Matplotlib", "matplotlib"),
        ("BioPython", "Bio"),
        ("Pandas", "pandas"),
        ("SciPy", "scipy"),
        ("Click", "click")
    ]
    
    available_modules = []
    missing_modules = []
    
    for name, module in modules_to_check:
        try:
            __import__(module)
            available_modules.append(name)
        except ImportError:
            missing_modules.append(name)
    
    print(f"\n✓ Available modules ({len(available_modules)}):")
    for module in available_modules:
        print(f"   - {module}")
    
    if missing_modules:
        print(f"\n⚠ Missing modules ({len(missing_modules)}):")
        for module in missing_modules:
            print(f"   - {module}")
        
        print(f"\nTo install missing dependencies:")
        print(f"   pip install -r requirements.txt")
    else:
        print(f"\n✓ All dependencies are installed!")
    
    # Show feature availability
    print(f"\nFeature availability:")
    print(f"   - Sequence utilities: ✓ Available")
    print(f"   - mtDNA reference: ✓ Available")
    print(f"   - File utilities: ✓ Available")
    print(f"   - Chromatogram reading: {'✓ Available' if 'BioPython' in available_modules else '⚠ Requires BioPython'}")
    print(f"   - Sequence alignment: {'✓ Available' if 'BioPython' in available_modules else '⚠ Requires BioPython'}")
    print(f"   - Visualization: {'✓ Available' if 'Matplotlib' in available_modules else '⚠ Requires Matplotlib'}")

if __name__ == "__main__":
    try:
        show_installation_status()
        test_basic_functionality()
        demonstrate_analysis_workflow()
        
        print("\n" + "=" * 50)
        print("Demo completed successfully!")
        print("=" * 50)
        print("\nNext steps:")
        print("1. Install missing dependencies: pip install -r requirements.txt")
        print("2. Place chromatogram files in data/chromatograms/")
        print("3. Run: texi-trace analyze your_file.ab1")
        print("4. Check the examples/ directory for more detailed usage")
        
    except Exception as e:
        print(f"\n✗ Error during demo: {e}")
        import traceback
        traceback.print_exc()