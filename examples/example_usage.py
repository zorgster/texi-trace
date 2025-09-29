#!/usr/bin/env python3
"""
Example usage of texi-trace for chromatogram analysis.

This script demonstrates how to use the texi-trace library programmatically
to analyze Sanger chromatogram data and perform mtDNA sequence alignment.
"""

import os
from pathlib import Path

# Import texi-trace modules
from texi_trace.chromatogram import ChromatogramReader
from texi_trace.alignment import SequenceAligner
from texi_trace.visualization import ChromatogramVisualizer


def analyze_single_chromatogram(chromatogram_path: str, output_dir: str = "output"):
    """
    Example: Analyze a single chromatogram file.
    
    Args:
        chromatogram_path: Path to chromatogram file (.ab1, .scf)
        output_dir: Directory to save results
    """
    print(f"Analyzing chromatogram: {chromatogram_path}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Read chromatogram data
    reader = ChromatogramReader(chromatogram_path)
    data = reader.read()
    
    print(f"Sequence length: {len(data['sequence'])} bases")
    print(f"Sequence: {data['sequence'][:50]}..." if len(data['sequence']) > 50 else f"Sequence: {data['sequence']}")
    
    if data['quality_scores']:
        avg_quality = sum(data['quality_scores']) / len(data['quality_scores'])
        print(f"Average quality score: {avg_quality:.2f}")
    
    # Create visualizations
    visualizer = ChromatogramVisualizer()
    
    # Plot chromatogram traces
    if data['traces']:
        fig = visualizer.plot_traces(
            traces=data['traces'],
            sequence=data['sequence'],
            peak_positions=data.get('peak_positions', [])
        )
        fig.savefig(os.path.join(output_dir, "chromatogram_traces.png"), 
                   dpi=300, bbox_inches='tight')
        print(f"Trace plot saved to: {output_dir}/chromatogram_traces.png")
    
    # Plot quality scores
    if data['quality_scores']:
        fig = visualizer.plot_quality_scores(
            quality_scores=data['quality_scores'],
            sequence=data['sequence']
        )
        fig.savefig(os.path.join(output_dir, "quality_scores.png"),
                   dpi=300, bbox_inches='tight')
        print(f"Quality plot saved to: {output_dir}/quality_scores.png")
    
    return data


def align_to_reference(sequence: str, reference_path: str):
    """
    Example: Align sequence to a reference genome.
    
    Args:
        sequence: Query sequence to align
        reference_path: Path to reference genome file (FASTA)
    """
    print(f"Aligning sequence to reference: {reference_path}")
    
    # Read reference sequence
    with open(reference_path, 'r') as f:
        lines = f.readlines()
        reference_seq = ''.join(line.strip() for line in lines if not line.startswith('>'))
    
    print(f"Reference length: {len(reference_seq)} bases")
    
    # Perform alignment
    aligner = SequenceAligner(reference_seq)
    
    # Try both global and local alignment
    global_result = aligner.align(sequence, local=False)
    local_result = aligner.align(sequence, local=True)
    
    print("\nGlobal alignment results:")
    print(f"  Score: {global_result['score']:.2f}")
    print(f"  Identity: {global_result['identity']:.2f}%")
    print(f"  Mismatches: {len(global_result['mismatches'])}")
    print(f"  Gaps: {global_result['gaps']}")
    
    print("\nLocal alignment results:")
    print(f"  Score: {local_result['score']:.2f}")
    print(f"  Identity: {local_result['identity']:.2f}%")
    print(f"  Mismatches: {len(local_result['mismatches'])}")
    print(f"  Gaps: {local_result['gaps']}")
    
    # Show first few mismatches
    if local_result['mismatches']:
        print("\nFirst 5 mismatches:")
        for i, mm in enumerate(local_result['mismatches'][:5]):
            print(f"  Position {mm['position']}: {mm['reference']} -> {mm['query']}")
    
    return local_result


def batch_analysis_example():
    """
    Example: Batch processing of multiple chromatogram files.
    """
    print("Batch analysis example")
    
    # This would work if you have actual chromatogram files
    data_dir = Path("data/chromatograms")
    
    if not data_dir.exists():
        print(f"Data directory {data_dir} not found. Create it and add .ab1/.scf files to test batch processing.")
        return
    
    # Find all chromatogram files
    chromatogram_files = list(data_dir.glob("*.ab1")) + list(data_dir.glob("*.scf"))
    
    if not chromatogram_files:
        print("No chromatogram files found in data directory.")
        return
    
    print(f"Found {len(chromatogram_files)} chromatogram files")
    
    results = []
    for file_path in chromatogram_files:
        print(f"Processing: {file_path.name}")
        try:
            data = analyze_single_chromatogram(str(file_path), f"output/{file_path.stem}")
            results.append({
                'file': file_path.name,
                'length': len(data['sequence']),
                'avg_quality': sum(data['quality_scores'])/len(data['quality_scores']) if data['quality_scores'] else 0
            })
        except Exception as e:
            print(f"Error processing {file_path.name}: {e}")
    
    # Summary report
    if results:
        print("\nBatch Analysis Summary:")
        print("-" * 50)
        for result in results:
            print(f"{result['file']:30} {result['length']:4d} bp  Q={result['avg_quality']:5.1f}")


def create_sample_reference():
    """Create a sample mtDNA reference sequence for testing."""
    # This is a partial human mtDNA reference (first 100 bases of rCRS)
    sample_mtdna = """GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTG"""
    
    os.makedirs("data/references", exist_ok=True)
    ref_path = "data/references/sample_mtdna.fasta"
    
    with open(ref_path, 'w') as f:
        f.write(">Sample_human_mtDNA_partial\n")
        f.write(sample_mtdna + "\n")
    
    print(f"Sample reference created: {ref_path}")
    return ref_path


def main():
    """Main example function demonstrating texi-trace usage."""
    print("texi-trace Example Usage")
    print("=" * 50)
    
    # Create sample reference if it doesn't exist
    ref_path = "data/references/sample_mtdna.fasta"
    if not os.path.exists(ref_path):
        ref_path = create_sample_reference()
    
    # Example 1: Create a sample sequence for demonstration
    print("\n1. Sample sequence analysis:")
    sample_sequence = "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTG"
    print(f"Sample sequence: {sample_sequence}")
    
    # Example 2: Alignment demonstration
    print("\n2. Alignment to reference:")
    try:
        alignment_result = align_to_reference(sample_sequence, ref_path)
    except FileNotFoundError:
        print(f"Reference file not found: {ref_path}")
        print("Create a reference file to test alignment functionality.")
    
    # Example 3: Instructions for chromatogram analysis
    print("\n3. Chromatogram analysis:")
    print("To analyze actual chromatogram files:")
    print("- Place .ab1 or .scf files in data/chromatograms/")
    print("- Run: analyze_single_chromatogram('path/to/file.ab1')")
    print("- Or use CLI: texi-trace analyze path/to/file.ab1")
    
    # Example 4: Batch processing info
    print("\n4. Batch processing:")
    print("For batch processing, run: batch_analysis_example()")
    
    print("\nExample complete! Check the 'output' directory for generated files.")


if __name__ == "__main__":
    main()