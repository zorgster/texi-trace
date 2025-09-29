# Usage Guide for texi-trace

## Quick Start

### Command Line Interface

The easiest way to use texi-trace is through the command line:

```bash
# Analyze a single chromatogram
texi-trace analyze sample.ab1

# Align to a reference genome  
texi-trace align sample.ab1 reference.fasta

# Batch process multiple files
texi-trace batch data/chromatograms/ --reference-file mtdna_ref.fasta
```

### Python API

For programmatic use:

```python
from texi_trace import ChromatogramReader, SequenceAligner, ChromatogramVisualizer

# Read chromatogram
reader = ChromatogramReader('sample.ab1')
data = reader.read()

# Analyze and visualize
visualizer = ChromatogramVisualizer()
fig = visualizer.plot_traces(data['traces'], data['sequence'])
fig.savefig('chromatogram.png')
```

## Detailed Usage

### 1. Chromatogram Analysis

#### Single File Analysis

```bash
# Basic analysis with default output
texi-trace analyze sample.ab1

# Specify output directory
texi-trace analyze sample.ab1 --output-dir results/

# Specify file format (if auto-detection fails)
texi-trace analyze sample.ab1 --format abi
```

#### Python API for Chromatogram Analysis

```python
from texi_trace.chromatogram import ChromatogramReader
from texi_trace.visualization import ChromatogramVisualizer
import os

# Read chromatogram file
reader = ChromatogramReader('data/chromatograms/sample.ab1')
data = reader.read()

# Access sequence information
print(f"Sequence: {data['sequence']}")
print(f"Length: {len(data['sequence'])} bases")
print(f"Quality scores: {data['quality_scores'][:10]}...")  # First 10

# Create visualizations
visualizer = ChromatogramVisualizer()

# Plot chromatogram traces
fig1 = visualizer.plot_traces(
    traces=data['traces'],
    sequence=data['sequence'],
    peak_positions=data.get('peak_positions', [])
)
fig1.savefig('output/traces.png', dpi=300)

# Plot quality scores
fig2 = visualizer.plot_quality_scores(
    quality_scores=data['quality_scores'],
    sequence=data['sequence']
)
fig2.savefig('output/quality.png', dpi=300)
```

### 2. Sequence Alignment

#### CLI Alignment

```bash
# Global alignment (default)
texi-trace align sample.ab1 reference.fasta --output-dir alignment_results/

# Local alignment
texi-trace align sample.ab1 reference.fasta --local-align --output-dir local_alignment/
```

#### Python API for Alignment

```python
from texi_trace.chromatogram import ChromatogramReader
from texi_trace.alignment import SequenceAligner, MtDNAReference
from texi_trace.visualization import AlignmentVisualizer

# Read chromatogram
reader = ChromatogramReader('sample.ab1')
query_sequence = reader.sequence

# Use built-in mtDNA reference
mtdna_ref = MtDNAReference()
aligner = SequenceAligner(mtdna_ref.reference)

# Perform alignment
alignment = aligner.align(query_sequence, local=False)

# Display results
print(f"Alignment score: {alignment['score']:.2f}")
print(f"Identity: {alignment['identity']:.2f}%")
print(f"Mismatches: {len(alignment['mismatches'])}")

# Visualize alignment
viz = AlignmentVisualizer()
fig = viz.plot_alignment_overview(alignment)
fig.savefig('alignment_overview.png', dpi=300)

# Find variants
variants = mtdna_ref.find_variants(query_sequence)
print(f"Found {len(variants)} variants")
```

### 3. Batch Processing

```bash
# Process all files in a directory
texi-trace batch data/chromatograms/

# With reference alignment
texi-trace batch data/chromatograms/ --reference-file data/references/human_mtdna.fasta

# Specify output location
texi-trace batch data/chromatograms/ --output-dir batch_results/
```

### 4. Quality Assessment

#### Using Quality Assessment Tools

```python
from texi_trace.chromatogram.quality import QualityAssessment
from texi_trace.chromatogram import ChromatogramReader

# Read chromatogram
reader = ChromatogramReader('sample.ab1')
data = reader.read()

# Assess quality
quality_stats = QualityAssessment.assess_overall_quality(data['quality_scores'])
print(f"Mean quality: {quality_stats['mean_quality']:.2f}")
print(f"Q30 percentage: {quality_stats['q30_percentage']:.1f}%")

# Find good quality regions
good_regions = QualityAssessment.find_quality_regions(
    data['quality_scores'], 
    min_quality=20, 
    min_length=50
)
print(f"Found {len(good_regions)} good quality regions")

# Trim low quality ends
trimmed_seq, trimmed_qual, start, end = QualityAssessment.trim_low_quality_ends(
    data['sequence'], 
    data['quality_scores']
)
print(f"Trimmed sequence: {len(trimmed_seq)} bases (was {len(data['sequence'])})")
```

### 5. mtDNA-Specific Analysis

#### Working with mtDNA References

```python
from texi_trace.alignment.mtdna import MtDNAReference

# Use built-in reference
mtdna = MtDNAReference()

# Or load custom reference
# mtdna = MtDNAReference.load_from_file('custom_mtdna.fasta')

# Get gene information
gene_regions = mtdna.get_gene_regions()
for gene, (start, end) in gene_regions.items():
    print(f"{gene}: {start}-{end}")

# Identify gene at specific position
gene = mtdna.identify_gene_region(1000)
print(f"Position 1000 is in: {gene}")

# Extract gene sequences
cox1_region = mtdna.get_region(5904, 7445)  # COX1 gene
print(f"COX1 length: {len(cox1_region)} bases")

# Translate to amino acids
protein = mtdna.translate_sequence(cox1_region)
print(f"COX1 protein: {protein[:50]}...")  # First 50 amino acids
```

### 6. Visualization Options

#### Chromatogram Visualization

```python
from texi_trace.visualization import ChromatogramVisualizer

viz = ChromatogramVisualizer()

# Basic trace plot
fig1 = viz.plot_traces(traces, sequence)

# Quality score plot
fig2 = viz.plot_quality_scores(quality_scores, sequence)

# Comprehensive overview
fig3 = viz.plot_sequence_overview(traces, sequence, quality_scores)

# Customized plot with specific region
fig4 = viz.plot_traces(traces, sequence, start_pos=100, end_pos=200)
```

#### Alignment Visualization

```python
from texi_trace.visualization import AlignmentVisualizer

viz = AlignmentVisualizer()

# Alignment overview
fig1 = viz.plot_alignment_overview(alignment_result)

# Mismatch distribution
fig2 = viz.plot_mismatch_distribution(alignment_result)

# Identity sliding window
fig3 = viz.plot_identity_sliding_window(alignment_result, window_size=50)
```

## File Organization

### Recommended Directory Structure

```
your_project/
├── data/
│   ├── chromatograms/          # Input .ab1/.scf files
│   │   ├── sample001.ab1
│   │   ├── sample002.ab1
│   │   └── ...
│   └── references/             # Reference sequences
│       ├── human_mtdna.fasta
│       └── custom_ref.fasta
├── results/                    # Analysis outputs
│   ├── sample001/
│   │   ├── traces.png
│   │   ├── quality.png
│   │   ├── alignment.txt
│   │   └── sequence.fasta
│   └── ...
└── scripts/                    # Custom analysis scripts
    └── my_analysis.py
```

## Advanced Usage

### Custom Analysis Pipeline

```python
import os
from pathlib import Path
from texi_trace import ChromatogramReader, SequenceAligner, MtDNAReference
from texi_trace.visualization import ChromatogramVisualizer, AlignmentVisualizer
from texi_trace.chromatogram.quality import QualityAssessment

def analyze_chromatogram_batch(input_dir, output_dir, reference_file=None):
    """Custom batch analysis function."""
    
    # Setup
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Load reference if provided
    if reference_file:
        mtdna_ref = MtDNAReference.load_from_file(reference_file)
        aligner = SequenceAligner(mtdna_ref.reference)
    else:
        mtdna_ref = MtDNAReference()
        aligner = SequenceAligner(mtdna_ref.reference)
    
    # Find chromatogram files
    files = list(input_path.glob('*.ab1')) + list(input_path.glob('*.scf'))
    
    results = []
    for file_path in files:
        print(f"Processing {file_path.name}...")
        
        try:
            # Read chromatogram
            reader = ChromatogramReader(str(file_path))
            data = reader.read()
            
            # Quality assessment
            quality_stats = QualityAssessment.assess_overall_quality(data['quality_scores'])
            
            # Sequence alignment
            alignment = aligner.align(data['sequence'])
            
            # Find variants
            variants = mtdna_ref.find_variants(data['sequence'])
            
            # Create output directory for this sample
            sample_output = output_path / file_path.stem
            sample_output.mkdir(exist_ok=True)
            
            # Generate visualizations
            viz = ChromatogramVisualizer()
            fig = viz.plot_sequence_overview(
                data['traces'], data['sequence'], data['quality_scores']
            )
            fig.savefig(sample_output / 'overview.png', dpi=300)
            
            # Save results
            results.append({
                'file': file_path.name,
                'sequence_length': len(data['sequence']),
                'mean_quality': quality_stats['mean_quality'],
                'alignment_identity': alignment['identity'],
                'num_variants': len(variants)
            })
            
        except Exception as e:
            print(f"Error processing {file_path.name}: {e}")
    
    return results

# Usage
results = analyze_chromatogram_batch(
    'data/chromatograms/',
    'results/',
    'data/references/human_mtdna.fasta'
)

# Print summary
for result in results:
    print(f"{result['file']}: {result['sequence_length']} bp, "
          f"Q={result['mean_quality']:.1f}, "
          f"ID={result['alignment_identity']:.1f}%, "
          f"{result['num_variants']} variants")
```

## Troubleshooting

### Common Issues

1. **File format not recognized**
   ```bash
   # Specify format explicitly
   texi-trace analyze sample.ab1 --format abi
   ```

2. **Low memory for large files**
   ```python
   # Process in chunks or reduce visualization resolution
   fig = viz.plot_traces(traces, sequence, start_pos=0, end_pos=1000)
   ```

3. **Missing dependencies**
   ```bash
   # Install specific dependency
   pip install biopython matplotlib
   ```

For more troubleshooting tips, see INSTALL.md.