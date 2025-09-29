# texi-trace

Built to scrutinise a Sanger chromatogram for mtDNA, check alignment to a small genome, and visualise the results.

## Overview

texi-trace is a comprehensive tool for analyzing Sanger chromatogram files, particularly optimized for mitochondrial DNA (mtDNA) sequencing analysis. It provides functionality for:

- **Chromatogram Analysis**: Read and analyze ABI (.ab1, .abi) and SCF (.scf) chromatogram files
- **Sequence Alignment**: Align sequences to reference genomes with optimized parameters for mtDNA
- **Quality Assessment**: Evaluate sequence quality and identify problematic regions
- **Visualization**: Generate publication-ready plots of chromatograms, alignments, and quality metrics
- **Batch Processing**: Analyze multiple files efficiently

## Installation

### Prerequisites
- Python 3.8 or higher
- pip package manager

### Install from source
```bash
git clone https://github.com/zorgster/texi-trace.git
cd texi-trace
pip install -r requirements.txt
pip install -e .
```

### Dependencies
The tool relies on several key bioinformatics and visualization libraries:
- BioPython for sequence analysis
- NumPy and SciPy for numerical computations
- Matplotlib and Seaborn for visualization
- Click for command-line interface

## Quick Start

### Command Line Interface

#### Analyze a single chromatogram file:
```bash
texi-trace analyze sample.ab1 --output-dir results/
```

#### Align to a reference genome:
```bash
texi-trace align sample.ab1 reference.fasta --output-dir alignment_results/
```

#### Batch process multiple files:
```bash
texi-trace batch data/chromatograms/ --reference-file mtdna_ref.fasta --output-dir batch_results/
```

### Python API

```python
from texi_trace import ChromatogramReader, SequenceAligner, ChromatogramVisualizer

# Read chromatogram data
reader = ChromatogramReader('sample.ab1')
data = reader.read()

# Analyze sequence quality
print(f"Sequence length: {len(data['sequence'])}")
print(f"Average quality: {sum(data['quality_scores'])/len(data['quality_scores']):.2f}")

# Create visualizations
visualizer = ChromatogramVisualizer()
fig = visualizer.plot_traces(data['traces'], data['sequence'])
fig.savefig('chromatogram_traces.png', dpi=300)

# Perform alignment
with open('reference.fasta', 'r') as f:
    ref_seq = ''.join(line.strip() for line in f if not line.startswith('>'))

aligner = SequenceAligner(ref_seq)
alignment = aligner.align(data['sequence'])
print(f"Identity: {alignment['identity']:.2f}%")
```

## File Formats Supported

### Input Formats
- **ABI files** (.ab1, .abi): Standard Applied Biosystems chromatogram format
- **SCF files** (.scf): Standard Chromatogram Format
- **FASTA files** (.fasta, .fa): Reference genome sequences

### Output Formats
- **PNG/JPEG**: High-resolution plots and visualizations
- **TXT**: Alignment results and analysis reports
- **FASTA**: Extracted sequences

## Features

### Chromatogram Analysis
- Parse trace data for all four bases (A, T, G, C)
- Extract called sequences and quality scores
- Identify peak positions and signal strengths
- Calculate signal-to-noise ratios

### Quality Assessment
- Phred quality score analysis
- Identification of low-quality regions
- Quality trimming recommendations
- Problem region detection

### Sequence Alignment
- Optimized alignment parameters for mtDNA
- Both global and local alignment options
- Mismatch and gap analysis
- Variant identification

### Visualization
- Chromatogram trace plots with base calls
- Quality score distributions
- Alignment overview plots
- Mismatch distribution analysis
- Identity sliding window plots

### mtDNA-Specific Features
- Human mtDNA reference sequences (rCRS)
- Gene region identification
- Mitochondrial genetic code translation
- Variant annotation

## Directory Structure

```
texi-trace/
├── src/texi_trace/          # Main package
│   ├── chromatogram/        # Chromatogram parsing
│   ├── alignment/           # Sequence alignment
│   ├── visualization/       # Plotting tools
│   └── cli.py              # Command-line interface
├── tests/                   # Unit tests
├── data/                    # Data storage
│   ├── chromatograms/      # Input chromatogram files
│   ├── references/         # Reference sequences
│   └── examples/           # Example data
├── examples/               # Usage examples
└── docs/                   # Documentation
```

## Usage Examples

### Basic Analysis Workflow

1. **Place your data**: Copy chromatogram files to `data/chromatograms/`
2. **Add references**: Place reference sequences in `data/references/`
3. **Run analysis**: Use CLI commands or Python API
4. **Review results**: Check output directory for plots and reports

### Example: Complete mtDNA Analysis

```python
import os
from texi_trace import ChromatogramReader, SequenceAligner, MtDNAReference
from texi_trace.visualization import ChromatogramVisualizer, AlignmentVisualizer

# Read chromatogram
reader = ChromatogramReader('data/chromatograms/sample.ab1')
data = reader.read()

# Use mtDNA reference
mtdna_ref = MtDNAReference()
aligner = SequenceAligner(mtdna_ref.reference)
alignment = aligner.align(data['sequence'])

# Create comprehensive visualizations
chrom_viz = ChromatogramVisualizer()
align_viz = AlignmentVisualizer()

# Generate plots
os.makedirs('output', exist_ok=True)

# Chromatogram overview
fig1 = chrom_viz.plot_sequence_overview(
    data['traces'], data['sequence'], data['quality_scores']
)
fig1.savefig('output/chromatogram_overview.png', dpi=300)

# Alignment results
fig2 = align_viz.plot_alignment_overview(alignment)
fig2.savefig('output/alignment_overview.png', dpi=300)

# Find variants
variants = mtdna_ref.find_variants(data['sequence'])
print(f"Found {len(variants)} variants")
```

## Contributing

Contributions are welcome! Please see our contributing guidelines for details on:
- Code style and standards
- Testing requirements
- Documentation updates
- Feature requests and bug reports

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Citation

If you use texi-trace in your research, please cite:

```
texi-trace: A tool for Sanger chromatogram analysis and mtDNA sequence alignment
Available at: https://github.com/zorgster/texi-trace
```

## Support

For questions, bug reports, or feature requests:
- Open an issue on GitHub
- Check the documentation in the `docs/` directory
- Review the examples in `examples/`

## Roadmap

- [ ] Support for additional file formats
- [ ] Integration with sequence databases
- [ ] Advanced variant calling
- [ ] Web-based interface
- [ ] Docker containerization
