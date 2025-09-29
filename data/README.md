# Data Directory

This directory is intended for storing:

## Input Data
- Chromatogram files (.ab1, .abi, .scf)
- Reference genome sequences (.fasta, .fa)
- mtDNA reference sequences

## Example Data Structure
```
data/
├── chromatograms/
│   ├── sample001.ab1
│   ├── sample002.ab1
│   └── ...
├── references/
│   ├── human_mtdna_reference.fasta
│   ├── rCRS.fasta
│   └── ...
└── examples/
    ├── test_chromatogram.ab1
    └── test_reference.fasta
```

## File Formats Supported

### Chromatogram Files
- **ABI format** (.ab1, .abi): Standard Sanger sequencing output
- **SCF format** (.scf): Standard Chromatogram Format

### Reference Files  
- **FASTA format** (.fasta, .fa): Standard sequence format for reference genomes

## Notes
- Place your chromatogram files in the `chromatograms/` subdirectory
- Reference sequences should go in `references/`
- Test data can be placed in `examples/`
- The tool can process files directly from these locations using relative paths