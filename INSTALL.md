# Installation Guide for texi-trace

## System Requirements

- Python 3.8 or higher
- pip package manager
- Git (for development)

## Quick Installation

### Option 1: Install from Source (Recommended)

```bash
# Clone the repository
git clone https://github.com/zorgster/texi-trace.git
cd texi-trace

# Install dependencies
pip install -r requirements.txt

# Install the package in development mode
pip install -e .
```

### Option 2: Install Dependencies Only

If you want to use the package without installing it:

```bash
cd texi-trace
pip install -r requirements.txt
```

## Dependencies

### Core Dependencies (Required)
- **biopython>=1.80**: For sequence I/O and chromatogram parsing
- **numpy>=1.24.0**: Numerical computations
- **click>=8.0.0**: Command-line interface

### Analysis Dependencies
- **pairwise2>=1.0.0**: Sequence alignment (part of BioPython)
- **scipy>=1.10.0**: Scientific computing
- **pandas>=1.5.0**: Data manipulation

### Visualization Dependencies
- **matplotlib>=3.6.0**: Plotting and visualization
- **seaborn>=0.12.0**: Statistical visualization

### Development Dependencies
- **pytest>=7.0.0**: Testing framework
- **pytest-cov>=4.0.0**: Coverage testing
- **black>=22.0.0**: Code formatting
- **flake8>=5.0.0**: Code linting

## Verification

After installation, verify that everything works:

```bash
# Test basic functionality
python examples/run_example.py

# Test CLI
texi-trace --help

# Run tests (if pytest is installed)
pytest tests/
```

## Platform-Specific Notes

### Linux/macOS
```bash
# May need to install system dependencies for some packages
# Ubuntu/Debian:
sudo apt-get update
sudo apt-get install python3-dev build-essential

# macOS (with Homebrew):
brew install python
```

### Windows
```cmd
# Use Windows Subsystem for Linux (WSL) for best compatibility
# Or install Python from python.org and use Command Prompt/PowerShell
pip install -r requirements.txt
```

## Troubleshooting

### Common Issues

1. **BioPython installation fails**
   ```bash
   # Try installing with conda instead of pip
   conda install -c conda-forge biopython
   ```

2. **Matplotlib backend issues**
   ```bash
   # For headless systems, set backend
   export MPLBACKEND=Agg
   ```

3. **Permission errors**
   ```bash
   # Use user installation
   pip install --user -r requirements.txt
   ```

4. **Missing system libraries**
   ```bash
   # Ubuntu/Debian
   sudo apt-get install python3-dev libhdf5-dev

   # CentOS/RHEL
   sudo yum install python3-devel hdf5-devel
   ```

### Minimal Installation

For basic functionality without full dependencies:

```bash
# Install only essential packages
pip install numpy matplotlib click

# Test with reduced functionality
python -c "
from src.texi_trace.utils.sequence_utils import reverse_complement
print('Basic functionality works:', reverse_complement('ATCG'))
"
```

## Development Installation

For development work:

```bash
# Clone with development branch
git clone -b develop https://github.com/zorgster/texi-trace.git
cd texi-trace

# Install development dependencies
pip install -r requirements.txt
pip install -e .

# Install pre-commit hooks (optional)
pre-commit install
```

## Docker Installation (Future)

A Docker container will be available in future releases:

```bash
# Pull and run (not yet available)
docker pull zorgster/texi-trace:latest
docker run -v /path/to/data:/data zorgster/texi-trace analyze /data/sample.ab1
```

## Next Steps

After successful installation:

1. **Read the documentation**: Check README.md and examples/
2. **Prepare your data**: Place chromatogram files in `data/chromatograms/`
3. **Run analysis**: `texi-trace analyze your_file.ab1`
4. **Explore examples**: Run scripts in `examples/` directory

## Getting Help

If you encounter issues:

1. Check this installation guide
2. Review the troubleshooting section
3. Run the verification steps
4. Open an issue on GitHub with error details
5. Include your Python version and OS information