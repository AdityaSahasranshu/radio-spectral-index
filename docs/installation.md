# Installation Guide

This guide provides detailed installation instructions for the radio spectral index analysis toolkit.

## System Requirements

### Operating System
- **Linux**: Ubuntu 18.04+, CentOS 7+, or similar
- **macOS**: 10.14+ (Mojave or later)
- **Windows**: Windows 10+ with Python support

### Python Version
- **Python 3.8** or higher (Python 3.9+ recommended)
- **64-bit installation** required for large FITS file handling

### Hardware Requirements
- **RAM**: Minimum 4 GB, 8 GB+ recommended for large datasets
- **Storage**: 500 MB for software, additional space for data files
- **CPU**: Modern multi-core processor recommended for convolution operations

## Installation Methods

### Method 1: Quick Installation (Recommended)

```bash
# Clone the repository
git clone https://github.com/yourusername/radio-spectral-index.git
cd radio-spectral-index

# Install dependencies
pip install -r requirements.txt

# Test installation
python run_example.py
```

### Method 2: Development Installation

```bash
# Clone repository
git clone https://github.com/yourusername/radio-spectral-index.git
cd radio-spectral-index

# Create virtual environment
python -m venv spectral_env
source spectral_env/bin/activate  # On Windows: spectral_env\Scripts\activate

# Install in development mode
pip install -e .
pip install -r requirements.txt

# Run tests
python tests/test_basic_functionality.py
```

### Method 3: Conda Installation

```bash
# Create conda environment
conda create -n spectral python=3.9
conda activate spectral

# Install astronomy packages via conda
conda install -c conda-forge astropy radio-beam reproject numpy

# Clone and install
git clone https://github.com/yourusername/radio-spectral-index.git
cd radio-spectral-index
```

## Dependency Installation

### Core Dependencies

The following packages are required:

```bash
pip install numpy>=1.21.0
pip install astropy>=5.0
pip install radio-beam>=0.3.3
pip install reproject>=0.8.0
```

### Optional Dependencies

For enhanced functionality:

```bash
# For plotting and visualization
pip install matplotlib>=3.5.0

# For advanced statistics
pip install scipy>=1.7.0

# For faster operations
pip install numba>=0.56.0

# For FITS file optimization
pip install fitsio>=1.1.0
```

### Platform-Specific Notes

#### Linux
```bash
# Ubuntu/Debian: Install system dependencies
sudo apt-get update
sudo apt-get install python3-dev gcc gfortran

# CentOS/RHEL: Install system dependencies  
sudo yum install python3-devel gcc gcc-gfortran
```

#### macOS
```bash
# Install Xcode command line tools
xcode-select --install

# Using Homebrew
brew install python
```

#### Windows
```bash
# Install Microsoft C++ Build Tools
# Download from: https://visualstudio.microsoft.com/visual-cpp-build-tools/

# Using conda is recommended on Windows
conda install -c conda-forge astropy radio-beam reproject
```

## Verification

### Quick Test
```bash
python -c "import numpy, astropy, radio_beam, reproject; print('All packages imported successfully!')"
```

### Full Test Suite
```bash
python tests/test_basic_functionality.py
```

Expected output:
```
BASIC FUNCTIONALITY TESTS
========================================
✓ Data Loading                   PASS
✓ Spectral Index Calculation     PASS
✓ Noise Handling                 PASS
✓ Beam Handling                  PASS
✓ Coordinate Handling            PASS
✓ Error Propagation              PASS
✓ File I/O Operations            PASS
✓ Integration Test               PASS
Tests passed: 8/8 (100%)
```

### Sample Data Test
```bash
python run_example.py
```

Expected output should include:
```
SPECTRAL INDEX RESULTS
=======================================
Valid pixels: [number]
Coverage: [percentage]% of image
Mean spectral index: [value] ± [error]
```

## Troubleshooting

### Common Issues

#### ImportError: No module named 'astropy'
```bash
# Solution: Install astropy
pip install astropy>=5.0
```

#### ImportError: No module named 'radio_beam'
```bash
# Solution: Install radio-beam
pip install radio-beam>=0.3.3
```

#### Memory errors with large files
```bash
# Solution: Increase virtual memory or use smaller cutouts
# Linux/macOS:
ulimit -v unlimited

# Or reduce image size in your analysis
```

#### Slow convolution operations
```bash
# Solution: Install additional optimization packages
pip install numba scipy

# Or consider using smaller beam matching
```

### Package Version Conflicts

Check package versions:
```bash
pip list | grep -E "(numpy|astropy|radio-beam|reproject)"
```

Create clean environment:
```bash
# Remove existing environment
pip uninstall numpy astropy radio-beam reproject

# Reinstall with specific versions
pip install -r requirements.txt
```

### Performance Optimization

#### For Large Datasets
```bash
# Install optimized linear algebra
pip install numpy[mkl]  # Intel MKL acceleration

# Use multiple cores
export OMP_NUM_THREADS=4
```

#### Memory Management
```python
# In your analysis scripts
import numpy as np
# Process data in chunks for memory efficiency
chunk_size = 1000
```

## Development Setup

### For Contributors

```bash
# Fork the repository on GitHub
git clone https://github.com/yourusername/radio-spectral-index.git
cd radio-spectral-index

# Create development environment
python -m venv dev_env
source dev_env/bin/activate

# Install in development mode
pip install -e .
pip install -r requirements.txt

# Install development tools
pip install pytest black flake8 mypy

# Run development tests
python -m pytest tests/
```

### Code Formatting
```bash
# Format code
black spectral_index_calculator.py run_example.py

# Check style
flake8 *.py

# Type checking
mypy spectral_index_calculator.py
```

## Docker Installation (Advanced)

Create Dockerfile:
```dockerfile
FROM python:3.9-slim

WORKDIR /app
COPY requirements.txt .
RUN pip install -r requirements.txt

COPY . .
CMD ["python", "run_example.py"]
```

Build and run:
```bash
docker build -t spectral-index .
docker run -v $(pwd)/data:/app/data spectral-index
```

## Updating

### Update from Git
```bash
git pull origin main
pip install -r requirements.txt --upgrade
```

### Update Dependencies
```bash
pip install --upgrade numpy astropy radio-beam reproject
```

## Getting Help

### Documentation
- **Main README**: [README.md](../README.md)
- **Methodology**: [methodology.md](methodology.md)
- **GitHub Issues**: Report bugs and request features

### Community Support
- **GitHub Discussions**: Ask questions and share experiences
- **Astropy Community**: For astronomy-specific Python questions
- **Stack Overflow**: Tag questions with `astropy` and `radio-astronomy`

### Professional Support
For institutional users requiring professional support, please contact the maintainers through the repository.

---

*For the most up-to-date installation instructions, always refer to the official repository README.*