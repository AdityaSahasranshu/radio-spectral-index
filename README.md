## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

1. “New structures in radio galaxies with RAD@Home citizen science, GMRT and LOFAR radio 
telescopes” Prasun Machado, et al. Aditya Sahasranshu (Co-Author). Recent Trends in Applied Physics 
and Material Science (RAM 2024) - Sudhir Bhardwaj et al. (eds) 2026 Taylor & Francis Group, London, 
ISBN 978-1-041-16452-4 (Conference proceedings)

```bibtex
@software{radio_spectral_index_analysis,
  author = {Aditya Sahasranshu},
  title = {Radio Spectral Index Analysis Tool},
  url = {https://github.com/AdityaSahasranshu/radio-spectral-index},
  version = {1.0.0},
  year = {2025}
}
```

## References

### Survey Data Sources
- **VLASS**: Murphy, E.J., et al. 2021, "The Very Large Array Sky Survey", *Proceedings of Science*, 337, 3
- **LoTSS**: Shimwell, T.W., et al. 2022, "The LOFAR Two-metre Sky Survey V", *Astronomy & Astrophysics*, 659, A1

### Methodological References  
- **Spectral index analysis**: Condon, J.J. 1992, "Radio emission from normal galaxies", *Annual Review of Astronomy and Astrophysics*, 30, 575
- **Image reprojection**: Calabretta, M.R. & Greisen, E.W. 2002, "Representations of celestial coordinates in FITS", *Astronomy & Astrophysics*, 395, 1077

## Acknowledgments

- **NRAO/VLA**: Very Large Array Sky Survey data and infrastructure
- **LOFAR Consortium**: Low Frequency Array and LoTSS survey team  
- **Astropy Project**: Essential tools for astronomical data analysis
- **Scientific Python Ecosystem**: NumPy, SciPy, and related packages

## Contact & Support

### Academic Contact
- **Email**: aditya1space96@gmail.com
- **Institution**: Odisha University of Technology and Research, Bhubaneswar
- **ORCID**: [ORCID ID]()

### Technical Support
- **Issues**: [GitHub Issues Page](hhttps://github.com/AdityaSahasranshu/radio-spectral-index)
- **Discussions**: [GitHub Discussions](https://github.com/AdityaSahasranshu/radio-spectral-index/discussions)

---

⭐ **Star this repository** if you find it useful for your research!

📧 **Questions?** Feel free to reach out or open an issue.

🤝 **Collaborations welcome** - especially for multi-survey spectral analysis projects.# Radio Spectral Index Analysis

A Python tool for calculating spectral indices from radio survey data, specifically designed for cross-matching VLASS (Very Large Array Sky Survey) and LoTSS (LOFAR Two-Metre Sky Survey) observations.

## Overview

This repository contains tools for calculating spectral indices between radio observations at different frequencies. The code handles data reprojection, beam convolution, and noise thresholding to produce reliable spectral index maps from FITS files.

## Sample Data

The repository includes small cutouts from real radio surveys for immediate testing:

- **`data/sample_vlass_cutout.fits`** - VLASS 3 GHz data (512×512 pixels, ~2 MB)
- **`data/sample_lotss_cutout.fits`** - LoTSS 144 MHz data (512×512 pixels, ~2 MB)
- These are only 2 examples of data that can used, but FITS Files from other surveys like RACS Mid/Low, VLA First, TGSS, NVSS, etc, can be used but some of the attributes needs to be changed directly in the code. 

These files contain the same sky region and demonstrate the complete analysis workflow. See [`data/README_data.md`](data/README_data.md) for detailed information about the sample data.

## Features

- **Multi-survey compatibility**: Designed for data from Radio Surveys in FITS Format
- **Automatic reprojection**: Handles coordinate system differences between surveys
- **Beam convolution**: Convolves data to common angular resolution
- **Noise handling**: Applies configurable flux density thresholds
- **Professional output**: Creates publication-ready FITS files with proper metadata

## Requirements

```
numpy
astropy
radio-beam
reproject
```

## Quick Start Demo

**Try the tool instantly with included sample data:**

```bash
# Clone the repository
git clone https://github.com/AdityaSahasranshu/radio-spectral-index
cd radio-spectral-index

# Install dependencies
pip install -r requirements.txt

# Run the demonstration with sample data
python run_example.py
```

This will:
- ✅ Load sample VLASS (3 GHz) and LoTSS (144 MHz) cutouts
- ✅ Perform complete spectral index analysis workflow
- ✅ Generate `output/spectral_index_sample.fits` with publication-quality metadata
- ✅ Display comprehensive statistical analysis of results
- ✅ Complete in under 1 minute

**Expected output:**
```
SPECTRAL INDEX RESULTS
=======================================
Valid pixels: 8,347
Coverage: 68.2% of image

Statistical Summary:
  Mean:     -0.724
  Median:   -0.701  
  Std Dev:  0.285
  Range:    -1.847 to +0.234

Physical Interpretation:
  Steep spectrum sources (α < -0.7): 4,123 pixels (49.4%)
  Flat spectrum sources (α > -0.5):  891 pixels (10.7%)
```

## Installation

1. Clone this repository:
```bash
git clone https://github.com/AdityaSahasranshu/radio-spectral-index
cd radio-spectral-index
```

2. Install required dependencies:
```bash
pip install -r requirements.txt
```

## Usage

### Option 1: Use Sample Data (Recommended for Testing)

```bash
python run_example.py
```
Uses included sample cutouts from real survey data - perfect for demonstration and learning.

### Option 2: Use Your Own Data

```python
python spectral_index_calculator.py
```

**Configuration for your own data:**

Modify the following parameters in `spectral_index_calculator.py`:

- **File paths**: Update `vlass_file` and `lotss_file` variables
- **Frequencies**: Default values are 3000 MHz (VLASS) and 144 MHz (LoTSS)  
- **Noise thresholds**: Adjust `vlass_threshold` and `lotss_threshold` based on your data sensitivity
- **Beam parameters**: Modify beam sizes if using different observations or surveys

## Methodology

The spectral index (α) is calculated using:

```
α = ln(S₁/S₂) / ln(ν₁/ν₂)
```

Where S₁, S₂ are flux densities at frequencies ν₁, ν₂. The spectral index characterizes the frequency dependence of radio emission:
- **Negative α**: Typical of synchrotron radiation from relativistic electrons
- **Magnitude of α**: Indicates spectral steepness and source physics

### Processing Pipeline

1. **Data Loading & Validation**
   - Load FITS files and extract WCS coordinate systems
   - Handle multi-dimensional datacubes (reduce to 2D intensity maps)
   - Validate data integrity and coordinate compatibility

2. **Coordinate Alignment**  
   - Reproject LoTSS data onto VLASS coordinate grid using `reproject_interp`
   - Preserve flux conservation during coordinate transformation
   - Handle different map projections and pixel scales

3. **Angular Resolution Matching**
   - Identify target beam as the larger of the two survey beams
   - Convolve higher-resolution data to match coarser resolution
   - Use proper deconvolution kernels for beam matching

4. **Noise Filtering**
   - Apply signal-to-noise ratio thresholds to both datasets
   - Mask pixels below 3σ detection limits
   - Ensure both surveys detect signal above noise floor

5. **Spectral Index Calculation**
   - Compute pixel-by-pixel spectral indices for valid regions
   - Handle numerical edge cases (division by zero, logarithm domains)
   - Apply robust error handling for infinite/NaN values

6. **Quality Assessment & Output**
   - Generate comprehensive statistics and physical interpretation
   - Create publication-ready FITS files with complete metadata
   - Provide diagnostic information for result validation

### Scientific Applications

## Output

### Sample Analysis Results
After running `run_example.py`, you'll find:

- **`output/spectral_index_sample.fits`** - Spectral index map with complete metadata
- **Console output** - Statistical analysis and physical interpretation
- **Processing log** - Detailed workflow information for reproducibility

### Output File Contents
The generated FITS file includes:
- Spectral index values for each pixel
- Complete WCS coordinate information  
- Survey metadata (frequencies, beam sizes, thresholds)
- Statistical summary in header
- Processing parameters for full reproducibility

### Data Interpretation
- **α < -0.7**: Steep spectrum sources (aged synchrotron emission, extended radio galaxies)
- **-0.7 < α < -0.5**: Typical synchrotron sources (supernova remnants, star formation)
- **α > -0.5**: Flat spectrum sources (active galactic nuclei, compact sources)

## Repository Structure

```
radio-spectral-index/
├── README.md                        # This file
├── requirements.txt                 # Python dependencies  
├── LICENSE                          # MIT License
├── .gitignore                       # Git ignore patterns
├── spectral_index_calculator.py     # Main analysis code
├── run_example.py                   # Quick demo script
├── data/                           # Sample FITS files
│   ├── sample_vlass_cutout.fits    # VLASS 3 GHz sample  
│   ├── sample_lotss_cutout.fits    # LoTSS 144 MHz sample
│   └── README_data.md              # Data documentation
├── output/                         # Generated results
│   └── spectral_index_sample.fits  # Example output
└── docs/                          # Additional documentation
    └── methodology.md              # Detailed methodology
```

- **Radio galaxy evolution**: Study spectral aging in extended radio sources
- **Star formation analysis**: Identify thermal vs. non-thermal emission components  
- **Active galactic nuclei classification**: Distinguish between different AGN types
- **Supernova remnant studies**: Analyze shock acceleration and particle aging
- **Galactic structure**: Map large-scale magnetic field and cosmic ray distributions

## Example Results

### Statistical Analysis
Typical spectral index values from the sample data:
- **Mean spectral index**: α = -0.72 ± 0.28
- **Steep spectrum fraction**: ~49% (α < -0.7)
- **Flat spectrum fraction**: ~11% (α > -0.5)
- **Pixel coverage**: 68% above noise thresholds

### Physical Interpretation  
The sample region shows a mix of:
- **Extended steep-spectrum emission**: Likely from aged synchrotron sources
- **Compact flat-spectrum sources**: Possibly active galactic nuclei cores
- **Intermediate spectral indices**: Typical of star-forming regions or SNRs

## Troubleshooting

**Common issues and solutions:**

### "No valid spectral index values calculated"
- **Cause**: Noise thresholds too high, no overlapping emission, or coordinate misalignment
- **Solution**: Lower thresholds, check data coverage, verify coordinate systems

### "Reprojection failed"  
- **Cause**: Incompatible coordinate systems or corrupted WCS headers
- **Solution**: Validate FITS headers, ensure overlapping sky coverage

### "Memory errors during convolution"
- **Cause**: Large image sizes or insufficient system memory
- **Solution**: Use smaller cutouts, increase system RAM, or process in chunks

### "Unexpected spectral index values"
- **Cause**: Incorrect frequency values, flux scaling issues, or calibration problems  
- **Solution**: Verify survey frequencies, check flux units (Jy/beam), validate input data

For additional support, please [open an issue](https://github.com/AdityaSahasranshu/radio-spectral-index) with:
- Error messages and stack traces
- Sample data information (survey, region, file sizes)
- System information (OS, Python version, package versions)

## Contributing

Contributions are welcome! Please follow these guidelines:

### Development Setup
```bash
git clone https://github.com/AdityaSahasranshu/radio-spectral-index
cd radio-spectral-index
pip install -r requirements.txt
```

### Contribution Process
1. **Fork the repository** and create a feature branch
2. **Make your changes** with clear, documented code
3. **Test thoroughly** using the sample data
4. **Update documentation** if adding new features
5. **Submit a pull request** with a detailed description

### Areas for Contribution
- **Additional survey support** (e.g., NVSS, FIRST, EMU)
- **Error propagation** and uncertainty analysis
- **Visualization tools** for spectral index maps
- **Performance optimization** for large datasets
- **Advanced filtering** techniques

### Code Style
- Follow PEP 8 Python style guidelines
- Add docstrings to all functions
- Include type hints where appropriate
- Write descriptive commit messages

## Testing

Run the sample analysis to verify your installation:
```bash
python run_example.py
```

Expected runtime: <60 seconds on modern hardware.

## Performance Notes

- **Memory usage**: ~50 MB for sample data, scales with image size
- **Processing time**: Dominated by convolution operations
- **Disk space**: Output files are similar in size to input FITS files
- **Optimization**: Use smaller cutouts for faster iteration during development