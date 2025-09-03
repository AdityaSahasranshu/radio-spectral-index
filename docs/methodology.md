# Detailed Methodology

This document provides an in-depth explanation of the spectral index analysis methodology implemented in this repository.

## Theoretical Background

### Spectral Index Definition

The spectral index (α) characterizes how radio flux density varies with frequency:

```
S(ν) = S₀ × (ν/ν₀)^α
```

Where:
- S(ν) = flux density at frequency ν
- S₀ = reference flux density at frequency ν₀
- α = spectral index

For two frequencies, the spectral index is calculated as:

```
α = ln(S₁/S₂) / ln(ν₁/ν₂)
```

### Physical Interpretation

**Spectral index values indicate different emission mechanisms:**

- **α < -1.0**: Very steep spectrum, often indicates old synchrotron emission
- **-1.0 < α < -0.5**: Typical synchrotron emission from cosmic ray electrons
- **-0.5 < α < 0**: Flat spectrum, may indicate self-absorption or AGN cores
- **α > 0**: Rising spectrum, could indicate thermal emission or unusual physics

### Synchrotron Radiation Theory

Most radio sources emit via synchrotron radiation from relativistic electrons spiraling in magnetic fields. For a power-law electron energy distribution N(E) ∝ E^(-p), the spectral index is:

```
α = -(p-1)/2
```

Typical values: p ≈ 2.2-2.8, giving α ≈ -0.6 to -0.9.

## Data Processing Pipeline

### 1. Data Loading and Validation

```python
# Load FITS files with error handling
vlass_data = fits.getdata(vlass_file)
vlass_header = fits.getheader(vlass_file)
vlass_wcs = WCS(vlass_header, naxis=2)
```

**Key considerations:**
- Handle different FITS formats (2D, 3D, 4D datacubes)
- Extract proper World Coordinate System (WCS) information
- Validate data integrity and units

### 2. Coordinate System Alignment

Radio surveys use different coordinate projections and pixel grids. We use the `reproject` library to align data:

```python
from reproject import reproject_interp

lotss_reprojected, footprint = reproject_interp(
    (lotss_data, lotss_wcs), 
    vlass_wcs, 
    vlass_data.shape
)
```

**Reprojection methods:**
- **Interpolation** (`reproject_interp`): Fast, preserves flux for smooth sources
- **Drizzling** (`reproject_exact`): More accurate for complex sources, slower
- **Adaptive** (`reproject_adaptive`): Best of both, computational expensive

### 3. Angular Resolution Matching

Different surveys have different angular resolutions (beam sizes). We convolve to a common resolution:

```python
from radio_beam import Beam
from astropy.convolution import convolve

# Define beams
vlass_beam = Beam(major=2.5*u.arcsec, minor=2.5*u.arcsec, pa=0*u.deg)
lotss_beam = Beam(major=6*u.arcsec, minor=6*u.arcsec, pa=0*u.deg)

# Convolve to larger beam
target_beam = max(vlass_beam, lotss_beam)
kernel = target_beam.deconvolve(vlass_beam).as_kernel(pixel_scale)
vlass_convolved = convolve(vlass_data, kernel)
```

**Mathematical basis:**
For Gaussian beams, the convolution kernel is also Gaussian:
```
σ_kernel = √(σ_target² - σ_original²)
```

### 4. Noise Analysis and Thresholding

**RMS Noise Estimation:**
We estimate noise in regions away from bright sources:

```python
# Use 90th percentile to exclude bright sources
clean_data = data[np.abs(data) < np.nanpercentile(np.abs(data), 90)]
rms_noise = np.nanstd(clean_data)
threshold = 3 * rms_noise  # 3-sigma detection limit
```

**Detection Masks:**
Only pixels detected above the noise threshold in both surveys are used:

```python
vlass_mask = (vlass_data > vlass_threshold) & np.isfinite(vlass_data)
lotss_mask = (lotss_data > lotss_threshold) & np.isfinite(lotss_data)
combined_mask = vlass_mask & lotss_mask
```

### 5. Spectral Index Calculation

For pixels meeting detection criteria:

```python
flux_ratio = vlass_flux[mask] / lotss_flux[mask]
frequency_ratio = vlass_freq / lotss_freq
spectral_index[mask] = np.log(flux_ratio) / np.log(frequency_ratio)
```

**Error handling:**
- Remove infinite values (division by zero)
- Remove NaN values (invalid operations)
- Apply reasonable bounds (e.g., -3 < α < +2)

## Error Analysis

### Sources of Uncertainty

1. **Thermal Noise**: Random noise in receivers
2. **Calibration Errors**: Systematic flux scale uncertainties (~5-10%)
3. **Resolution Effects**: Beam size differences affect extended sources
4. **Coordinate Alignment**: Reprojection introduces interpolation errors
5. **Variability**: Sources may change between observations

### Error Propagation

For spectral index error δα:

```
δα = √[(δS₁/S₁)² + (δS₂/S₂)²] / ln(ν₁/ν₂)
```

Where δS₁, δS₂ are flux density uncertainties.

### Quality Metrics

- **Signal-to-noise ratio**: S/σ > 3 for reliable detection
- **Coverage fraction**: Percentage of image with valid spectral indices
- **Statistical consistency**: Compare with known source populations

## Survey-Specific Considerations

### VLASS (Very Large Array Sky Survey)

- **Frequency**: 2-4 GHz (S-band)
- **Resolution**: 2.5" typical
- **Sensitivity**: ~120 μJy/beam RMS
- **Polarization**: Full Stokes
- **Coverage**: Entire sky δ > -40°

### LoTSS (LOFAR Two-metre Sky Survey)

- **Frequency**: 120-168 MHz (centered at 144 MHz)
- **Resolution**: 6" typical
- **Sensitivity**: ~83 μJy/beam RMS
- **Polarization**: Stokes I
- **Coverage**: Northern sky δ > +27°

### Cross-matching Challenges

1. **Frequency difference**: Large lever arm (3000/144 ≈ 21) is good for spectral index precision
2. **Resolution difference**: Factor of ~2.4 difference requires careful convolution
3. **Sensitivity difference**: Similar depths but different systematics
4. **Temporal baseline**: Observations separated by years may show source evolution

## Validation and Quality Control

### Statistical Validation

1. **Source Count Consistency**: Compare with known radio source populations
2. **Spectral Index Distribution**: Should peak around α ≈ -0.7 for synchrotron sources
3. **Spatial Coherence**: Extended sources should show coherent spectral structure

### Visual Inspection

- **Overlay on optical images**: Identify counterparts and source types
- **Compare with catalogs**: Cross-match with known radio sources
- **Check for artifacts**: Look for systematic patterns indicating processing errors

### Known Issues and Limitations

1. **Extended Source Recovery**: Large-scale structure may be filtered differently
2. **Point Source Assumption**: Analysis assumes sources are unresolved
3. **Primary Beam Correction**: Edge effects from telescope sensitivity patterns
4. **Bandpass Variations**: Frequency-dependent instrumental effects

## Advanced Techniques

### Error Weighting

For improved precision, weight pixels by their signal-to-noise ratio:

```python
weights = 1 / (error_vlass**2 + error_lotss**2)
weighted_mean_alpha = np.average(spectral_indices, weights=weights)
```

### Resolved Source Analysis

For extended sources, consider:
- **Region-based analysis**: Average over source components
- **Radial profiles**: Study spectral steepening with distance from core
- **Spectral tomography**: Map spectral variations within sources

### Multi-frequency Extensions

This framework can be extended to additional frequencies:
- **Power-law fitting**: Fit α across multiple frequencies
- **Curvature detection**: Identify deviations from simple power laws
- **Physical modeling**: Fit synchrotron aging or absorption models

## Implementation Notes

### Performance Optimization

- **Memory management**: Process large images in chunks
- **Parallel processing**: Use multiprocessing for independent pixels
- **Caching**: Store intermediate results for iterative analysis

### Code Structure

```python
class SpectralIndexAnalyzer:
    def __init__(self, freq1, freq2):
        self.freq1 = freq1
        self.freq2 = freq2
    
    def load_data(self, file1, file2):
        # Load and validate FITS files
        
    def align_data(self):
        # Reproject to common grid
        
    def match_resolution(self):
        # Convolve to common beam
        
    def calculate_spectral_index(self):
        # Main calculation with error handling
        
    def save_results(self, output_file):
        # Write FITS with metadata
```

## References

### Key Papers

1. **Condon & Ransom (2016)**: "Essential Radio Astronomy" - Comprehensive theoretical background
2. **Shimwell et al. (2022)**: LoTSS DR2 - Survey description and data products
3. **Lacy et al. (2020)**: VLASS survey paper - Technical specifications
4. **Hardcastle et al. (2016)**: "Radio-loud AGN spectral index studies" - Interpretation guidelines

### Technical Documentation

- **Astropy**: https://docs.astropy.org/ - FITS handling and WCS
- **reproject**: https://reproject.readthedocs.io/ - Image reprojection
- **radio-beam**: https://radio-beam.readthedocs.io/ - Beam manipulation

### Data Access

- **VLASS**: https://science.nrao.edu/vlass - Official data releases
- **LoTSS**: https://lofar-surveys.org/ - Survey data and catalogs

## Future Developments

### Planned Enhancements

1. **Automated source finding**: Integrate with PyBDSF or similar
2. **Error propagation**: Full uncertainty analysis
3. **Catalog cross-matching**: Automated identification of known sources
4. **Visualization tools**: Interactive spectral index maps

### Research Applications

- **Galaxy evolution**: Spectral aging in radio galaxies
- **Star formation**: Thermal/non-thermal emission separation  
- **Cosmology**: High-redshift radio source studies
- **Galactic astronomy**: Pulsar and SNR characterization

---

*This methodology forms the scientific foundation for reliable spectral index analysis from multi-frequency radio survey data.*