# Sample Data Documentation

This folder contains small sample cutouts from real radio survey data for demonstration and testing purposes.

## Files

### `UNK.AUTH-P-VLASS3.fits`
- **Survey**: Very Large Array Sky Survey 3.1 (VLASS)  
- **Frequency**: 3 GHz (S-band)
- **Angular Resolution**: ~2.5 arcsec
- **Sensitivity**: ~120 μJy/beam RMS
- **3σ Value**: 0.00021 Jy/beam
- **Data Type**: Stokes I intensity map
- **Size**: 60"

### `astron.fits`
- **Survey**: LOFAR Two-metre Sky Survey (LoTSS)
- **Frequency**: 144 MHz (L-band)  
- **Angular Resolution**: ~6 arcsec
- **Sensitivity**: ~83 μJy/beam RMS
- **3σ Value**: 0.00024 Jy/beam
- **Data Type**: Stokes I intensity map
- **Size**: 60"

## Data Sources

- **VLASS**: National Radio Astronomy Observatory (NRAO)
  - Public data release available at: https://science.nrao.edu/vlass
  - Citation: [Include proper VLASS citation]

- **LoTSS**: LOFAR Surveys Key Science Project  
  - Public data release available at: https://lotss.org/
  - Citation: [Include proper LoTSS citation]

## Usage Notes

1. **Coordinate Systems**: Both cutouts cover the same sky region but may use different coordinate projections
2. **Flux Scaling**: Data are in Jy/beam units
3. **Noise Properties**: RMS noise varies across the images due to primary beam correction
4. **Missing Data**: NaN values indicate regions outside survey coverage

## Creating Your Own Cutouts

If you want to create custom cutouts from full survey images:
- **Easiest Approach**- I will recomend using software like Aladin and DS9 to create cutcouts from larger survey data. 

## Expected Results

When running the spectral index analysis on this sample data, you should expect:

- **Typical spectral indices**: -0.8 to -0.5 (steep synchrotron spectrum)
- **Coverage**: ~60-80% of pixels above noise threshold
- **Source types**: Mix of extended emission and compact sources
- **Processing time**: < 1 minute on standard hardware

## Troubleshooting

**Common issues and solutions:**

1. **"No valid spectral index values"**
   - Check that both FITS files are in the correct directory
   - Verify files aren't corrupted (try opening in DS9)
   - Noise thresholds might be too high

2. **Reprojection errors**
   - Ensure both files cover overlapping sky regions
   - Check coordinate system compatibility

3. **Memory issues**
   - Sample files should be small (<10 MB each)
   - If using larger cutouts, consider reducing image size

## File Formats

All FITS files follow the standard format:
- **Header**: Contains WCS coordinate information, survey metadata
- **Data**: 2D numpy arrays with flux values in Jy/beam
- **Units**: Brightness units specified in BUNIT header keyword