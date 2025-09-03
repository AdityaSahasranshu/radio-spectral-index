# -*- coding: utf-8 -*-
"""
Created on Wed Sep  3 16:22:14 2025

@author: adity
"""
"""
Example script to demonstrate spectral index analysis using sample data.

This script uses small sample cutouts from VLASS and LoTSS surveys to showcase
the spectral index calculation workflow. Perfect for testing and demonstration.

Usage:
    python run_example.py

Output:
    - output/spectral_index_sample.fits: Calculated spectral index map
    - Console output with analysis statistics and interpretation
    
Author: [Your Name]
Date: 2025
"""

import os
import sys
import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from radio_beam import Beam
from astropy.convolution import convolve
from reproject import reproject_interp
import warnings

# Suppress minor warnings for cleaner output
warnings.filterwarnings('ignore', category=RuntimeWarning)

def get_pixel_scale(header, wcs):
    """
    Extract pixel scale from FITS header, handling different WCS formats
    
    Parameters:
    -----------
    header : astropy.io.fits.Header
        FITS header containing WCS information
    wcs : astropy.wcs.WCS
        World Coordinate System object
        
    Returns:
    --------
    astropy.units.Quantity
        Pixel scale in degrees
    """
    try:
        # Try CDELT2 first (traditional format)
        if 'CDELT2' in header:
            return abs(header['CDELT2']) * u.deg
        # Try CD matrix format
        elif 'CD2_2' in header:
            return abs(header['CD2_2']) * u.deg
        # Try PC matrix with CDELT
        elif 'PC2_2' in header and 'CDELT2' in header:
            return abs(header['PC2_2'] * header['CDELT2']) * u.deg
        # Use WCS pixel scale calculation
        else:
            pixel_scale = wcs.proj_plane_pixel_scales()[1]  # Get Y-axis pixel scale
            return pixel_scale.to(u.deg)
    except Exception as e:
        print(f"Warning: Could not determine pixel scale from header: {e}")
        # Fallback: calculate from WCS
        pixel_scale = wcs.proj_plane_pixel_scales()[1]
        return pixel_scale.to(u.deg)

def convolve_if_needed(data, original_beam, target_beam, pixel_scale):
    """
    Convolve data to target beam if necessary
    
    Parameters:
    -----------
    data : numpy.ndarray
        Input data array
    original_beam : radio_beam.Beam
        Original beam of the data
    target_beam : radio_beam.Beam
        Target beam for convolution
    pixel_scale : astropy.units.Quantity
        Pixel scale of the data
        
    Returns:
    --------
    numpy.ndarray
        Convolved data array
    """
    if target_beam > original_beam:
        print(f"  Convolving from {original_beam} to {target_beam}")
        try:
            kernel = target_beam.deconvolve(original_beam).as_kernel(pixel_scale)
            return convolve(data, kernel, preserve_nan=True)
        except Exception as e:
            print(f"  Warning: Convolution failed ({e}), using original data")
            return data
    else:
        print("  No convolution needed")
    return data

def check_data_files():
    """Check if required sample data files exist"""
    vlass_file = os.path.join("data", "sample_vlass_cutout.fits")
    lotss_file = os.path.join("data", "sample_lotss_cutout.fits")
    
    missing_files = []
    if not os.path.exists(vlass_file):
        missing_files.append(vlass_file)
    if not os.path.exists(lotss_file):
        missing_files.append(lotss_file)
    
    if missing_files:
        print("ERROR: Required sample data files not found:")
        for file in missing_files:
            print(f"  - {file}")
        print("\nPlease ensure sample FITS files are in the 'data/' folder.")
        print("Expected files:")
        print("  - data/sample_vlass_cutout.fits")
        print("  - data/sample_lotss_cutout.fits")
        return False
    
    return True

def main():
    """
    Main function to calculate spectral index using sample data
    """
    
    print("=" * 70)
    print(" RADIO SPECTRAL INDEX ANALYSIS - SAMPLE DEMONSTRATION")
    print("=" * 70)
    print("This demo showcases spectral index calculation between VLASS (3 GHz)")
    print("and LoTSS (144 MHz) using real survey data cutouts.\n")
    
    # Check if sample data exists
    if not check_data_files():
        return 1
    
    # Define paths to sample data
    vlass_file = os.path.join("data", "sample_vlass_cutout.fits")
    lotss_file = os.path.join("data", "sample_lotss_cutout.fits")
    
    print(f"Loading sample FITS files...")
    print(f"  VLASS: {vlass_file}")
    print(f"  LoTSS: {lotss_file}")
    
    # Load FITS files
    try:
        print("  Reading VLASS data...")
        vlass_data = fits.getdata(vlass_file)
        vlass_header = fits.getheader(vlass_file)
        vlass_wcs = WCS(vlass_header, naxis=2)
        
        print("  Reading LoTSS data...")
        lotss_data = fits.getdata(lotss_file)
        lotss_header = fits.getheader(lotss_file)
        lotss_wcs = WCS(lotss_header, naxis=2)
        
    except Exception as e:
        print(f"ERROR loading FITS files: {e}")
        return 1
    
    print(f"\nOriginal data shapes:")
    print(f"  VLASS: {vlass_data.shape}")
    print(f"  LoTSS: {lotss_data.shape}")
    
    # Handle extra dimensions - reduce to 2D
    original_vlass_shape = vlass_data.shape
    original_lotss_shape = lotss_data.shape
    
    if vlass_data.ndim > 2:
        if vlass_data.ndim == 4:
            vlass_data = vlass_data[0, 0]  # Remove frequency and Stokes axes
        else:
            vlass_data = vlass_data[0]  # Remove one axis
        print(f"  Reduced VLASS from {original_vlass_shape} to {vlass_data.shape}")
    
    if lotss_data.ndim > 2:
        if lotss_data.ndim == 4:
            lotss_data = lotss_data[0, 0]  # Remove frequency and Stokes axes
        else:
            lotss_data = lotss_data[0]  # Remove one axis
        print(f"  Reduced LoTSS from {original_lotss_shape} to {lotss_data.shape}")
    
    # Data quality checks
    print(f"\nData quality assessment:")
    vlass_valid = np.isfinite(vlass_data).sum()
    lotss_valid = np.isfinite(lotss_data).sum()
    print(f"  VLASS valid pixels: {vlass_valid:,} ({100*vlass_valid/vlass_data.size:.1f}%)")
    print(f"  LoTSS valid pixels: {lotss_valid:,} ({100*lotss_valid/lotss_data.size:.1f}%)")
    
    # Reproject LoTSS to VLASS grid
    print(f"\nCoordinate reprojection:")
    print(f"  Reprojecting LoTSS to VLASS coordinate grid...")
    try:
        lotss_reprojected, footprint = reproject_interp((lotss_data, lotss_wcs), 
                                                       vlass_wcs, vlass_data.shape)
        overlap_fraction = np.sum(footprint > 0) / footprint.size
        print(f"  Reprojection successful - {100*overlap_fraction:.1f}% overlap")
    except Exception as e:
        print(f"ERROR in reprojection: {e}")
        return 1
    
    # Set beam information (typical values for these surveys)
    vlass_beam = Beam(major=2.5*u.arcsec, minor=2.5*u.arcsec, pa=0*u.deg)
    lotss_beam = Beam(major=6*u.arcsec, minor=6*u.arcsec, pa=0*u.deg)
    
    # Determine target beam (larger of the two)
    target_beam = max(vlass_beam, lotss_beam)
    
    print(f"\nBeam matching:")
    print(f"  VLASS beam: {vlass_beam}")
    print(f"  LoTSS beam: {lotss_beam}")
    print(f"  Target beam: {target_beam}")
    
    # Get pixel scale
    vlass_pixel_scale = get_pixel_scale(vlass_header, vlass_wcs)
    print(f"  Pixel scale: {vlass_pixel_scale.to(u.arcsec):.2f}")
    
    print(f"\nConvolving to common resolution...")
    vlass_conv = convolve_if_needed(vlass_data, vlass_beam, target_beam, vlass_pixel_scale)
    lotss_conv = convolve_if_needed(lotss_reprojected, lotss_beam, target_beam, vlass_pixel_scale)
    
    # Calculate noise thresholds (3-sigma detection limits)
    print(f"\nNoise analysis and thresholding:")
    
    # Calculate RMS in regions away from bright sources
    vlass_rms = np.nanstd(vlass_conv[np.abs(vlass_conv) < np.nanpercentile(np.abs(vlass_conv), 90)])
    lotss_rms = np.nanstd(lotss_conv[np.abs(lotss_conv) < np.nanpercentile(np.abs(lotss_conv), 90)])
    
    vlass_threshold = 3 * vlass_rms
    lotss_threshold = 3 * lotss_rms
    
    print(f"  VLASS RMS: {vlass_rms:.2e} Jy/beam, 3σ threshold: {vlass_threshold:.2e} Jy/beam")
    print(f"  LoTSS RMS: {lotss_rms:.2e} Jy/beam, 3σ threshold: {lotss_threshold:.2e} Jy/beam")
    
    # Create masks for valid detections
    vlass_mask = (vlass_conv > vlass_threshold) & np.isfinite(vlass_conv)
    lotss_mask = (lotss_conv > lotss_threshold) & np.isfinite(lotss_conv)
    combined_mask = vlass_mask & lotss_mask
    
    n_vlass_det = np.sum(vlass_mask)
    n_lotss_det = np.sum(lotss_mask)
    n_combined_det = np.sum(combined_mask)
    
    print(f"  VLASS detections: {n_vlass_det:,} pixels")
    print(f"  LoTSS detections: {n_lotss_det:,} pixels")
    print(f"  Combined detections: {n_combined_det:,} pixels")
    
    if n_combined_det == 0:
        print("\nERROR: No overlapping detections found!")
        print("This could be due to:")
        print("  - Very different noise levels between surveys")
        print("  - No real sources in this field")
        print("  - Coordinate alignment issues")
        return 1
    
    # Calculate spectral index
    vlass_freq = 3000  # MHz
    lotss_freq = 144   # MHz
    freq_ratio = vlass_freq / lotss_freq
    
    print(f"\nSpectral index calculation:")
    print(f"  VLASS frequency: {vlass_freq} MHz")
    print(f"  LoTSS frequency: {lotss_freq} MHz")
    print(f"  Frequency ratio: {freq_ratio:.1f}")
    
    # Initialize spectral index array
    spectral_index = np.full_like(vlass_conv, np.nan, dtype=np.float32)
    
    # Calculate only for pixels with valid detections in both surveys
    if n_combined_det > 0:
        flux_ratio = vlass_conv[combined_mask] / lotss_conv[combined_mask]
        spectral_index[combined_mask] = np.log(flux_ratio) / np.log(freq_ratio)
    
    # Remove any remaining problematic values
    spectral_index[~np.isfinite(spectral_index)] = np.nan
    
    # Final statistics
    valid_alpha = spectral_index[np.isfinite(spectral_index)]
    
    print(f"\n" + "=" * 50)
    print(" SPECTRAL INDEX RESULTS")
    print("=" * 50)
    
    if len(valid_alpha) > 0:
        # Basic statistics
        mean_alpha = np.mean(valid_alpha)
        median_alpha = np.median(valid_alpha)
        std_alpha = np.std(valid_alpha)
        min_alpha = np.min(valid_alpha)
        max_alpha = np.max(valid_alpha)
        
        print(f"Valid spectral index pixels: {len(valid_alpha):,}")
        print(f"Image coverage: {100*len(valid_alpha)/spectral_index.size:.1f}%")
        print(f"\nStatistical Summary:")
        print(f"  Mean:      {mean_alpha:+.3f}")
        print(f"  Median:    {median_alpha:+.3f}")
        print(f"  Std Dev:   {std_alpha:.3f}")
        print(f"  Range:     {min_alpha:+.3f} to {max_alpha:+.3f}")
        
        # Physical classification
        steep_spectrum = np.sum(valid_alpha < -0.7)
        typical_spectrum = np.sum((valid_alpha >= -0.7) & (valid_alpha <= -0.5))
        flat_spectrum = np.sum(valid_alpha > -0.5)
        
        print(f"\nSpectral Classification:")
        print(f"  Steep spectrum (α < -0.7):    {steep_spectrum:5,} pixels ({100*steep_spectrum/len(valid_alpha):5.1f}%)")
        print(f"  Typical synchrotron (-0.7≤α≤-0.5): {typical_spectrum:5,} pixels ({100*typical_spectrum/len(valid_alpha):5.1f}%)")
        print(f"  Flat spectrum (α > -0.5):     {flat_spectrum:5,} pixels ({100*flat_spectrum/len(valid_alpha):5.1f}%)")
        
        # Physical interpretation
        print(f"\nPhysical Interpretation:")
        if mean_alpha < -0.8:
            print("→ Dominated by steep-spectrum sources (aged synchrotron emission)")
        elif mean_alpha < -0.6:
            print("→ Typical synchrotron-dominated region (star formation, SNRs)")
        elif mean_alpha < -0.4:
            print("→ Mix of synchrotron and flat-spectrum sources")
        else:
            print("→ Unusual: dominated by flat-spectrum sources (AGN cores?)")
            
    else:
        print("WARNING: No valid spectral index values calculated!")
        print("Possible issues:")
        print("  - Noise thresholds too conservative")
        print("  - No overlapping bright sources")
        print("  - Data calibration problems")
        return 1
    
    # Create output directory
    os.makedirs("output", exist_ok=True)
    
    # Save results
    output_file = os.path.join("output", "spectral_index_sample.fits")
    print(f"\nSaving results:")
    print(f"  Output file: {output_file}")
    
    # Create comprehensive FITS header
    hdu = fits.PrimaryHDU(spectral_index, header=vlass_header)
    
    # Add extensive metadata
    hdu.header['OBJECT'] = 'Sample Spectral Index Map'
    hdu.header['BUNIT'] = 'dimensionless'
    hdu.header['COMMENT'] = 'Spectral index map: α = log(S₁/S₂)/log(ν₁/ν₂)'
    hdu.header['COMMENT'] = f'Between {vlass_freq} MHz (VLASS) and {lotss_freq} MHz (LoTSS)'
    hdu.header['COMMENT'] = 'Created from sample data for demonstration'
    
    # Survey information
    hdu.header['VFREQ'] = (vlass_freq, 'VLASS frequency [MHz]')
    hdu.header['LFREQ'] = (lotss_freq, 'LoTSS frequency [MHz]')
    hdu.header['VTHRESH'] = (vlass_threshold, 'VLASS 3-sigma threshold [Jy/beam]')
    hdu.header['LTHRESH'] = (lotss_threshold, 'LoTSS 3-sigma threshold [Jy/beam]')
    
    # Analysis results
    hdu.header['NPIXVAL'] = (len(valid_alpha), 'Number of valid spectral index pixels')
    hdu.header['COVERAGE'] = (100*len(valid_alpha)/spectral_index.size, 'Percentage coverage of image')
    
    if len(valid_alpha) > 0:
        hdu.header['ALPHAMN'] = (mean_alpha, 'Mean spectral index')
        hdu.header['ALPHAMD'] = (median_alpha, 'Median spectral index')
        hdu.header['ALPHASTD'] = (std_alpha, 'Spectral index standard deviation')
        hdu.header['ALPHAMIN'] = (min_alpha, 'Minimum spectral index')
        hdu.header['ALPHAMAX'] = (max_alpha, 'Maximum spectral index')
        hdu.header['NSTEEP'] = (steep_spectrum, 'Number of steep-spectrum pixels')
        hdu.header['NFLAT'] = (flat_spectrum, 'Number of flat-spectrum pixels')
    
    # Processing metadata
    hdu.header['SOFTWAR'] = 'radio-spectral-index v1.0'
    hdu.header['AUTHOR'] = '[Your Name]'
    
    try:
        hdu.writeto(output_file, overwrite=True)
        print(f"  ✓ FITS file saved successfully")
        
        # File size
        file_size = os.path.getsize(output_file) / 1024 / 1024  # MB
        print(f"  File size: {file_size:.1f} MB")
        
    except Exception as e:
        print(f"ERROR writing output file: {e}")
        return 1
    
    # Summary and recommendations
    print(f"\n" + "=" * 70)
    print(" ANALYSIS COMPLETE!")
    print("=" * 70)
    
    if len(valid_alpha) > 0:
        print(f"✓ Successfully calculated spectral indices for {len(valid_alpha):,} pixels")
        print(f"✓ Mean spectral index: {mean_alpha:+.3f} ± {std_alpha:.3f}")
        print(f"✓ Results saved to: {output_file}")
        
        print(f"\nNext Steps:")
        print(f"• View the spectral index map using DS9, CASA, or Python")
        print(f"• Typical α = -0.7 indicates synchrotron emission")
        print(f"• Compare with source catalogs for validation")
        print(f"• Consider error analysis for quantitative studies")
        
        # Quality assessment
        if len(valid_alpha) > 1000:
            quality = "Excellent"
        elif len(valid_alpha) > 100:
            quality = "Good"
        else:
            quality = "Limited"
            
        print(f"\nData Quality Assessment: {quality}")
        print(f"({len(valid_alpha):,} valid pixels for spectral analysis)")
        
    return 0

if __name__ == "__main__":
    sys.exit(main())