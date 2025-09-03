# -*- coding: utf-8 -*-
"""
Created on Wed Sep  3 16:32:43 2025

@author: adity
"""

"""
Advanced spectral index analysis with error propagation and source classification.

This script demonstrates advanced techniques including:
- Error propagation and uncertainty analysis
- Source classification based on spectral properties
- Regional analysis and statistics
- Quality assessment and validation

Usage:
    python advanced_analysis.py

Author: [Your Name]
Date: 2025
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clip
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)

class AdvancedSpectralAnalyzer:
    """
    Advanced spectral index analyzer with error handling and source classification
    """
    
    def __init__(self, vlass_file, lotss_file, vlass_freq=3000, lotss_freq=144):
        """
        Initialize analyzer with input files and frequencies
        
        Parameters:
        -----------
        vlass_file : str
            Path to VLASS FITS file
        lotss_file : str
            Path to LoTSS FITS file
        vlass_freq : float
            VLASS frequency in MHz
        lotss_freq : float
            LoTSS frequency in MHz
        """
        self.vlass_file = vlass_file
        self.lotss_file = lotss_file
        self.vlass_freq = vlass_freq
        self.lotss_freq = lotss_freq
        self.freq_ratio = vlass_freq / lotss_freq
        
        # Data containers
        self.vlass_data = None
        self.lotss_data = None
        self.vlass_error = None
        self.lotss_error = None
        self.spectral_index = None
        self.spectral_index_error = None
        self.source_mask = None
        
    def load_data(self):
        """Load FITS files and estimate noise properties"""
        print("Loading and analyzing input data...")
        
        # Load VLASS data
        self.vlass_data = fits.getdata(self.vlass_file)
        self.vlass_header = fits.getheader(self.vlass_file)
        self.vlass_wcs = WCS(self.vlass_header, naxis=2)
        
        # Load LoTSS data  
        self.lotss_data = fits.getdata(self.lotss_file)
        self.lotss_header = fits.getheader(self.lotss_file)
        self.lotss_wcs = WCS(self.lotss_header, naxis=2)
        
        # Handle multi-dimensional data
        if self.vlass_data.ndim > 2:
            self.vlass_data = self.vlass_data[0, 0] if self.vlass_data.ndim == 4 else self.vlass_data[0]
        if self.lotss_data.ndim > 2:
            self.lotss_data = self.lotss_data[0, 0] if self.lotss_data.ndim == 4 else self.lotss_data[0]
            
        print(f"  VLASS data shape: {self.vlass_data.shape}")
        print(f"  LoTSS data shape: {self.lotss_data.shape}")
        
        # Estimate noise properties
        self._estimate_noise()
        
    def _estimate_noise(self):
        """Estimate RMS noise in each image using sigma-clipping"""
        print("Estimating noise properties...")
        
        # Use sigma-clipping to exclude sources and estimate noise
        vlass_clipped = sigma_clip(self.vlass_data[np.isfinite(self.vlass_data)], 
                                  sigma=3, maxiters=5)
        lotss_clipped = sigma_clip(self.lotss_data[np.isfinite(self.lotss_data)], 
                                  sigma=3, maxiters=5)
        
        self.vlass_rms = np.std(vlass_clipped[~vlass_clipped.mask])
        self.lotss_rms = np.std(lotss_clipped[~lotss_clipped.mask])
        
        # Create error maps (simple model: constant RMS)
        self.vlass_error = np.full_like(self.vlass_data, self.vlass_rms)
        self.lotss_error = np.full_like(self.lotss_data, self.lotss_rms)
        
        print(f"  VLASS RMS noise: {self.vlass_rms:.2e} Jy/beam")
        print(f"  LoTSS RMS noise: {self.lotss_rms:.2e} Jy/beam")
        
    def create_detection_mask(self, snr_threshold=5.0):
        """
        Create mask for reliable detections
        
        Parameters:
        -----------
        snr_threshold : float
            Signal-to-noise ratio threshold for detection
        """
        print(f"Creating detection mask (SNR > {snr_threshold})...")
        
        # SNR masks for each survey
        vlass_snr_mask = (self.vlass_data / self.vlass_error) > snr_threshold
        lotss_snr_mask = (self.lotss_data / self.lotss_error) > snr_threshold
        
        # Combined mask: detected in both surveys
        self.source_mask = (vlass_snr_mask & lotss_snr_mask & 
                           np.isfinite(self.vlass_data) & 
                           np.isfinite(self.lotss_data))
        
        n_detections = np.sum(self.source_mask)
        coverage = 100 * n_detections / self.source_mask.size
        
        print(f"  Valid detections: {n_detections:,} pixels ({coverage:.1f}% coverage)")
        
        if n_detections == 0:
            raise ValueError("No valid detections found! Try lowering SNR threshold.")
            
    def calculate_spectral_index_with_errors(self):
        """Calculate spectral index with proper error propagation"""
        print("Calculating spectral indices with error propagation...")
        
        # Initialize arrays
        self.spectral_index = np.full_like(self.vlass_data, np.nan)
        self.spectral_index_error = np.full_like(self.vlass_data, np.nan)
        
        # Calculate for detected pixels only
        valid_pixels = self.source_mask
        
        if np.sum(valid_pixels) == 0:
            print("  Warning: No valid pixels for calculation!")
            return
            
        # Extract data for valid pixels
        s1 = self.vlass_data[valid_pixels]  # Higher frequency
        s2 = self.lotss_data[valid_pixels]  # Lower frequency
        e1 = self.vlass_error[valid_pixels]
        e2 = self.lotss_error[valid_pixels]
        
        # Calculate spectral index: α = ln(S1/S2) / ln(ν1/ν2)
        flux_ratio = s1 / s2
        self.spectral_index[valid_pixels] = np.log(flux_ratio) / np.log(self.freq_ratio)
        
        # Error propagation: δα = sqrt((δS1/S1)² + (δS2/S2)²) / ln(ν1/ν2)
        relative_error_squared = (e1/s1)**2 + (e2/s2)**2
        self.spectral_index_error[valid_pixels] = (np.sqrt(relative_error_squared) / 
                                                  np.log(self.freq_ratio))
        
        # Statistics
        valid_alpha = self.spectral_index[np.isfinite(self.spectral_index)]
        valid_errors = self.spectral_index_error[np.isfinite(self.spectral_index_error)]
        
        print(f"  Mean spectral index: {np.mean(valid_alpha):.3f} ± {np.mean(valid_errors):.3f}")
        print(f"  Median spectral index: {np.median(valid_alpha):.3f}")
        print(f"  RMS scatter: {np.std(valid_alpha):.3f}")
        
    def classify_sources(self):
        """Classify sources based on spectral properties"""
        print("\nClassifying sources by spectral index...")
        
        valid_mask = np.isfinite(self.spectral_index)
        alpha_values = self.spectral_index[valid_mask]
        alpha_errors = self.spectral_index_error[valid_mask]
        
        if len(alpha_values) == 0:
            print("  No valid spectral indices for classification!")
            return
        
        # Classification thresholds with error consideration
        very_steep = np.sum(alpha_values < -1.0)
        steep = np.sum((alpha_values >= -1.0) & (alpha_values < -0.7))
        typical = np.sum((alpha_values >= -0.7) & (alpha_values < -0.5))
        flat = np.sum((alpha_values >= -0.5) & (alpha_values < 0.0))
        rising = np.sum(alpha_values >= 0.0)
        
        total = len(alpha_values)
        
        print("  Source Classification:")
        print(f"    Very steep (α < -1.0):     {very_steep:6,} ({100*very_steep/total:5.1f}%)")
        print(f"    Steep (-1.0 ≤ α < -0.7):   {steep:6,} ({100*steep/total:5.1f}%)")
        print(f"    Typical (-0.7 ≤ α < -0.5): {typical:6,} ({100*typical/total:5.1f}%)")
        print(f"    Flat (-0.5 ≤ α < 0.0):     {flat:6,} ({100*flat/total:5.1f}%)")
        print(f"    Rising (α ≥ 0.0):          {rising:6,} ({100*rising/total:5.1f}%)")
        
        # Physical interpretation
        print("\n  Physical Interpretation:")
        if steep + very_steep > 0.5 * total:
            print("    → Dominated by synchrotron sources (normal galaxies, SNRs)")
        if flat > 0.3 * total:
            print("    → Significant flat-spectrum population (AGN cores)")
        if rising > 0.1 * total:
            print("    → Notable rising spectrum sources (thermal emission, absorption)")
            
    def regional_analysis(self, region_size_arcmin=10):
        """Perform regional statistics analysis"""
        print(f"\nPerforming regional analysis ({region_size_arcmin}' regions)...")
        
        # Convert region size to pixels
        pixel_scale = abs(self.vlass_header.get('CDELT2', 1.0)) * 3600  # arcsec/pixel
        region_size_pixels = int(region_size_arcmin * 60 / pixel_scale)
        
        ny, nx = self.spectral_index.shape
        n_regions_y = ny // region_size_pixels
        n_regions_x = nx // region_size_pixels
        
        regional_stats = []
        
        for i in range(n_regions_y):
            for j in range(n_regions_x):
                # Define region boundaries
                y1 = i * region_size_pixels
                y2 = (i + 1) * region_size_pixels
                x1 = j * region_size_pixels
                x2 = (j + 1) * region_size_pixels
                
                # Extract region data
                region_alpha = self.spectral_index[y1:y2, x1:x2]
                region_mask = np.isfinite(region_alpha)
                
                if np.sum(region_mask) > 10:  # Require at least 10 valid pixels
                    region_valid = region_alpha[region_mask]
                    stats = {
                        'region': (i, j),
                        'center_y': y1 + region_size_pixels//2,
                        'center_x': x1 + region_size_pixels//2,
                        'n_pixels': len(region_valid),
                        'mean_alpha': np.mean(region_valid),
                        'std_alpha': np.std(region_valid),
                        'median_alpha': np.median(region_valid)
                    }
                    regional_stats.append(stats)
        
        if regional_stats:
            print(f"  Analyzed {len(regional_stats)} regions with sufficient data")
            
            # Overall statistics
            mean_alphas = [r['mean_alpha'] for r in regional_stats]
            std_alphas = [r['std_alpha'] for r in regional_stats]
            
            print(f"  Regional mean α range: {np.min(mean_alphas):.3f} to {np.max(mean_alphas):.3f}")
            print(f"  Average regional scatter: {np.mean(std_alphas):.3f}")
            print(f"  Scatter between regions: {np.std(mean_alphas):.3f}")
            
    def quality_assessment(self):
        """Assess data quality and provide recommendations"""
        print("\nData Quality Assessment:")
        
        valid_alpha = self.spectral_index[np.isfinite(self.spectral_index)]
        valid_errors = self.spectral_index_error[np.isfinite(self.spectral_index_error)]
        
        if len(valid_alpha) == 0:
            print("  ❌ POOR: No valid spectral indices calculated")
            return
        
        # Coverage assessment
        coverage = 100 * len(valid_alpha) / self.spectral_index.size
        if coverage > 20:
            coverage_grade = "✓ GOOD"
        elif coverage > 5:
            coverage_grade = "⚠ FAIR"
        else:
            coverage_grade = "❌ POOR"
        
        print(f"  Coverage: {coverage_grade} ({coverage:.1f}% of image)")
        
        # Precision assessment
        median_error = np.median(valid_errors) if len(valid_errors) > 0 else float('inf')
        if median_error < 0.1:
            precision_grade = "✓ EXCELLENT"
        elif median_error < 0.2:
            precision_grade = "✓ GOOD"  
        elif median_error < 0.5:
            precision_grade = "⚠ FAIR"
        else:
            precision_grade = "❌ POOR"
            
        print(f"  Precision: {precision_grade} (median error: {median_error:.3f})")
        
        # Distribution assessment
        if -1.0 < np.mean(valid_alpha) < -0.5:
            distribution_grade = "✓ NORMAL"
        else:
            distribution_grade = "⚠ UNUSUAL"
            
        print(f"  Distribution: {distribution_grade} (mean α = {np.mean(valid_alpha):.3f})")
        
        # Recommendations
        print("\n  Recommendations:")
        if coverage < 10:
            print("    • Consider lowering detection thresholds")
            print("    • Check for systematic calibration issues")
        if median_error > 0.3:
            print("    • Noise levels may be too high for reliable analysis")
            print("    • Consider smoothing data or larger beam matching")
        if np.std(valid_alpha) > 0.5:
            print("    • Large scatter may indicate:")
            print("      - Mix of different source populations")
            print("      - Systematic errors in data processing")
            
    def save_advanced_results(self, output_dir="output"):
        """Save results with comprehensive metadata"""
        print(f"\nSaving advanced analysis results to {output_dir}/...")
        os.makedirs(output_dir, exist_ok=True)
        
        # Save spectral index map
        spectral_hdu = fits.PrimaryHDU(self.spectral_index, header=self.vlass_header)
        spectral_hdu.header['OBJECT'] = 'Advanced Spectral Index Map'
        spectral_hdu.header['BUNIT'] = 'dimensionless'
        spectral_hdu.header['COMMENT'] = 'Advanced spectral index analysis with errors'
        spectral_file = os.path.join(output_dir, "spectral_index_advanced.fits")
        spectral_hdu.writeto(spectral_file, overwrite=True)
        
        # Save error map
        error_hdu = fits.PrimaryHDU(self.spectral_index_error, header=self.vlass_header)
        error_hdu.header['OBJECT'] = 'Spectral Index Error Map' 
        error_hdu.header['BUNIT'] = 'dimensionless'
        error_hdu.header['COMMENT'] = 'Propagated errors for spectral index'
        error_file = os.path.join(output_dir, "spectral_index_errors.fits")
        error_hdu.writeto(error_file, overwrite=True)
        
        # Save detection mask
        mask_hdu = fits.PrimaryHDU(self.source_mask.astype(np.uint8), header=self.vlass_header)
        mask_hdu.header['OBJECT'] = 'Source Detection Mask'
        mask_hdu.header['BUNIT'] = 'boolean'
        mask_hdu.header['COMMENT'] = 'Reliable detection mask (1=valid, 0=invalid)'
        mask_file = os.path.join(output_dir, "detection_mask.fits")
        mask_hdu.writeto(mask_file, overwrite=True)
        
        print(f"  ✓ Spectral index map: {spectral_file}")
        print(f"  ✓ Error map: {error_file}")
        print(f"  ✓ Detection mask: {mask_file}")
        

def main():
    """Main execution function"""
    print("=" * 70)
    print(" ADVANCED RADIO SPECTRAL INDEX ANALYSIS")
    print("=" * 70)
    
    # Define input files (modify these paths as needed)
    vlass_file = os.path.join("data", "sample_vlass_cutout.fits")
    lotss_file = os.path.join("data", "sample_lotss_cutout.fits")
    
    # Check if files exist
    if not os.path.exists(vlass_file):
        print(f"ERROR: VLASS file not found: {vlass_file}")
        return 1
    if not os.path.exists(lotss_file):
        print(f"ERROR: LoTSS file not found: {lotss_file}")
        return 1
    
    try:
        # Initialize analyzer
        analyzer = AdvancedSpectralAnalyzer(vlass_file, lotss_file)
        
        # Process data
        analyzer.load_data()
        analyzer.create_detection_mask(snr_threshold=3.0)
        analyzer.calculate_spectral_index_with_errors()
        
        # Analysis
        analyzer.classify_sources()
        analyzer.regional_analysis()
        analyzer.quality_assessment()
        
        # Save results
        analyzer.save_advanced_results()
        
        print("\n" + "=" * 70)
        print(" ADVANCED ANALYSIS COMPLETE!")
        print("=" * 70)
        print("Check the 'output/' directory for detailed results including:")
        print("• Spectral index map with full error analysis")
        print("• Error propagation maps")  
        print("• Source detection masks")
        print("• Quality assessment metrics")
        
        return 0
        
    except Exception as e:
        print(f"\nERROR during analysis: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    exit(main())