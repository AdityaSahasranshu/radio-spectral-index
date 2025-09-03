# -*- coding: utf-8 -*-
"""
Created on Fri May 23 23:04:11 2025

@author: adity
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 00:18:14 2024

@author: adity
"""

import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from radio_beam import Beam
from astropy.convolution import convolve
from reproject import reproject_interp

def get_pixel_scale(header, wcs):
    """
    Extract pixel scale from FITS header, handling different WCS formats
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

# Load FITS files
vlass_file = r"D:\SOMU\SALTY LOTSS\FITS 1\P225+47\New Folder (11)\UNK.AUTH-P-VLASS3.fits"
lotss_file = r"D:\SOMU\SALTY LOTSS\FITS 1\P225+47\New Folder (11)\astron.fits"

print("Loading FITS files...")
vlass_data = fits.getdata(vlass_file)
vlass_header = fits.getheader(vlass_file)
vlass_wcs = WCS(vlass_header, naxis=2)

lotss_data = fits.getdata(lotss_file)
lotss_header = fits.getheader(lotss_file)
lotss_wcs = WCS(lotss_header, naxis=2)

print(f"VLASS data shape: {vlass_data.shape}")
print(f"LoTSS data shape: {lotss_data.shape}")

# Handle extra dimensions
if vlass_data.ndim > 2:
    vlass_data = vlass_data[0, 0]
    print("Reduced VLASS data to 2D")

if lotss_data.ndim > 2:
    lotss_data = lotss_data[0]
    print("Reduced LoTSS data to 2D")

print(f"Final VLASS data shape: {vlass_data.shape}")
print(f"Final LoTSS data shape: {lotss_data.shape}")

# Reproject LoTSS to VLASS grid
print("Reprojecting LoTSS to VLASS grid...")
lotss_reprojected, _ = reproject_interp((lotss_data, lotss_wcs), vlass_wcs, vlass_data.shape)

# Set beam information
vlass_beam = Beam(major=2.5*u.arcsec, minor=2.5*u.arcsec, pa=0*u.deg)
lotss_beam = Beam(major=5*u.arcsec, minor=5*u.arcsec, pa=0*u.deg)

# Convolve to common resolution
target_beam = max(vlass_beam, lotss_beam)
print(f"Target beam: {target_beam}")

def convolve_if_needed(data, original_beam, target_beam, pixel_scale):
    if target_beam > original_beam:
        print(f"Convolving from {original_beam} to {target_beam}")
        kernel = target_beam.deconvolve(original_beam).as_kernel(pixel_scale)
        return convolve(data, kernel)
    else:
        print("No convolution needed")
    return data

# Get pixel scale using the improved function
vlass_pixel_scale = get_pixel_scale(vlass_header, vlass_wcs)
print(f"VLASS pixel scale: {vlass_pixel_scale}")

print("Convolving data to common resolution...")
vlass_conv = convolve_if_needed(vlass_data, vlass_beam, target_beam, vlass_pixel_scale)
lotss_conv = convolve_if_needed(lotss_reprojected, lotss_beam, target_beam, vlass_pixel_scale)

# Apply noise threshold
vlass_threshold = 0.00021
lotss_threshold = 0.00024

print("Applying noise thresholds...")
vlass_conv[vlass_conv < vlass_threshold] = 0
lotss_conv[lotss_conv < lotss_threshold] = 0

# Calculate spectral index
vlass_freq = 3000  # MHz
lotss_freq = 144  # MHz

print("Calculating spectral index...")
with np.errstate(divide='ignore', invalid='ignore'):
    spectral_index = np.log(vlass_conv / lotss_conv) / np.log(vlass_freq / lotss_freq)
    spectral_index[np.isnan(spectral_index)] = 0
    spectral_index[np.isinf(spectral_index)] = 0  # Also handle infinite values

# Print some statistics
valid_pixels = spectral_index[spectral_index != 0]
if len(valid_pixels) > 0:
    print(f"Spectral index statistics:")
    print(f"  Valid pixels: {len(valid_pixels)}")
    print(f"  Mean: {np.mean(valid_pixels):.3f}")
    print(f"  Median: {np.median(valid_pixels):.3f}")
    print(f"  Std: {np.std(valid_pixels):.3f}")
    print(f"  Range: {np.min(valid_pixels):.3f} to {np.max(valid_pixels):.3f}")
else:
    print("Warning: No valid spectral index values calculated!")

# Create a new FITS file
print("Creating output FITS file...")
hdu = fits.PrimaryHDU(spectral_index, header=vlass_header)

# Add some metadata to the header
hdu.header['OBJECT'] = 'Spectral Index Map'
hdu.header['BUNIT'] = 'dimensionless'
hdu.header['COMMENT'] = f'Spectral index between {vlass_freq} MHz and {lotss_freq} MHz'
hdu.header['VFREQ'] = (vlass_freq, 'VLASS frequency in MHz')
hdu.header['LFREQ'] = (lotss_freq, 'LoTSS frequency in MHz')
hdu.header['VTHRESH'] = (vlass_threshold, 'VLASS noise threshold')
hdu.header['LTHRESH'] = (lotss_threshold, 'LoTSS noise threshold')

hdu.writeto('spectral_index.fits', overwrite=True)

print("Spectral index FITS file created: spectral_index.fits")