
import os
import numpy as np
import warnings
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from astropy.stats import mad_std
from radio_beam import Beam
from reproject import reproject_interp

# Suppress warnings
warnings.filterwarnings('ignore')

class RadioMap:
    def __init__(self, filepath, name="Map"):
        self.filepath = filepath
        self.name = name
        self.data = None
        self.header = None
        self.wcs = None
        self.beam = None
        self.freq = None
        
        self.load_data()

    def load_data(self):
        if not os.path.exists(self.filepath):
            raise FileNotFoundError(f"File not found: {self.filepath}")
        
        with fits.open(self.filepath) as hdul:
            hdu = hdul[0] if len(hdul) > 0 and hdul[0].data is not None else hdul[1]
            self.header = hdu.header
            self.data = np.squeeze(hdu.data)
            self.wcs = WCS(self.header).celestial
            
        print(f"\n--- Loaded {self.name} ---")
        self._get_frequency()
        self._get_beam()
        self._check_units()

    def _get_frequency(self):
        keys = ['RESTFRQ', 'FREQ', 'CRVAL3'] 
        found = False
        for k in keys:
            if k in self.header:
                val = self.header[k]
                if val > 1e6:
                    self.freq = val * u.Hz
                    print(f"  Freq: {self.freq.to(u.MHz):.2f}")
                    found = True
                    break
        if not found:
            val = float(input(f"  [!] Enter Frequency for {self.name} (MHz): "))
            self.freq = val * u.MHz

    def _get_beam(self):
        try:
            self.beam = Beam.from_fits_header(self.header)
            print(f"  Beam: {self.beam}")
        except:
            print(f"  [!] Beam info missing in {self.name}.")
            bmaj = float(input(f"  >>> Enter Major Axis (arcsec): "))
            bmin = float(input(f"  >>> Enter Minor Axis (arcsec): "))
            bpa = float(input(f"  >>> Enter PA (deg): "))
            self.beam = Beam(major=bmaj*u.arcsec, minor=bmin*u.arcsec, pa=bpa*u.deg)

    def _check_units(self):
        """Detects mJy vs Jy and converts to Jy if needed."""
        unit_str = self.header.get('BUNIT', '').lower()
        
        if 'mjy' in unit_str:
            print("  [Correction] Converting mJy -> Jy")
            self.data = self.data / 1000.0
            self.header['BUNIT'] = 'Jy/beam'
        elif 'jy' in unit_str:
            pass 
        else:
            print("  [!] Unit unknown. Assuming Jy/beam.")
            is_mjy = input("  Is this map in mJy/beam? (y/n): ").lower()
            if is_mjy == 'y':
                self.data = self.data / 1000.0

def convolve_to_common(source_map, target_beam):
    """Convolves map AND scales pixel flux for the new beam size."""
    # Check if convolution is needed
    if source_map.beam == target_beam:
        print(f"  {source_map.name} already matches target beam.")
        return source_map.data

    print(f"\n--- Convolving {source_map.name} ---")
    print(f"  From: {source_map.beam}")
    print(f"  To:   {target_beam}")
    
    try:
        # 1. Calculate Kernel
        kernel = target_beam.deconvolve(source_map.beam)
        kernel_pix = kernel.as_kernel(source_map.wcs.proj_plane_pixel_area()**0.5)
        
        # 2. Convolve (FFT)
        from astropy.convolution import convolve_fft
        convolved_data = convolve_fft(source_map.data, kernel_pix, allow_huge=True)
        
        # 3. Apply Beam Area Scaling
        source_area = source_map.beam.sr.value
        target_area = target_beam.sr.value
        scale_factor = target_area / source_area
        
        print(f"  [Physics] Scaling Factor: {scale_factor:.4f}")
        return convolved_data * scale_factor
        
    except ValueError:
        print("  [CRITICAL ERROR] Deconvolution failed despite Common Beam logic.")
        return source_map.data

def calculate_spectral_index_workflow():
    print("========================================================")
    print("   Spectral Index Pipeline (Common Resolution Method)   ")
    print("========================================================")

    # 1. INPUTS
    path1 = input("Enter path to FITS File 1: ").strip().replace('"', '')
    path2 = input("Enter path to FITS File 2: ").strip().replace('"', '')

    map1 = RadioMap(path1, "Map 1")
    map2 = RadioMap(path2, "Map 2")

    # 2. DEFINE COMMON BEAM (Worst Case)
    # We find the largest axis in EITHER map and create a circular beam of that size.
    # This guarantees both maps can be convolved to this resolution.
    max_axis_1 = max(map1.beam.major, map1.beam.minor)
    max_axis_2 = max(map2.beam.major, map2.beam.minor)
    
    # Add 1% buffer to ensure kernel is well-defined
    common_size = max(max_axis_1, max_axis_2) * 1.01 
    
    common_beam = Beam(major=common_size, minor=common_size, pa=0*u.deg)
    print(f"\n[Resolution] Defined Common Circular Beam: {common_beam.major.to(u.arcsec):.2f}")

    # 3. CONVOLVE BOTH MAPS
    data1_conv = convolve_to_common(map1, common_beam)
    data2_conv = convolve_to_common(map2, common_beam)

    # 4. REGRIDDING (Map 1 -> Map 2 Grid)
    # Since both are now at the same resolution (common_beam), we can regrid safely.
    print(f"\n--- Regridding Map 1 to Map 2 Grid (Order=3) ---")
    
    data1_aligned, footprint = reproject_interp(
        (data1_conv, map1.wcs),
        map2.wcs,
        shape_out=map2.data.shape,
        order=3
    )
    
    # Now we use data1_aligned (Map 1 on Grid 2) and data2_conv (Map 2 on Grid 2)

    # 5. NOISE & MASKING
    print("\n--- Noise Characterization (MAD) ---")
    rms_1 = mad_std(data1_aligned, ignore_nan=True)
    rms_2 = mad_std(data2_conv, ignore_nan=True)
    print(f"  Map 1 (Proc) RMS: {rms_1:.4e}")
    print(f"  Map 2 (Proc) RMS: {rms_2:.4e}")

    sigma_thresh = float(input("\n>>> Enter Sigma Threshold (e.g., 3.0): "))
    mask = (data1_aligned > sigma_thresh*rms_1) & (data2_conv > sigma_thresh*rms_2)
    
    # 6. CALCULATION
    print("\n--- Calculating Alpha ---")
    S1 = data1_aligned[mask]
    S2 = data2_conv[mask]
    v1 = map1.freq.to(u.Hz).value
    v2 = map2.freq.to(u.Hz).value
    
    alpha_map = np.full_like(map2.data, np.nan)
    
    with np.errstate(invalid='ignore'):
        alpha_vals = np.log10(S1 / S2) / np.log10(v1 / v2)
        alpha_map[mask] = alpha_vals

    # 7. SAVE
    print("\n--- Saving ---")
    out_header = map2.header.copy()
    
    # Update header to reflect new common beam
    try:
        out_header.update(common_beam.to_header_keywords())
        out_header['HISTORY'] = f'Convolved to common beam: {common_size.to(u.arcsec):.2f}'
    except: pass
    
    fits.writeto('spidx_common.fits', alpha_map, out_header, overwrite=True)
    print("  Saved: spidx_common.fits")
    print("  Done.")

if __name__ == "__main__":
    calculate_spectral_index_workflow()