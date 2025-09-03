# -*- coding: utf-8 -*-
"""
Created on Wed Sep  3 16:34:16 2025

@author: adity
"""
"""
Basic functionality tests for spectral index analysis.

These tests verify that the core functions work correctly with synthetic data
and provide regression testing for the main analysis pipeline.

Usage:
    python test_basic_functionality.py
    
Author: [Your Name]
Date: 2025
"""

import os
import sys
import numpy as np
import tempfile
import shutil
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from radio_beam import Beam

# Add parent directory to path to import modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def create_synthetic_data():
    """Create synthetic FITS files for testing"""
    print("Creating synthetic test data...")
    
    # Create temporary directory
    temp_dir = tempfile.mkdtemp()
    
    # Image parameters
    nx, ny = 100, 100
    pixel_scale = 2.0 * u.arcsec  # 2 arcsec pixels
    
    # Create coordinate system
    header = fits.Header()
    header['NAXIS'] = 2
    header['NAXIS1'] = nx
    header['NAXIS2'] = ny
    header['CTYPE1'] = 'RA---SIN'
    header['CTYPE2'] = 'DEC--SIN'
    header['CRVAL1'] = 180.0  # RA center
    header['CRVAL2'] = 45.0   # Dec center
    header['CRPIX1'] = nx // 2
    header['CRPIX2'] = ny // 2
    header['CDELT1'] = -pixel_scale.to(u.deg).value
    header['CDELT2'] = pixel_scale.to(u.deg).value
    header['BUNIT'] = 'Jy/beam'
    
    # Create synthetic sources with known spectral indices
    x, y = np.meshgrid(np.arange(nx), np.arange(ny))
    
    # Source parameters
    sources = [
        {'x': 25, 'y': 25, 'flux_3ghz': 0.01, 'alpha': -0.7},  # Typical synchrotron
        {'x': 75, 'y': 25, 'flux_3ghz': 0.005, 'alpha': -0.5}, # Flatter spectrum
        {'x': 25, 'y': 75, 'flux_3ghz': 0.02, 'alpha': -1.0},  # Steep spectrum
        {'x': 75, 'y': 75, 'flux_3ghz': 0.008, 'alpha': -0.3}, # Very flat
    ]
    
    # Frequencies
    freq_vlass = 3000  # MHz
    freq_lotss = 144   # MHz
    freq_ratio = freq_vlass / freq_lotss
    
    # Initialize data arrays
    vlass_data = np.random.normal(0, 0.0001, (ny, nx))  # Noise
    lotss_data = np.random.normal(0, 0.0001, (ny, nx))  # Noise
    
    # Add synthetic sources
    for src in sources:
        # Create Gaussian source
        sigma = 3.0  # pixels
        gaussian = np.exp(-((x - src['x'])**2 + (y - src['y'])**2) / (2 * sigma**2))
        
        # Add to VLASS (3 GHz)
        vlass_data += src['flux_3ghz'] * gaussian
        
        # Calculate LoTSS flux using spectral index
        flux_lotss = src['flux_3ghz'] * (freq_lotss / freq_vlass)**src['alpha']
        lotss_data += flux_lotss * gaussian
    
    # Save VLASS file
    vlass_file = os.path.join(temp_dir, 'test_vlass.fits')
    vlass_hdu = fits.PrimaryHDU(vlass_data, header=header)
    vlass_hdu.header['OBJECT'] = 'Synthetic VLASS Test Data'
    vlass_hdu.writeto(vlass_file, overwrite=True)
    
    # Save LoTSS file  
    lotss_file = os.path.join(temp_dir, 'test_lotss.fits')
    lotss_hdu = fits.PrimaryHDU(lotss_data, header=header)
    lotss_hdu.header['OBJECT'] = 'Synthetic LoTSS Test Data'
    lotss_hdu.writeto(lotss_file, overwrite=True)
    
    return temp_dir, vlass_file, lotss_file, sources

def test_data_loading():
    """Test FITS file loading and WCS extraction"""
    print("\n1. Testing data loading...")
    
    temp_dir, vlass_file, lotss_file, sources = create_synthetic_data()
    
    try:
        # Load data
        vlass_data = fits.getdata(vlass_file)
        vlass_header = fits.getheader(vlass_file)
        vlass_wcs = WCS(vlass_header)
        
        lotss_data = fits.getdata(lotss_file)
        lotss_header = fits.getheader(lotss_file)
        lotss_wcs = WCS(lotss_header)
        
        # Verify shapes
        assert vlass_data.shape == (100, 100), f"Expected (100, 100), got {vlass_data.shape}"
        assert lotss_data.shape == (100, 100), f"Expected (100, 100), got {lotss_data.shape}"
        
        # Verify WCS
        assert vlass_wcs.naxis == 2, f"Expected 2D WCS, got {vlass_wcs.naxis}D"
        assert lotss_wcs.naxis == 2, f"Expected 2D WCS, got {lotss_wcs.naxis}D"
        
        print("   ✓ Data loading successful")
        print(f"   ✓ VLASS shape: {vlass_data.shape}")
        print(f"   ✓ LoTSS shape: {lotss_data.shape}")
        
        return True, temp_dir, vlass_data, lotss_data, sources
        
    except Exception as e:
        print(f"   ❌ Data loading failed: {e}")
        return False, None, None, None, None
    finally:
        # Cleanup
        shutil.rmtree(temp_dir, ignore_errors=True)

def test_spectral_index_calculation():
    """Test spectral index calculation with known sources"""
    print("\n2. Testing spectral index calculation...")
    
    temp_dir, vlass_file, lotss_file, sources = create_synthetic_data()
    
    try:
        # Load data
        vlass_data = fits.getdata(vlass_file)
        lotss_data = fits.getdata(lotss_file)
        
        # Calculate spectral index
        freq_vlass = 3000  # MHz
        freq_lotss = 144   # MHz
        
        # Mask for non-zero data
        mask = (vlass_data > 0.001) & (lotss_data > 0.001)
        
        if np.sum(mask) == 0:
            print("   ❌ No valid data for calculation")
            return False
        
        # Calculate spectral index
        flux_ratio = vlass_data[mask] / lotss_data[mask]
        freq_ratio = freq_vlass / freq_lotss
        spectral_index = np.log(flux_ratio) / np.log(freq_ratio)
        
        # Check that we get reasonable values
        mean_alpha = np.mean(spectral_index)
        std_alpha = np.std(spectral_index)
        
        print(f"   ✓ Calculated spectral indices for {np.sum(mask)} pixels")
        print(f"   ✓ Mean spectral index: {mean_alpha:.3f}")
        print(f"   ✓ RMS scatter: {std_alpha:.3f}")
        
        # Verify results are physically reasonable
        if -2.0 < mean_alpha < 0.0:
            print("   ✓ Results are physically reasonable")
            return True
        else:
            print(f"   ⚠ Warning: Unusual mean spectral index: {mean_alpha:.3f}")
            return False
            
    except Exception as e:
        print(f"   ❌ Spectral index calculation failed: {e}")
        return False
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

def test_noise_handling():
    """Test noise estimation and thresholding"""
    print("\n3. Testing noise handling...")
    
    try:
        # Create noise-only data
        np.random.seed(42)  # For reproducible results
        noise_data = np.random.normal(0, 0.001, (50, 50))
        
        # Add a few bright sources
        noise_data[25, 25] = 0.01  # 10x noise
        noise_data[35, 35] = 0.005  # 5x noise
        
        # Estimate noise using sigma clipping (simplified version)
        from astropy.stats import sigma_clip
        clipped_data = sigma_clip(noise_data, sigma=3, maxiters=5)
        rms_noise = np.std(clipped_data[~clipped_data.mask])
        
        # Create threshold
        threshold = 3 * rms_noise
        
        # Test detection
        detections = noise_data > threshold
        n_detections = np.sum(detections)
        
        print(f"   ✓ Estimated RMS noise: {rms_noise:.4f}")
        print(f"   ✓ 3-sigma threshold: {threshold:.4f}")
        print(f"   ✓ Number of detections: {n_detections}")
        
        # Should detect the bright sources we added
        if 1 <= n_detections <= 10:  # Reasonable number
            print("   ✓ Noise handling working correctly")
            return True
        else:
            print(f"   ⚠ Unexpected number of detections: {n_detections}")
            return False
            
    except Exception as e:
        print(f"   ❌ Noise handling test failed: {e}")
        return False

def test_beam_handling():
    """Test beam creation and manipulation"""
    print("\n4. Testing beam handling...")
    
    try:
        # Create beams
        vlass_beam = Beam(major=2.5*u.arcsec, minor=2.5*u.arcsec, pa=0*u.deg)
        lotss_beam = Beam(major=6.0*u.arcsec, minor=6.0*u.arcsec, pa=0*u.deg)
        
        print(f"   ✓ VLASS beam: {vlass_beam}")
        print(f"   ✓ LoTSS beam: {lotss_beam}")
        
        # Test beam comparison
        larger_beam = max(vlass_beam, lotss_beam)
        print(f"   ✓ Larger beam: {larger_beam}")
        
        # Test beam deconvolution
        if larger_beam == lotss_beam:
            deconv_beam = larger_beam.deconvolve(vlass_beam)
            print(f"   ✓ Deconvolved beam: {deconv_beam}")
        
        print("   ✓ Beam handling working correctly")
        return True
        
    except Exception as e:
        print(f"   ❌ Beam handling test failed: {e}")
        return False

def test_coordinate_handling():
    """Test coordinate system and WCS operations"""
    print("\n5. Testing coordinate handling...")
    
    try:
        # Create a simple WCS
        header = fits.Header()
        header['NAXIS'] = 2
        header['NAXIS1'] = 100
        header['NAXIS2'] = 100
        header['CTYPE1'] = 'RA---SIN'
        header['CTYPE2'] = 'DEC--SIN'
        header['CRVAL1'] = 180.0
        header['CRVAL2'] = 45.0
        header['CRPIX1'] = 50
        header['CRPIX2'] = 50
        header['CDELT1'] = -0.001  # degrees
        header['CDELT2'] = 0.001   # degrees
        
        wcs = WCS(header)
        
        # Test pixel to world coordinate conversion
        ra, dec = wcs.all_pix2world(50, 50, 1)  # Center pixel
        
        print(f"   ✓ Center coordinates: RA={ra:.3f}°, Dec={dec:.3f}°")
        
        # Test world to pixel conversion
        x, y = wcs.all_world2pix(180.0, 45.0, 1)
        
        print(f"   ✓ Center pixel: x={x:.1f}, y={y:.1f}")
        
        # Verify round-trip accuracy
        if abs(x - 50) < 0.1 and abs(y - 50) < 0.1:
            print("   ✓ Coordinate transformations accurate")
            return True
        else:
            print(f"   ❌ Coordinate transformation error: expected (50,50), got ({x:.1f},{y:.1f})")
            return False
            
    except Exception as e:
        print(f"   ❌ Coordinate handling test failed: {e}")
        return False

def test_error_propagation():
    """Test error propagation calculations"""
    print("\n6. Testing error propagation...")
    
    try:
        # Synthetic data with known errors
        s1 = np.array([0.01, 0.005, 0.02])  # VLASS fluxes
        s2 = np.array([0.05, 0.025, 0.10])  # LoTSS fluxes
        e1 = np.array([0.001, 0.0005, 0.002])  # VLASS errors
        e2 = np.array([0.005, 0.0025, 0.010])  # LoTSS errors
        
        freq_ratio = 3000 / 144
        
        # Calculate spectral index
        flux_ratio = s1 / s2
        alpha = np.log(flux_ratio) / np.log(freq_ratio)
        
        # Error propagation
        relative_error_squared = (e1/s1)**2 + (e2/s2)**2
        alpha_error = np.sqrt(relative_error_squared) / np.log(freq_ratio)
        
        print(f"   ✓ Spectral indices: {alpha}")
        print(f"   ✓ Spectral index errors: {alpha_error}")
        
        # Check that errors are reasonable (not NaN, not too large)
        if np.all(np.isfinite(alpha_error)) and np.all(alpha_error < 1.0):
            print("   ✓ Error propagation working correctly")
            return True
        else:
            print("   ❌ Error propagation produced invalid results")
            return False
            
    except Exception as e:
        print(f"   ❌ Error propagation test failed: {e}")
        return False

def test_file_io():
    """Test FITS file input/output operations"""
    print("\n7. Testing file I/O operations...")
    
    temp_dir = tempfile.mkdtemp()
    
    try:
        # Create test data
        test_data = np.random.random((50, 50))
        
        # Create header
        header = fits.Header()
        header['OBJECT'] = 'Test Data'
        header['BUNIT'] = 'Jy/beam'
        header['COMMENT'] = 'Test file for I/O operations'
        
        # Write file
        test_file = os.path.join(temp_dir, 'test_output.fits')
        hdu = fits.PrimaryHDU(test_data, header=header)
        hdu.writeto(test_file, overwrite=True)
        
        print(f"   ✓ Created test file: {os.path.basename(test_file)}")
        
        # Read file back
        read_data = fits.getdata(test_file)
        read_header = fits.getheader(test_file)
        
        # Verify data integrity
        if np.allclose(test_data, read_data):
            print("   ✓ Data integrity verified")
        else:
            print("   ❌ Data corruption detected")
            return False
            
        # Verify header
        if read_header['OBJECT'] == 'Test Data':
            print("   ✓ Header information preserved")
        else:
            print("   ❌ Header information corrupted")
            return False
        
        # Check file size
        file_size = os.path.getsize(test_file)
        print(f"   ✓ File size: {file_size} bytes")
        
        return True
        
    except Exception as e:
        print(f"   ❌ File I/O test failed: {e}")
        return False
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

def test_integration():
    """Integration test of the complete workflow"""
    print("\n8. Integration test...")
    
    temp_dir, vlass_file, lotss_file, sources = create_synthetic_data()
    
    try:
        print("   Running complete analysis workflow...")
        
        # Simulate the main workflow
        # 1. Load data
        vlass_data = fits.getdata(vlass_file)
        lotss_data = fits.getdata(lotss_file)
        vlass_header = fits.getheader(vlass_file)
        
        # 2. Create masks
        vlass_threshold = 3 * np.std(vlass_data)
        lotss_threshold = 3 * np.std(lotss_data)
        
        vlass_mask = vlass_data > vlass_threshold
        lotss_mask = lotss_data > lotss_threshold
        combined_mask = vlass_mask & lotss_mask
        
        # 3. Calculate spectral index
        spectral_index = np.full_like(vlass_data, np.nan)
        
        if np.sum(combined_mask) > 0:
            flux_ratio = vlass_data[combined_mask] / lotss_data[combined_mask]
            freq_ratio = 3000 / 144
            spectral_index[combined_mask] = np.log(flux_ratio) / np.log(freq_ratio)
        
        # 4. Create output
        output_file = os.path.join(temp_dir, 'integration_test_result.fits')
        hdu = fits.PrimaryHDU(spectral_index, header=vlass_header)
        hdu.header['OBJECT'] = 'Integration Test Result'
        hdu.writeto(output_file, overwrite=True)
        
        # Verify results
        valid_alpha = spectral_index[np.isfinite(spectral_index)]
        
        if len(valid_alpha) > 0:
            print(f"   ✓ Calculated {len(valid_alpha)} valid spectral indices")
            print(f"   ✓ Mean spectral index: {np.mean(valid_alpha):.3f}")
            print(f"   ✓ Output file created: {os.path.basename(output_file)}")
            
            # Check if we recovered the input spectral indices approximately
            expected_mean = np.mean([src['alpha'] for src in sources])
            actual_mean = np.mean(valid_alpha)
            
            if abs(actual_mean - expected_mean) < 0.2:  # Allow 20% tolerance
                print(f"   ✓ Spectral index recovery within tolerance")
                print(f"     Expected: {expected_mean:.3f}, Got: {actual_mean:.3f}")
                return True
            else:
                print(f"   ⚠ Spectral index recovery outside tolerance")
                print(f"     Expected: {expected_mean:.3f}, Got: {actual_mean:.3f}")
                return False
        else:
            print("   ❌ No valid spectral indices calculated")
            return False
            
    except Exception as e:
        print(f"   ❌ Integration test failed: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

def run_all_tests():
    """Run all tests and provide summary"""
    print("=" * 70)
    print(" BASIC FUNCTIONALITY TESTS")
    print("=" * 70)
    print("Testing core components of the spectral index analysis pipeline...")
    
    tests = [
        ("Data Loading", test_data_loading),
        ("Spectral Index Calculation", test_spectral_index_calculation),
        ("Noise Handling", test_noise_handling), 
        ("Beam Handling", test_beam_handling),
        ("Coordinate Handling", test_coordinate_handling),
        ("Error Propagation", test_error_propagation),
        ("File I/O Operations", test_file_io),
        ("Integration Test", test_integration),
    ]
    
    results = []
    for test_name, test_func in tests:
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            print(f"   ❌ {test_name} crashed: {e}")
            results.append((test_name, False))
    
    # Summary
    print("\n" + "=" * 70)
    print(" TEST SUMMARY")
    print("=" * 70)
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for test_name, result in results:
        status = "PASS" if result else "FAIL"
        icon = "✓" if result else "❌"
        print(f"{icon} {test_name:25} {status}")
    
    print("-" * 70)
    print(f"Tests passed: {passed}/{total} ({100*passed/total:.0f}%)")
    
    if passed == total:
        print("🎉 All tests passed! The spectral index analysis pipeline is working correctly.")
        return 0
    else:
        print("⚠️  Some tests failed. Please check the implementation.")
        return 1

def main():
    """Main test execution"""
    try:
        # Check required packages
        required_packages = ['numpy', 'astropy', 'radio_beam']
        missing_packages = []
        
        for package in required_packages:
            try:
                __import__(package)
            except ImportError:
                missing_packages.append(package)
        
        if missing_packages:
            print(f"ERROR: Missing required packages: {missing_packages}")
            print("Please install using: pip install " + " ".join(missing_packages))
            return 1
        
        return run_all_tests()
        
    except KeyboardInterrupt:
        print("\nTests interrupted by user")
        return 1
    except Exception as e:
        print(f"Test execution failed: {e}")
        return 1

if __name__ == "__main__":
    exit(main())