"""
===================================================================================
|                         FIREFLY — Full Spectral Fitting                         |
===================================================================================
 <> Module: Local_fits_create.py

 <> Author:
	- Samuel Helps ~ <samuel.helps.sh__at__gmail.com>

 <> Purpose:
	- Creates FITS files of DESI EDR spectra for a specified range of galaxies.
___________________________________________________________________________________
|Institute of Cosmology and Gravitation - University of Portsmouth, Portsmouth, UK|
‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
"""
import os
import sys
import numpy as np
from astropy.io import fits
from astropy.table import Table
import requests
from tqdm import tqdm
import time
from datetime import timedelta

# =======================================
#              SETTINGS
# =======================================

# Directory of this script
script_dir = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.abspath(os.path.join(script_dir, "..", "..", ".."))

# Data source URLs
BASE_URL = "https://data.desi.lbl.gov/public/edr"
ZCATALOG_URL = f"{BASE_URL}/spectro/redux/fuji/zcatalog/zall-fuji.fits"
FASTSPEC_URL = f"{BASE_URL}/vac/edr/fastspecfit/fuji/v3.2/catalogs/zall-pix-fuji.fits"

# Output folder for downloads
OUTPUT_DIR = os.path.join(ROOT_DIR, "Data", "DESI", "DESI_EDR_raw(downloads)")
os.makedirs(OUTPUT_DIR, exist_ok=True)

#########################################
START_SPECTRUM = 445000  # Starting index
END_SPECTRUM = 450000   # Ending index 
#########################################
# =======================================

def download_file(url, output_path, desc=None):
    """Download required DESI data files with live progress bar"""
    try:
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            total_size = int(r.headers.get('content-length', 0))
            with open(output_path, 'wb') as f:
                with tqdm(total=total_size, unit='iB', unit_scale=True, desc=desc) as pbar:
                    for chunk in r.iter_content(chunk_size=8192):
                        size = f.write(chunk)
                        pbar.update(size)
        return True
    except requests.exceptions.RequestException as e:
        print(f"Error downloading {url}: {e}")
        return False

def get_snr(row):
    """
    Calculate SNR using template-based (S/N)² values summed over B,R,Z bands.
    Returns squared root of SNR² as final SNR value.
    """
    try:
        if row['BGS_TARGET'] > 0:          # BGS target
            snr2 = row['TSNR2_BGS']
        elif row['DESI_TARGET'] & (2**1):  # ELG bit
            snr2 = row['TSNR2_ELG']
        elif row['DESI_TARGET'] & (2**2):  # LRG bit
            snr2 = row['TSNR2_LRG']
        else:
            snr2_values = []
            for key in ['TSNR2_BGS', 'TSNR2_ELG', 'TSNR2_LRG']:
                if key in row.dtype.names and np.isfinite(row[key]):
                    snr2_values.append(row[key])
            if snr2_values:
                snr2 = np.mean(snr2_values)
            else:
                snr2 = row['TSNR2_GPBDARK']
        if np.isfinite(snr2) and snr2 > 0:
            return np.sqrt(snr2)
        else:
            return None
    except Exception as e:
        print(f"Error calculating SNR for TARGETID {row['TARGETID']}: {e}")
        return None

def get_survey_type(row):
    """
    Determine survey type from target flags and program name.
    Returns a string identifying the specific program and target type
    """
    program = row['PROGRAM'].strip()
    if row['BGS_TARGET'] > 0:
        return f"BGS_{program}".ljust(15)
    elif row['DESI_TARGET'] > 0:
        if row['DESI_TARGET'] & (2**1):
            return f"ELG_{program}".ljust(15)
        elif row['DESI_TARGET'] & (2**2):
            return f"LRG_{program}".ljust(15)
        elif row['DESI_TARGET'] & (2**3):
            return f"QSO_{program}".ljust(15)
    elif row['MWS_TARGET'] > 0:
        return f"MWS_{program}".ljust(15)
    elif row['SCND_TARGET'] > 0:
        return f"SCND_{program}".ljust(15)
    return f"OTHER_{program}".ljust(15)

# === Load fastspec catalog and index by TARGETID (dict for quick future use) ===
print("\nAttempting to load fastspec catalog...")
fastspec_path = os.path.join(OUTPUT_DIR, "fastspec-fuji.fits")
fastspec_dict_path = os.path.join(OUTPUT_DIR, "fastspec_dict.npy")

try:
    # Try to load pre-computed dictionary if it exists
    if os.path.exists(fastspec_dict_path):
        print("Loading pre-computed fastspec dictionary...")
        fastspec_dict = np.load(fastspec_dict_path, allow_pickle=True).item()
        print(f"Loaded fastspec dictionary with {len(fastspec_dict)} entries")
    else:
        print("Reading fastspec catalog and creating dictionary...")
        # Read the whole table but only extract needed columns
        with fits.open(fastspec_path) as hdul:
            data = hdul[1].data
            targetids = data['TARGETID']
            vdisps = data['VDISP']
            print("Creating fastspec dictionary...")
            fastspec_dict = {tid: v for tid, v in zip(targetids, vdisps)}
            print(f"Saving fastspec dictionary for future use...")
            np.save(fastspec_dict_path, fastspec_dict)
            print(f"Created fastspec dictionary with {len(fastspec_dict)} entries")

except Exception as e:
    print(f"Error reading fastspec catalog: {e}")
    sys.exit(1)

def get_velocity_dispersion_from_fastspec(targetid):
    """Retrieve velocity dispersion from fastspec catalog if available"""
    vdisp = fastspec_dict.get(targetid, None)
    return float(vdisp) if vdisp is not None and np.isfinite(vdisp) else None

def get_snr_from_fastspec(row):
    """(FAILSAFE) Calculate SNR using template-based SNR² values from fastspec"""
    try:
        # Get relevant SNR² values based on target type
        snr2_values = []
        if row['BGS_TARGET'] > 0:
            if np.isfinite(row['TSNR2_BGS']):
                snr2_values.append(row['TSNR2_BGS'])
        elif row['DESI_TARGET'] & (2**1):  # ELG bit
            if np.isfinite(row['TSNR2_ELG']):
                snr2_values.append(row['TSNR2_ELG'])
        elif row['DESI_TARGET'] & (2**2):  # LRG bit
            if np.isfinite(row['TSNR2_LRG']):
                snr2_values.append(row['TSNR2_LRG'])
        
        # Use mean of available SNR² values
        if snr2_values:
            snr2 = np.mean(snr2_values)
            return np.sqrt(snr2)  # Convert SNR² to SNR
            
    except Exception as e:
        print(f"Error calculating fastspec SNR for TARGETID {row['TARGETID']}: {e}")
    return None

def calculate_vdisp(wavelength, resolution):
    """(FAILSAFE) Calculate velocity dispersion from resolution matrix"""
    try:
        # Resolution matrix is now 2D (wavelength x resolution)
        diag_elements = resolution
        fwhm = np.median(diag_elements)
        
        # Convert FWHM to sigma (FWHM = 2.355 * sigma)
        sigma_lambda = fwhm / 2.355
        
        # Convert wavelength resolution to velocity resolution
        # dv/c = dλ/λ
        lambda_central = np.median(wavelength)
        vdisp = 3e5 * (sigma_lambda / lambda_central)  # km/s
        
        # Only check for unreasonably low or invalid vdisp
        if not np.isfinite(vdisp) or vdisp < 50.0:
            return 125.0  # Default value for invalid measurements
            
        return vdisp
    except Exception as e:
        print(f"Error calculating velocity dispersion: {e}")
        return 125.0

def process_galaxy(row, output_dir):
    """Process individual galaxy row to extract spectrum and metadata"""
    try:
        targetid = row['TARGETID']
        ra, dec = row['TARGET_RA'], row['TARGET_DEC']
        z = row['Z']
        pix = row['HEALPIX']
        pixgroup = pix // 100
        program = row['PROGRAM'].strip().lower()
        survey = row['SURVEY'].strip().lower()
        
        # Download and process coadd file
        coadd_url = (f"{BASE_URL}/spectro/redux/fuji/healpix/"
                     f"{survey}/{program}/{pixgroup}/{pix}/"
                     f"coadd-{survey}-{program}-{pix}.fits")
        coadd_file = os.path.join(output_dir, f"coadd-{pix}.fits")
        if not os.path.exists(coadd_file):
            if not download_file(coadd_url, coadd_file, f"Downloading {pix}"):
                return None
                
        with fits.open(coadd_file) as hdul:
            fibermap = hdul['FIBERMAP'].data
            match = np.where(fibermap['TARGETID'] == targetid)[0]
            if len(match) == 0:
                print(f"No matching fiber found for TARGETID {targetid}")
                return None
                
            fiber_idx = match[0]
            
            # Load wavelength data first
            b_wave = hdul['B_WAVELENGTH'].data
            r_wave = hdul['R_WAVELENGTH'].data
            z_wave = hdul['Z_WAVELENGTH'].data
            
            # Get velocity dispersion from fastspec first
            vdisp = get_velocity_dispersion_from_fastspec(targetid)
        
            # If fastspec value is invalid, calculate from resolution matrices
            if vdisp is None or not np.isfinite(vdisp):
                b_vdisp = calculate_vdisp(b_wave, hdul['B_RESOLUTION'].data[fiber_idx][5])
                r_vdisp = calculate_vdisp(r_wave, hdul['R_RESOLUTION'].data[fiber_idx][5])
                z_vdisp = calculate_vdisp(z_wave, hdul['Z_RESOLUTION'].data[fiber_idx][5])
                vdisp = np.median([b_vdisp, r_vdisp, z_vdisp])
            
            # Load flux and ivar data
            b_flux = hdul['B_FLUX'].data[fiber_idx]
            r_flux = hdul['R_FLUX'].data[fiber_idx]
            z_flux = hdul['Z_FLUX'].data[fiber_idx]
            
            b_ivar = hdul['B_IVAR'].data[fiber_idx]
            r_ivar = hdul['R_IVAR'].data[fiber_idx]
            z_ivar = hdul['Z_IVAR'].data[fiber_idx]
            
            # Combine and sort wavelength, flux, and ivar
            wave = np.concatenate([b_wave, r_wave, z_wave])
            wave_sort = np.argsort(wave)
            wave = wave[wave_sort]
            flux = np.concatenate([b_flux, r_flux, z_flux])[wave_sort]
            ivar = np.concatenate([b_ivar, r_ivar, z_ivar])[wave_sort]
            
            survey_type = get_survey_type(row)
            
            # Try template-based SNR first
            #snr = get_snr(row)
            snr = get_snr_from_fastspec(row)
            
            # Fall back to flux-based SNR if template method failed
            if snr is None:
                valid = (ivar > 0) & np.isfinite(flux) & np.isfinite(ivar)
                if np.any(valid):
                    snr = np.median(np.abs(flux[valid] * np.sqrt(ivar[valid])))
                else:
                    print(f"Both SNR calculations failed for TARGETID {targetid}")
                    return None

            return {
                'ID': targetid,
                'WAVE': wave,
                'FLUX': flux,
                'IVAR': ivar,
                'REDSHIFT': z,
                'RA': ra,
                'DEC': dec,
                'VDISP': vdisp,
                'SNR': snr,
                'SURVEY_TYPE': survey_type
            }

    except Exception as e:
        print(f"Error processing TARGETID {targetid}: {e}")
        return None

def main():
    print("Starting main execution...")
    edr_dir = os.path.join(ROOT_DIR, "Data", "DESI", "DESI_EDR_data")
    os.makedirs(edr_dir, exist_ok=True)
    base_name = 'DESI_EDR'
    output_file = os.path.join(edr_dir, f'{base_name}_{START_SPECTRUM}-{END_SPECTRUM}.fits')
    
    # Check for existing output file
    if os.path.exists(output_file):
        print(f"\nFound existing output file: {output_file}")
        response = input("Do you want to reprocess the data? (y/n): ").lower()
        if response != 'y':
            print("Using existing file. Exiting...")
            sys.exit(0)
    start_time = time.time()

    # Download catalog if not present
    print("Downloading catalog...")
    zcat_path = os.path.join(OUTPUT_DIR, "zall-fuji.fits")
    if not os.path.exists(zcat_path):
        if not download_file(ZCATALOG_URL, zcat_path, "Downloading catalog"):
            sys.exit(1)

    # Load catalog and filter for galaxies
    print("Loading catalog...")
    zcat = Table.read(zcat_path, hdu=1)
    galaxies = zcat[zcat['SPECTYPE'] == 'GALAXY']
    
    # Select the specified range
    galaxies = galaxies[START_SPECTRUM:END_SPECTRUM]
    print(f"Processing galaxies from index {START_SPECTRUM} to {END_SPECTRUM}")
    print(f"Total galaxies in selection: {len(galaxies)}")

    # Process each galaxy and collect results
    combined_rows = []
    for row in tqdm(galaxies, desc="Processing galaxies"):
        result = process_galaxy(row, OUTPUT_DIR)
        if result:
            combined_rows.append(result)

    if not combined_rows:
        print("No valid spectra found!")
        sys.exit(1)

    print(f"Writing {len(combined_rows)} spectra...")
    n_wave = len(combined_rows[0]['WAVE'])

    # Convert lists to numpy arrays with correct data types
    all_ids = np.array([r['ID'] for r in combined_rows], dtype=np.int64)
    all_waves = np.array([r['WAVE'] for r in combined_rows], dtype=np.float64)
    all_fluxes = np.array([r['FLUX'] for r in combined_rows], dtype=np.float64)
    all_ivars = np.array([r['IVAR'] for r in combined_rows], dtype=np.float64)
    all_redshifts = np.array([r['REDSHIFT'] for r in combined_rows], dtype=np.float64)
    all_ras = np.array([r['RA'] for r in combined_rows], dtype=np.float64)
    all_decs = np.array([r['DEC'] for r in combined_rows], dtype=np.float64)
    all_vdisps = np.array([r['VDISP'] for r in combined_rows], dtype=np.float64)
    all_snrs = np.array([r['SNR'] for r in combined_rows], dtype=np.float64)
    all_survey_types = np.array([r['SURVEY_TYPE'].ljust(15) for r in combined_rows], dtype='S15')

    # Create FITS columns with the numpy arrays
    cols = [
        fits.Column(name='ID', format='K', array=all_ids),
        fits.Column(name='WAVE', format=f'{n_wave}D', array=all_waves),
        fits.Column(name='FLUX', format=f'{n_wave}D', array=all_fluxes),
        fits.Column(name='IVAR', format=f'{n_wave}D', array=all_ivars),
        fits.Column(name='REDSHIFT', format='D', array=all_redshifts),
        fits.Column(name='RA', format='D', array=all_ras),
        fits.Column(name='DEC', format='D', array=all_decs),
        fits.Column(name='VDISP', format='D', array=all_vdisps),
        fits.Column(name='SNR', format='D', array=all_snrs),
        fits.Column(name='SURVEY_TYPE', format='15A', array=all_survey_types)
    ]

    # Create and write final FITS file
    primary_hdu = fits.PrimaryHDU()
    table_hdu = fits.BinTableHDU.from_columns(cols)
    table_hdu.name = 'SPECTRA'
    hdul = fits.HDUList([primary_hdu, table_hdu])
    hdul.writeto(output_file, overwrite=True)
    hdul.close()

    end_time = time.time()
    execution_time = end_time - start_time
    time_str = str(timedelta(seconds=int(execution_time)))
    
    print(f"\n✅ Successfully saved to {output_file}")
    print(f"Total execution time: {time_str}")
    print(f"Processed {len(combined_rows)} galaxies ({len(combined_rows)/execution_time:.1f} galaxies/second)")

if __name__ == "__main__":
    print("Script starting...")
    main()
