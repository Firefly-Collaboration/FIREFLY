"""
===================================================================================
|                         FIREFLY — Full Spectral Fitting                         |
===================================================================================
 <> Module: NERSC_fits_create.py

 <> Author:
	- Samuel Helps ~ <samuel.helps.sh__at__gmail.com>

 <> Contributors:
	- Kieran Graham ~ <kieran.graham__at__port.ac.uk>

 <> Purpose:
	- Creates DESI EDR galaxy FITS files from the data repository on NERSC.
___________________________________________________________________________________
|Institute of Cosmology and Gravitation - University of Portsmouth, Portsmouth, UK|
‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
"""
import os
import numpy as np
from astropy.io import fits
from astropy.table import Table
from tqdm import tqdm
from desitarget.targetmask import desi_mask

# =======================================
#              SETTINGS
# =======================================

# Directory of this script
script_dir = os.path.dirname(os.path.abspath(__file__))
# Path to DESI folder (assumes script is in firefly\Data\DESI\fits_creator)
desi_dir = os.path.abspath(os.path.join(script_dir, ".."))  # one level up
# Output folder
OUTPUT_DIR = os.path.join(desi_dir, "DESI_EDR_data")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Data Sources (on NERSC)
BASE_DIR = "/global/cfs/cdirs/desi/public/edr"
ZCATALOG_PATH = os.path.join(BASE_DIR, "spectro/redux/fuji/zcatalog/zall-pix-fuji.fits")
FASTSPEC_PATH = os.path.join(BASE_DIR, "vac/edr/fastspecfit/fuji/v3.2/catalogs/fastspec-fuji.fits")

# Spectrum range
START_SPECTRUM = 1720000
END_SPECTRUM = 1730000

# Optional Filters
FILTER_BY_REDSHIFT = False
REDSHIFT_MIN = 0.4
REDSHIFT_MAX = 0.8

FILTER_BY_SURVEY = False
ALLOWED_SURVEYS = ['main', 'sv1', 'sv3']  # lowercase

FILTER_BY_PROGRAM = False
ALLOWED_PROGRAMS = ['bright', 'dark']  # lowercase

FILTER_BY_SPECTYPE = True
SPECTYPE_FILTER = 'GALAXY'

FILTER_BY_SUBSURVEY = False
SUBSURVEY = 'LRG'  # Options: 'LRG', 'ELG', 'BGS'
# =======================================

# Load fastspec dictionary (in memory only)
print("Reading fastspec catalog...")
with fits.open(FASTSPEC_PATH) as hdul:
    data = hdul[1].data
    fastspec_dict = {tid: v for tid, v in zip(data['TARGETID'], data['VDISP'])}

def get_velocity_dispersion_from_fastspec(targetid):
    """Retrieve velocity dispersion from fastspec catalog if available"""
    vdisp = fastspec_dict.get(targetid, None)
    return float(vdisp) if vdisp is not None and np.isfinite(vdisp) else None

def get_snr_from_fastspec(row):
    """(FAILSAFE) Calculate SNR using template-based SNR² values from fastspec"""
    try:
        snr2_values = []
        if row['BGS_TARGET'] > 0:
            if np.isfinite(row.get('TSNR2_BGS', np.nan)):
                snr2_values.append(row['TSNR2_BGS'])
        elif row['DESI_TARGET'] & (2**1):  # ELG bit
            if np.isfinite(row.get('TSNR2_ELG', np.nan)):
                snr2_values.append(row['TSNR2_ELG'])
        elif row['DESI_TARGET'] & (2**2):  # LRG bit
            if np.isfinite(row.get('TSNR2_LRG', np.nan)):
                snr2_values.append(row['TSNR2_LRG'])

        if snr2_values:
            snr2 = np.mean(snr2_values)
            return np.sqrt(snr2)  # Convert SNR² to SNR
    except Exception:
        pass
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
    return f"OTHER_{program}".ljust(15)

def calculate_vdisp(wavelength, resolution):
    """(FAILSAFE) Calculate velocity dispersion from resolution matrix"""
    try:
        fwhm = np.median(resolution)
        sigma_lambda = fwhm / 2.355
        lambda_central = np.median(wavelength)
        vdisp = 3e5 * (sigma_lambda / lambda_central)
        return vdisp if np.isfinite(vdisp) and vdisp >= 50 else 125.0
    except:
        return 125.0

def process_galaxy(row):
    """Process individual galaxy row to extract spectrum and metadata"""
    targetid = row['TARGETID']
    ra, dec = row['TARGET_RA'], row['TARGET_DEC']
    z = row['Z']
    pix = row['HEALPIX']
    pixgroup = pix // 100
    program = row['PROGRAM'].strip().lower()
    survey = row['SURVEY'].strip().lower()

    coadd_path = f"{BASE_DIR}/spectro/redux/fuji/healpix/{survey}/{program}/{pixgroup}/{pix}/coadd-{survey}-{program}-{pix}.fits"
    if not os.path.exists(coadd_path):
        print(f"Missing file: {coadd_path}")
        return None

    with fits.open(coadd_path) as hdul:
        fibermap = hdul['FIBERMAP'].data
        match = np.where(fibermap['TARGETID'] == targetid)[0]
        if len(match) == 0:
            return None
        i = match[0]

        b_wave = hdul['B_WAVELENGTH'].data
        r_wave = hdul['R_WAVELENGTH'].data
        z_wave = hdul['Z_WAVELENGTH'].data
        b_flux = hdul['B_FLUX'].data[i]
        r_flux = hdul['R_FLUX'].data[i]
        z_flux = hdul['Z_FLUX'].data[i]
        b_ivar = hdul['B_IVAR'].data[i]
        r_ivar = hdul['R_IVAR'].data[i]
        z_ivar = hdul['Z_IVAR'].data[i]

        vdisp = get_velocity_dispersion_from_fastspec(targetid)
        if vdisp is None:
            b_vdisp = calculate_vdisp(b_wave, hdul['B_RESOLUTION'].data[i][5])
            r_vdisp = calculate_vdisp(r_wave, hdul['R_RESOLUTION'].data[i][5])
            z_vdisp = calculate_vdisp(z_wave, hdul['Z_RESOLUTION'].data[i][5])
            vdisp = np.median([b_vdisp, r_vdisp, z_vdisp])

        wave = np.concatenate([b_wave, r_wave, z_wave])
        flux = np.concatenate([b_flux, r_flux, z_flux])
        ivar = np.concatenate([b_ivar, r_ivar, z_ivar])
        wave_sort = np.argsort(wave)
        wave = wave[wave_sort]
        flux = flux[wave_sort]
        ivar = ivar[wave_sort]

        # First attempt to get SNR from fastspec catalog row info
        snr = get_snr_from_fastspec(row)

        # Fall back to flux and ivar calculation if above failed
        if snr is None:
            valid = (ivar > 0) & np.isfinite(flux) & np.isfinite(ivar)
            if np.any(valid):
                snr = np.median(np.abs(flux[valid] * np.sqrt(ivar[valid])))

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
            'SURVEY_TYPE': get_survey_type(row)
        }

def main():
    # Load full catalog
    zcat = Table.read(ZCATALOG_PATH, hdu=1)
    print(f"Loaded full catalog: {len(zcat)} entries")

    # Apply filters
    if FILTER_BY_SPECTYPE:
        zcat = zcat[zcat['SPECTYPE'] == SPECTYPE_FILTER]
    if FILTER_BY_REDSHIFT:
        zcat = zcat[(zcat['Z'] >= REDSHIFT_MIN) & (zcat['Z'] <= REDSHIFT_MAX)]
    if FILTER_BY_SURVEY:
        zcat = zcat[np.isin(np.char.lower(zcat['SURVEY']), ALLOWED_SURVEYS)]
    if FILTER_BY_PROGRAM:
        zcat = zcat[np.isin(np.char.lower(zcat['PROGRAM']), ALLOWED_PROGRAMS)]
    if FILTER_BY_SUBSURVEY:
        if SUBSURVEY.upper() == 'LRG':
            mask = desi_mask.LRG.mask
        elif SUBSURVEY.upper() == 'ELG':
            mask = desi_mask.ELG.mask
        elif SUBSURVEY.upper() == 'BGS':
            mask = desi_mask.BGS_ANY.mask
        else:
            raise ValueError(f"Unsupported SUBSURVEY: {SUBSURVEY}")
        zcat = zcat[(zcat['DESI_TARGET'] & mask) != 0]

    print(f"Total galaxies matching filters: {len(zcat)}")

    # Select range
    galaxies = zcat[START_SPECTRUM:END_SPECTRUM]
    print(f"Processing {len(galaxies)} galaxies")

    # Process galaxies
    results = []
    for row in tqdm(galaxies):
        r = process_galaxy(row)
        if r:
            results.append(r)

    if not results:
        print("No valid results found.")
        return

    # Create FITS columns
    n_wave = len(results[0]['WAVE'])
    cols = [
        fits.Column(name='ID', format='K', array=np.array([r['ID'] for r in results])),
        fits.Column(name='WAVE', format=f'{n_wave}D', array=np.array([r['WAVE'] for r in results])),
        fits.Column(name='FLUX', format=f'{n_wave}D', array=np.array([r['FLUX'] for r in results])),
        fits.Column(name='IVAR', format=f'{n_wave}D', array=np.array([r['IVAR'] for r in results])),
        fits.Column(name='REDSHIFT', format='D', array=np.array([r['REDSHIFT'] for r in results])),
        fits.Column(name='RA', format='D', array=np.array([r['RA'] for r in results])),
        fits.Column(name='DEC', format='D', array=np.array([r['DEC'] for r in results])),
        fits.Column(name='VDISP', format='D', array=np.array([r['VDISP'] for r in results])),
        fits.Column(name='SNR', format='D', array=np.array([r['SNR'] for r in results])),
        fits.Column(name='SURVEY_TYPE', format='15A', array=np.array([r['SURVEY_TYPE'].encode() for r in results]))
    ]

    # Write to the new Firefly-compatible FITS file
    output_file = os.path.join(OUTPUT_DIR, f"DESI_EDR_{START_SPECTRUM}-{END_SPECTRUM}.fits")
    hdul = fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU.from_columns(cols)])
    hdul.writeto(output_file, overwrite=True)
    print(f"Saved output to {output_file}")

if __name__ == "__main__":
    main()
