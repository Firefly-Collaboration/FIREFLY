"""
===================================================================================
|                         FIREFLY — Full Spectral Fitting                         |
===================================================================================
 <> Module: read_firefly_DESI-SDSS.py

 <> Author:
	- Samuel Helps ~ <samuel.helps.sh__at__gmail.com>

 <> Purpose:
	- Read, analyse and plot FIREFLY outputs for matched DESI and SDSS spectra.

 <> Acknowledgements:
    - Spectra retrieved via SPARCL (sparcl.client). Please cite SPARCL developers
      and data providers when publishing.
___________________________________________________________________________________
|Institute of Cosmology and Gravitation - University of Portsmouth, Portsmouth, UK|
‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
"""
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import fcntl
import os
import sys
from matplotlib.colors import Normalize
from matplotlib.cm import get_cmap
from matplotlib import gridspec

if len(sys.argv) != 5:
    print("Usage: python read_firefly_DESI-SDSS.py <spectrum_index> <file_range> <desi_temp_file> <sdss_temp_file>")
    sys.exit(1)

spec_index = int(sys.argv[1])
file_range = sys.argv[2]  # e.g. "0-116" for greatwall sample
temp_file_desi = sys.argv[3]
temp_file_sdss = sys.argv[4]

# If needed, override path building for temp_files here
hdul_desi = fits.open(temp_file_desi)
hdul_sdss = fits.open(temp_file_sdss)

# Set up output directories
sample = os.environ.get('SAMPLE_NAME', 'sample')
# Derive repo root from script location and place outputs under firefly/Results
script_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.abspath(os.path.join(script_dir, "..", "..", "..", ".."))
base_dir = os.path.join(repo_root, 'Results', f"{sample}_compare_{file_range}")
master_dir = os.path.join(base_dir, 'master_files')
spectra_dir = os.path.join(base_dir, f'spectra_plots_{file_range}')
spectra_comp_dir = os.path.join(base_dir, f'spectra_comparison_{file_range}')

for directory in [master_dir, spectra_dir, spectra_comp_dir]:
    os.makedirs(directory, exist_ok=True)

# -------------------------
# Robust extraction helpers
# -------------------------
def get_column_array(hdu, keys, default=None):
    """Return a numpy array for the first matching key in keys found in hdu.data or None."""
    if hdu is None:
        return default
    data = hdu.data
    if data is None:
        return default
    for k in keys:
        try:
            arr = data[k]
            return np.array(arr)
        except Exception:
            pass
        try:
            names = []
            try:
                names = list(data.columns.names)
            except Exception:
                names = list(data.dtype.names) if hasattr(data, 'dtype') and data.dtype.names else []
            lower_map = {n.lower(): n for n in names}
            if k.lower() in lower_map:
                arr = data[lower_map[k.lower()]]
                return np.array(arr)
        except Exception:
            pass
    return default

def get_header_value(hdul, keys, default=None):
    """Try multiple header locations and keys to find a scalar header value."""
    for h in [hdul[1] if len(hdul) > 1 else None, hdul[0]]:
        if h is None:
            continue
        for k in keys:
            if k in h.header:
                return h.header.get(k)
    return default

# ---------------------------------------
# Extract data from DESI and SDSS outputs
# ---------------------------------------
desi_hdu = hdul_desi[1] if len(hdul_desi) > 1 else hdul_desi[0]
wave_desi = get_column_array(desi_hdu, ['wavelength', 'WAVE', 'wave'], default=None)
flux_desi = get_column_array(desi_hdu, ['original_data', 'original_flux', 'FLUX', 'flux'], default=None)
model_desi = get_column_array(desi_hdu, ['firefly_model', 'model', 'FIREFLY_MODEL'], default=None)
flux_err_desi = get_column_array(desi_hdu, ['flux_error', 'flux_err', 'ERROR', 'error'], default=None)

sdss_hdu = hdul_sdss[1] if len(hdul_sdss) > 1 else hdul_sdss[0]
wave_sdss = get_column_array(sdss_hdu, ['wavelength', 'WAVE', 'wave'], default=None)
flux_sdss = get_column_array(sdss_hdu, ['original_data', 'original_flux', 'FLUX', 'flux'], default=None)
model_sdss = get_column_array(sdss_hdu, ['firefly_model', 'model', 'FIREFLY_MODEL'], default=None)
flux_err_sdss = get_column_array(sdss_hdu, ['flux_error', 'flux_err', 'ERROR', 'error'], default=None)

# Diagnostics
if wave_desi is None:
    print(f"Warning: DESI wavelength column not found in {temp_file_desi}.")
if flux_desi is None:
    print(f"Warning: DESI flux column not found in {temp_file_desi}.")
if wave_sdss is None:
    print(f"Warning: SDSS wavelength column not found in {temp_file_sdss}.")
if flux_sdss is None:
    print(f"Warning: SDSS flux column not found in {temp_file_sdss}.")

if wave_desi is None or flux_desi is None:
    raise KeyError(f"DESI data missing required columns in {temp_file_desi}. Found columns: {getattr(desi_hdu,'columns',None)}")
if wave_sdss is None or flux_sdss is None:
    raise KeyError(f"SDSS data missing required columns in {temp_file_sdss}. Found columns: {getattr(sdss_hdu,'columns',None)}")

# Ensure arrays are 1D numpy arrays
wave_desi = np.asarray(wave_desi).astype(float)
flux_desi = np.asarray(flux_desi).astype(float)
flux_err_desi = np.asarray(flux_err_desi).astype(float) if flux_err_desi is not None else np.full_like(flux_desi, np.inf)
model_desi = np.asarray(model_desi).astype(float) if model_desi is not None else None

wave_sdss = np.asarray(wave_sdss).astype(float)
flux_sdss = np.asarray(flux_sdss).astype(float)
flux_err_sdss = np.asarray(flux_err_sdss).astype(float) if flux_err_sdss is not None else np.full_like(flux_sdss, np.inf)
model_sdss = np.asarray(model_sdss).astype(float) if model_sdss is not None else None

# Retrive SNR, redshift, IDs, metadata etc.
# Try to fetch SNR from the table data first (per-spectrum), then fallback to primary header
snr_desi = hdul_desi[0].header.get('SNR', -1)
try:
    # If the table has an 'SNR' column, read its first element (temp files are single-spectrum)
    if hasattr(desi_hdu, 'data') and desi_hdu.data is not None:
        # Data may be a 1-row table, or arrays; try to access SNR robustly
        if 'SNR' in getattr(desi_hdu.data, 'dtype', {}).names if hasattr(desi_hdu.data, 'dtype') else False:
            snr_arr = desi_hdu.data['SNR']
            snr_desi = float(np.atleast_1d(snr_arr)[0])
        else:
            # try columns.names if available
            try:
                names = list(desi_hdu.columns.names)
                lower_map = {n.lower(): n for n in names}
                if 'snr' in lower_map:
                    snr_arr = desi_hdu.data[lower_map['snr']]
                    snr_desi = float(np.atleast_1d(snr_arr)[0])
            except Exception:
                pass
except Exception:
    pass

snr_sdss = hdul_sdss[0].header.get('SNR', -1)
try:
    if hasattr(sdss_hdu, 'data') and sdss_hdu.data is not None:
        if 'SNR' in getattr(sdss_hdu.data, 'dtype', {}).names if hasattr(sdss_hdu.data, 'dtype') else False:
            snr_arr_s = sdss_hdu.data['SNR']
            snr_sdss = float(np.atleast_1d(snr_arr_s)[0])
        else:
            try:
                names_s = list(sdss_hdu.columns.names)
                lower_map_s = {n.lower(): n for n in names_s}
                if 'snr' in lower_map_s:
                    snr_arr_s = sdss_hdu.data[lower_map_s['snr']]
                    snr_sdss = float(np.atleast_1d(snr_arr_s)[0])
            except Exception:
                pass
except Exception:
    pass

# Redshift detection
try:
    redshift_desi = float(desi_hdu.header.get('redshift', hdul_desi[0].header.get('REDSHIFT', -1)))
except Exception:
    try:
        redshift_desi = float(desi_hdu.data['redshift'][0])
    except Exception:
        redshift_desi = -1.0

try:
    redshift_sdss = float(sdss_hdu.header.get('redshift', hdul_sdss[0].header.get('REDSHIFT', -1)))
except Exception:
    try:
        redshift_sdss = float(sdss_hdu.data['redshift'][0])
    except Exception:
        redshift_sdss = -1.0

redshift_desi = np.round(redshift_desi, 3)
redshift_sdss = np.round(redshift_sdss, 3)

id_desi = hdul_desi[0].header.get('ID', 'N/A')
id_sdss = hdul_sdss[0].header.get('ID', 'N/A')

model_type_desi = desi_hdu.header.get('MODELS', '').strip() if 'MODELS' in desi_hdu.header else ''
model_type_sdss = sdss_hdu.header.get('MODELS', '').strip() if 'MODELS' in sdss_hdu.header else ''
full_model_name_desi = model_type_desi if model_type_desi else 'MaStar'
full_model_name_sdss = model_type_sdss if model_type_sdss else 'MaStar'

# Title strings: SDSS title keeps only z in brackets
title_desi = f"{sample} DESI Galaxy #{spec_index} (ID: {id_desi}, z={redshift_desi}, S/N={snr_desi:.1f})"
title_sdss = f"{sample} SDSS Galaxy #{spec_index} (z={redshift_sdss})"

# ------------------------------------------------------
# Function to extract parameters from a FIREFLY FITS HDU
# ------------------------------------------------------
def extract_parameters(hdul):
    header = hdul[1].header
    params = {
        'spec_index': header.get('SPEC_IDX', -1),
        'id': hdul[0].header.get('ID', -1),
        'ra': hdul[0].header.get('RA', -1),
        'dec': hdul[0].header.get('DEC', -1),
        'redshift': header.get('redshift', -1),
        'snr': hdul[0].header.get('SNR', -1),
        'survey_type': hdul[0].header.get('SURVEY_TYPE', 'Unknown'),
        'model_type': header.get('MODEL', 'Unknown'),
        'imf': header.get('IMF', 'Unknown'),
        'converged': header.get('converged', False),
        'age_lightW': header.get('age_lightW', -1),
        'age_lightW_up_1sig': header.get('age_lightW_up_1sig', -1),
        'age_lightW_low_1sig': header.get('age_lightW_low_1sig', -1),
        'metallicity_lightW': header.get('metallicity_lightW', -1),
        'metallicity_lightW_up_1sig': header.get('metallicity_lightW_up_1sig', -1),
        'metallicity_lightW_low_1sig': header.get('metallicity_lightW_low_1sig', -1),
        'age_massW': header.get('age_massW', -1),
        'age_massW_up_1sig': header.get('age_massW_up_1sig', -1),
        'age_massW_low_1sig': header.get('age_massW_low_1sig', -1),
        'metallicity_massW': header.get('metallicity_massW', -1),
        'metallicity_massW_up_1sig': header.get('metallicity_massW_up_1sig', -1),
        'metallicity_massW_low_1sig': header.get('metallicity_massW_low_1sig', -1),
        'total_mass': header.get('total_mass', -1),
        'stellar_mass': header.get('stellar_mass', -1),
        'EBV': header.get('EBV', -1),
        'ssp_number': header.get('ssp_number', 0)
    }
    ssp_number = int(params['ssp_number']) if params['ssp_number'] is not None else 0
    for i in range(ssp_number):
        params[f'total_mass_ssp_{i}'] = header.get(f'total_mass_ssp_{i}', -1)
        params[f'stellar_mass_ssp_{i}'] = header.get(f'stellar_mass_ssp_{i}', -1)
        params[f'log_age_ssp_{i}'] = header.get(f'log_age_ssp_{i}', -1)
        params[f'metal_ssp_{i}'] = header.get(f'metal_ssp_{i}', -1)
        params[f'SFR_ssp_{i}'] = header.get(f'SFR_ssp_{i}', -1)
        params[f'weightMass_ssp_{i}'] = header.get(f'weightMass_ssp_{i}', -1)
        params[f'weightLight_ssp_{i}'] = header.get(f'weightLight_ssp_{i}', -1)
    return params

params_desi = extract_parameters(hdul_desi)
params_sdss = extract_parameters(hdul_sdss)

# ------------------
# Update master CSVs
# ------------------
cut_tag = 'Cut' if os.environ.get('CUT_MODE', '').lower() == 'cut' else 'NoCut'
master_file_desi = os.path.join(master_dir, f"{sample}_DESI_master_{cut_tag}.csv")
master_file_sdss = os.path.join(master_dir, f"{sample}_SDSS_master_{cut_tag}.csv")

def update_master_file(params, master_file):
    os.makedirs(os.path.dirname(master_file), exist_ok=True)
    column_names = {
        'spec_index': 'Spectrum_Index',
        'id': 'ID',
        'ra': 'RA',
        'dec': 'DEC',
        'redshift': 'Redshift',
        'snr': 'Signal_to_Noise_Ratio',
        'survey_type': 'Survey_Type',
        'model_type': 'Model_Type',
        'converged': 'Converged',
        'age_lightW': 'Light_Weighted_Age_log_Gyr',
        'metallicity_lightW': 'Light_Weighted_Z_H',
        'stellar_mass': 'log_M_Msun',
        'EBV': 'E_B_V',
    }
    with open(master_file + '.lock', 'w') as lockfile:
        try:
            fcntl.flock(lockfile, fcntl.LOCK_EX)
            try:
                df = pd.read_csv(master_file)
                reverse_map = {v: k for k, v in column_names.items()}
                df = df.rename(columns=reverse_map)
                new_row = pd.DataFrame([params])
                if params['spec_index'] in df['spec_index'].values:
                    idx = df.index[df['spec_index'] == params['spec_index']].item()
                    df.loc[idx, new_row.columns] = new_row.iloc[0]
                else:
                    df = pd.concat([df, new_row], ignore_index=True)
            except (FileNotFoundError, pd.errors.EmptyDataError):
                df = pd.DataFrame([params])
            df = df.sort_values('spec_index')
            df = df.rename(columns=column_names)
            df.to_csv(master_file, index=False)
        finally:
            fcntl.flock(lockfile, fcntl.LOCK_UN)

update_master_file(params_desi, master_file_desi)
update_master_file(params_sdss, master_file_sdss)

# -----------------------------
# Helper for safe header values
# -----------------------------
def safe_header_val(hdu, key, fallback=-1):
    try:
        return hdu[1].header.get(key, fallback)
    except Exception:
        return fallback

# ---------------------------------
# Create combined comparison figure
# ---------------------------------
fig, axes = plt.subplots(3, 2, figsize=(16, 10))
plt.subplots_adjust(hspace=0.3)

# DESI spectrum plot
axes[0,0].plot(wave_desi, flux_desi, linewidth=0.6, label='DESI data', alpha=0.8)
if model_desi is not None:
    axes[0,0].plot(wave_desi, model_desi, linewidth=1.2, label=f'Model: {full_model_name_desi}')
axes[0,0].set_ylabel('Flux [$10^{-17}$ erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$]')
axes[0,0].set_title(title_desi)
axes[0,0].legend(loc='best', fontsize='small')

# DESI params text
age_lightW_desi = safe_header_val(hdul_desi, 'age_lightW', -99)
metallicity_lightW_desi = safe_header_val(hdul_desi, 'metallicity_lightW', -99)
stellar_mass_desi = safe_header_val(hdul_desi, 'stellar_mass', -99)
ebv_desi = safe_header_val(hdul_desi, 'EBV', -99)
params_text_desi = (f'Age: {10**age_lightW_desi:.2f} Gyr\n'
                    f'[Z/H]: {metallicity_lightW_desi:.2f} dex\n'
                    f'log M/M☉: {stellar_mass_desi:.2f}\n'
                    f'E(B-V): {ebv_desi:.2f}')
axes[0,0].annotate(params_text_desi, xy=(0.02, 0.95), xycoords='axes fraction',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8), verticalalignment='top')

# SDSS spectrum plot
axes[0,1].plot(wave_sdss, flux_sdss, linewidth=0.6, label='SDSS data', alpha=0.8)
if model_sdss is not None:
    axes[0,1].plot(wave_sdss, model_sdss, linewidth=1.2, label=f'Model: {full_model_name_sdss}')
axes[0,1].set_title(title_sdss)
axes[0,1].legend(loc='best', fontsize='small')

age_lightW_sdss = safe_header_val(hdul_sdss, 'age_lightW', -99)
metallicity_lightW_sdss = safe_header_val(hdul_sdss, 'metallicity_lightW', -99)
stellar_mass_sdss = safe_header_val(hdul_sdss, 'stellar_mass', -99)
ebv_sdss = safe_header_val(hdul_sdss, 'EBV', -99)
params_text_sdss = (f'Age: {10**age_lightW_sdss:.2f} Gyr\n'
                    f'[Z/H]: {metallicity_lightW_sdss:.2f} dex\n'
                    f'log M/M☉: {stellar_mass_sdss:.2f}\n'
                    f'E(B-V): {ebv_sdss:.2f}')
axes[0,1].annotate(params_text_sdss, xy=(0.02, 0.95), xycoords='axes fraction',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8), verticalalignment='top')

# Lookback time & metallicity
ssp_n_desi = int(params_desi.get('ssp_number', 0) or 0)
ssp_n_sdss = int(params_sdss.get('ssp_number', 0) or 0)

if ssp_n_desi > 0:
    csp_age_desi = np.array([hdul_desi[1].header.get(f'log_age_ssp_{k}', -99) for k in range(ssp_n_desi)], dtype=float)
    csp_light_desi = np.array([hdul_desi[1].header.get(f'weightLight_ssp_{k}', 0.0) for k in range(ssp_n_desi)], dtype=float)
else:
    csp_age_desi = np.array([0.0])
    csp_light_desi = np.array([0.0])

if ssp_n_sdss > 0:
    csp_age_sdss = np.array([hdul_sdss[1].header.get(f'log_age_ssp_{k}', -99) for k in range(ssp_n_sdss)], dtype=float)
    csp_light_sdss = np.array([hdul_sdss[1].header.get(f'weightLight_ssp_{k}', 0.0) for k in range(ssp_n_sdss)], dtype=float)
else:
    csp_age_sdss = np.array([0.0])
    csp_light_sdss = np.array([0.0])

axes[1,0].bar(10**(csp_age_desi), csp_light_desi, width=1.0, alpha=0.5, color='#779ECB')
axes[1,0].set_xlabel('Lookback time (Gyr)')
axes[1,0].set_ylabel('Frequency')
axes[1,0].set_title('DESI Lookback Time')
axes[1,0].set_xlim(0, 14)
axes[1,0].set_ylim(0, max(csp_light_desi.max(), 1) * 1.1)
axes[1,0].annotate(params_text_desi, xy=(0.02, 0.95), xycoords='axes fraction',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8), verticalalignment='top')

axes[1,1].bar(10**(csp_age_sdss), csp_light_sdss, width=1.0, alpha=0.5, color='#779ECB')
axes[1,1].set_xlabel('Lookback time (Gyr)')
axes[1,1].set_title('SDSS Lookback Time')
axes[1,1].set_xlim(0, 14)
axes[1,1].set_ylim(0, max(csp_light_sdss.max(), 1) * 1.1)
axes[1,1].annotate(params_text_sdss, xy=(0.02, 0.95), xycoords='axes fraction',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8), verticalalignment='top')

if ssp_n_desi > 0:
    csp_Z_desi = np.array([hdul_desi[1].header.get(f'metal_ssp_{k}', 0.0) for k in range(ssp_n_desi)], dtype=float)
else:
    csp_Z_desi = np.array([0.0])

if ssp_n_sdss > 0:
    csp_Z_sdss = np.array([hdul_sdss[1].header.get(f'metal_ssp_{k}', 0.0) for k in range(ssp_n_sdss)], dtype=float)
else:
    csp_Z_sdss = np.array([0.0])

axes[2,0].bar(csp_Z_desi, csp_light_desi, width=0.1, alpha=0.5, color='#B0171F')
axes[2,0].set_xlabel('[Z/H] (dex)')
axes[2,0].set_ylabel('Frequency')
axes[2,0].set_title('DESI Metallicity')
axes[2,0].set_xlim(-2, 2)
axes[2,0].set_ylim(0, max(csp_light_desi.max(), 1) * 1.1)
axes[2,0].annotate(params_text_desi, xy=(0.02, 0.95), xycoords='axes fraction',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8), verticalalignment='top')

axes[2,1].bar(csp_Z_sdss, csp_light_sdss, width=0.1, alpha=0.5, color='#B0171F')
axes[2,1].set_xlabel('[Z/H] (dex)')
axes[2,1].set_title('SDSS Metallicity')
axes[2,1].set_xlim(-2, 2)
axes[2,1].set_ylim(0, max(csp_light_sdss.max(), 1) * 1.1)
axes[2,1].annotate(params_text_sdss, xy=(0.02, 0.95), xycoords='axes fraction',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8), verticalalignment='top')

# Save combined figure
fig_path = os.path.join(spectra_dir, f'{spec_index}_comparison_{file_range}.png')
plt.savefig(fig_path, bbox_inches='tight', dpi=300)
plt.close(fig)

# ------------------------------------------------------
# Create spectra-only comparison (DESI left, SDSS right)
# ------------------------------------------------------
fig2 = plt.figure(figsize=(20, 4.7))
gs = gridspec.GridSpec(3, 2, height_ratios=[3, 1, 0.6], hspace=0.0, wspace=0.3)

# Left column (DESI)
ax1 = fig2.add_subplot(gs[0, 0])  # top: spectrum + model (height ratio 3)
ax2 = fig2.add_subplot(gs[1, 0], sharex=ax1)  # middle: residuals (height ratio 1)
ax3 = fig2.add_subplot(gs[2, 0], sharex=ax1)  # bottom: colour bar (height ratio 0.6)

# Right column (SDSS)
ax1_r = fig2.add_subplot(gs[0, 1])  # top
ax2_r = fig2.add_subplot(gs[1, 1], sharex=ax1_r)  # middle
ax3_r = fig2.add_subplot(gs[2, 1], sharex=ax1_r)  # bottom

# DESI top panel (spectrum + model)
ax1.plot(wave_desi, flux_desi, color='#006699', linewidth=0.6, label='Original Data', alpha=0.8)
if model_desi is not None:
    ax1.plot(wave_desi, model_desi, color='#DA9C00', linewidth=1.4, label=f'Model: {full_model_name_desi}')
ax1.set_ylabel('Flux [$10^{-17}$ erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$]')
ax1.set_title(title_desi)
ax1.legend(loc='upper right', bbox_to_anchor=(0.98, 0.98))
# build params_text_desi if not already built:
age_lightW_desi = safe_header_val(hdul_desi, 'age_lightW', -99)
metallicity_lightW_desi = safe_header_val(hdul_desi, 'metallicity_lightW', -99)
stellar_mass_desi = safe_header_val(hdul_desi, 'stellar_mass', -99)
ebv_desi = safe_header_val(hdul_desi, 'EBV', -99)
params_text_desi = (f'age: {np.around(10**age_lightW_desi,decimals=2)} Gyr\n'
                    f'[Z/H]: {np.around(metallicity_lightW_desi,decimals=2)} dex\n'
                    f'log M/Msun: {np.around(stellar_mass_desi,decimals=2)}\n'
                    f'E(B-V): {np.around(ebv_desi,decimals=2)} mag')
ax1.annotate(params_text_desi, xy=(0.02, 0.95), xycoords='axes fraction',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
             verticalalignment='top')

# DESI middle panel (residuals)
if model_desi is not None:
    residuals_desi = flux_desi - model_desi
else:
    residuals_desi = flux_desi * 0.0
ax2.plot(wave_desi, residuals_desi, 'k-', linewidth=0.6, label='Residual')
ax2.axhline(y=0, color='r', linestyle='--', alpha=0.5)
ax2.legend()

# DESI bottom panel (flux-colour bar)
flux_norm = (flux_desi - flux_desi.min()) / (flux_desi.max() - flux_desi.min() + 1e-30)
cmap = get_cmap('rainbow')
wave_norm = Normalize(vmin=wave_desi.min(), vmax=wave_desi.max())
colored_spectrum = cmap(wave_norm(wave_desi))[:, :3]
height = 50
y = np.linspace(0, 1, height + 1)
x = np.append(wave_desi, wave_desi[-1])
flux_grid = np.tile(flux_norm, (height, 1))
rgb_array = np.zeros((height, len(wave_desi), 3))
for c in range(3):
    rgb_array[:, :, c] = flux_grid * colored_spectrum[:, c]
ax3.set_yticks([])
ax3.set_xlabel('Wavelength [Å]')
# Draw pcolormesh onto ax3
ax3.pcolormesh(x, y, flux_grid, color=rgb_array.reshape(-1, 3), shading='auto')
ax3.set_xlim(wave_desi.min(), wave_desi.max())
ax3.set_ylim(0, 1)

# SDSS top panel
ax1_r.plot(wave_sdss, flux_sdss, color='#006699', linewidth=0.6, label='Original Data', alpha=0.8)
if model_sdss is not None:
    ax1_r.plot(wave_sdss, model_sdss, color='#DA9C00', linewidth=1.4, label=f'Model: {full_model_name_sdss}')
ax1_r.set_title(title_sdss)
ax1_r.legend(loc='upper right', bbox_to_anchor=(0.98, 0.98))
age_lightW_sdss = safe_header_val(hdul_sdss, 'age_lightW', -99)
metallicity_lightW_sdss = safe_header_val(hdul_sdss, 'metallicity_lightW', -99)
stellar_mass_sdss = safe_header_val(hdul_sdss, 'stellar_mass', -99)
ebv_sdss = safe_header_val(hdul_sdss, 'EBV', -99)
params_text_sdss = (f'age: {np.around(10**age_lightW_sdss,decimals=2)} Gyr\n'
                    f'[Z/H]: {np.around(metallicity_lightW_sdss,decimals=2)} dex\n'
                    f'log M/Msun: {np.around(stellar_mass_sdss,decimals=2)}\n'
                    f'E(B-V): {np.around(ebv_sdss,decimals=2)} mag')
ax1_r.annotate(params_text_sdss, xy=(0.02, 0.95), xycoords='axes fraction',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
               verticalalignment='top')

# SDSS middle panel (residuals)
if model_sdss is not None:
    residuals_sdss = flux_sdss - model_sdss
else:
    residuals_sdss = flux_sdss * 0.0
ax2_r.plot(wave_sdss, residuals_sdss, 'k-', linewidth=0.6, label='Residual')
ax2_r.axhline(y=0, color='r', linestyle='--', alpha=0.5)
ax2_r.legend()

# SDSS bottom panel (colour bar)
flux_norm_s = (flux_sdss - flux_sdss.min()) / (flux_sdss.max() - flux_sdss.min() + 1e-30)
cmap_s = get_cmap('rainbow')
wave_norm_s = Normalize(vmin=wave_sdss.min(), vmax=wave_sdss.max())
colored_spectrum_s = cmap_s(wave_norm_s(wave_sdss))[:, :3]
flux_grid_s = np.tile(flux_norm_s, (height, 1))
rgb_array_s = np.zeros((height, len(wave_sdss), 3))
for c in range(3):
    rgb_array_s[:, :, c] = flux_grid_s * colored_spectrum_s[:, c]
ax3_r.set_yticks([])
ax3_r.set_xlabel('Wavelength [Å]')
ax3_r.pcolormesh(np.append(wave_sdss, wave_sdss[-1]), y, flux_grid_s, color=rgb_array_s.reshape(-1, 3), shading='auto')
ax3_r.set_xlim(wave_sdss.min(), wave_sdss.max())
ax3_r.set_ylim(0, 1)

# match layout exactly: remove spacing between stacked panels
fig2.subplots_adjust(hspace=0)

# Save spectra-only figure
spectra_only_path = os.path.join(spectra_comp_dir, f'{spec_index}_spectra_only_{file_range}.png')
plt.savefig(spectra_only_path, bbox_inches='tight', dpi=300)
plt.close(fig2)

print(f"\nProcessed spectrum #{spec_index}. Combined plot saved to {fig_path}. Spectra-only saved to {spectra_only_path}.")
