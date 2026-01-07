"""
===================================================================================
|                         FIREFLY — Full Spectral Fitting                         |
===================================================================================
 <> Module: read_firefly_DESI_EDR.py

 <> Author:
	- Samuel Helps ~ <samuel.helps.sh__at__gmail.com>

 <> Purpose:
	- Read, analyse and plot FIREFLY outputs for DESI EDR spectra into per file
      result directories with a master CSV, fitted spectra, lookback time and
      metallicity plots.
___________________________________________________________________________________
|Institute of Cosmology and Gravitation - University of Portsmouth, Portsmouth, UK|
‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
"""
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import get_cmap
import numpy as np
import pandas as pd
import fcntl
import os
import sys


if len(sys.argv) != 4:
    print("Usage: python read_firefly_DESI_EDR.py <spectrum_index> <file_range> <temp_file>")
    sys.exit(1)

spec_index = int(sys.argv[1])
file_range = sys.argv[2]  # e.g., "0-10000"
temp_file = sys.argv[3]   # Full path to temp file

hdul = fits.open(temp_file)

# Define output directories
##################### ONLY CHANGE GALAXY INDEX ##############################
start_index, end_index = map(int, file_range.split('-'))
#############################################################################

script_dir = os.path.dirname(os.path.abspath(__file__))  # .../Firefly/Launch/DESI/SCIAMA/core
repo_root = os.path.abspath(os.path.join(script_dir, "..", "..", "..", ".."))  # <repo>/Firefly

remote_base_dir = os.path.join(repo_root, "Results", f"DESI_EDR_{file_range}")

master_file = os.path.join(remote_base_dir, f'master_files/EDR_master_{file_range}.csv')
spectrum_output_dir = os.path.join(remote_base_dir, f'spectra_plots_{file_range}')
lookback_time_output_dir = os.path.join(remote_base_dir, f'lookback_time_plots_{file_range}')
metallicity_output_dir = os.path.join(remote_base_dir, f'metallicity_plots_{file_range}')

# Create all required directories
for directory in [
    os.path.dirname(master_file),
    spectrum_output_dir,
    lookback_time_output_dir,
    metallicity_output_dir
]:
    os.makedirs(directory, exist_ok=True)

# Load data from temp FITS file
data=hdul[1].data
wave = data['wavelength']
flux = data['original_data']
model = data['firefly_model']
flux_error = data['flux_error'] 
snr = hdul[0].header['SNR']
redshift = np.round(hdul[1].header['redshift'], 3)
id = hdul[0].header.get('ID', -1)
classification = "Galaxy"

# Get model type and library with error handling
try:
    model_type = hdul[1].header['MODELS'].strip()  # Try to get MODELS first
except KeyError:
    try:
        model_type = 'M11' if 'MILES' in hdul[1].header['MODEL'] else 'MaStar'  # Determine from MODEL keyword
    except KeyError:
        model_type = 'Unknown'  # Default if neither keyword exists
try:
    model_lib = hdul[1].header['MODEL'].strip()
except KeyError:
    model_lib = 'Unknown'

# Create descriptive model name based on model type
if model_type == 'M11':
    if model_lib in ['MILES', 'STELIB', 'ELODIE', 'MARCS']:
        full_model_name = f"M11-{model_lib}"
    else:
        full_model_name = "M11-Unknown"
elif model_type == 'MaStar':
    full_model_name = "MaStar"
else:
    full_model_name = "Unknown model"
try:
    spec_index = hdul[1].header['SPEC_IDX']  # Get the stored index
except KeyError:
    spec_index = "Unknown"  # Fallback if index not found

# Universal title string for all plots
title_string = f"EDR {classification} #{spec_index} (ID: {id}, z = {redshift}, S/N = {snr:.2f})"

# Universal derived stellar parameters text box for all plots
params_text = (f'age: {np.around(10**hdul[1].header["age_lightW"],decimals=2)} Gyr\n'
              f'[Z/H]: {np.around(hdul[1].header["metallicity_lightW"],decimals=2)} dex\n'
              f'log M/Msun: {np.around(hdul[1].header["stellar_mass"],decimals=2)}\n'
              f'E(B-V): {np.around(hdul[1].header["EBV"],decimals=2)} mag')

csp_age=np.ndarray(hdul[1].header['ssp_number'])
csp_Z=np.ndarray(hdul[1].header['ssp_number'])
csp_light=np.ndarray(hdul[1].header['ssp_number'])
csp_mass=np.ndarray(hdul[1].header['ssp_number'])
for i in range(len(csp_age)):
	csp_age[i]=hdul[1].header['log_age_ssp_'+str(i)]
	csp_Z[i]=hdul[1].header['metal_ssp_'+str(i)]
	csp_light[i]=hdul[1].header['weightLight_ssp_'+str(i)]
	csp_mass[i]=hdul[1].header['weightMass_ssp_'+str(i)]

############################
# Spectra and model fit plot
############################
fig1, (ax1, ax2, ax3) = plt.subplots(3, 1, height_ratios=[3, 1, 0.6], sharex=True, figsize=(10, 4.7))
fig1.subplots_adjust(hspace=0)
x_min, x_max = wave.min(), wave.max()
ax1.plot(wave, flux, color='#006699', linewidth=0.6, label='Original Data', alpha=0.8)
if model is not None:
    ax1.plot(wave, model, color='#DA9C00', linewidth=1.4, label=f'Model: {full_model_name}')
#ax1.plot(wave, 1/np.sqrt(flux_error), 'k-', alpha=0.5, label='Flux Error')
ax1.set_ylabel('Flux [$10^{-17}$ erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$]')
ax1.set_title(title_string)
ax1.legend(loc='upper right', bbox_to_anchor=(0.98, 0.98))
ax1.annotate(params_text, xy=(0.02, 0.95), xycoords='axes fraction',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
             verticalalignment='top')
# Plot residuals in middle panel
residuals = flux - model
ax2.plot(wave, residuals, 'k-', linewidth=0.6, label='Residual')
ax2.axhline(y=0, color='r', linestyle='--', alpha=0.5)
#ax2.set_ylabel('Residuals')
ax2.legend()

# Normalise flux for brightness scaling on spectrum responsive colour bar
flux_norm = (flux - flux.min()) / (flux.max() - flux.min())  # Scale between 0 and 1

# Define the color mapping based on wavelength
cmap = get_cmap('rainbow')  # Use a full spectrum colormap
wave_norm = Normalize(vmin=wave.min(), vmax=wave.max())  # Normalizs wavelength for color mapping
colored_spectrum = cmap(wave_norm(wave))[:, :3]  # Extract RGB values (ignore alpha)

# Create a 2D grid for pcolormesh
height = 50  # Number of pixels in height for the spectrum image
y = np.linspace(0, 1, height + 1)  # Y grid for pcolormesh
x = np.append(wave, wave[-1])  # Extend X grid to match edges

# Create the flux grid
flux_grid = np.tile(flux_norm, (height, 1))  # Tile flux for a full image

# Convert color data into a format usable by pcolormesh
rgb_array = np.zeros((height, len(wave), 3))  # Initialize RGB array
for i in range(3):  # Iterate over R, G, B channels
    rgb_array[:, :, i] = flux_grid * colored_spectrum[:, i]  # Apply color without distorting brightness

# Use pcolormesh to plot the image
c = ax3.pcolormesh(x, y, flux_grid, 
                   color=rgb_array.reshape(-1, 3),  # Flatten RGB for correct mapping
                   shading='auto')

# Set the axis properties
ax3.set_xlim(wave.min(), wave.max())
ax3.set_ylim(0, 1)
ax3.set_yticks([])
ax3.set_xlabel('Wavelength [Å]')

fig1.savefig(f'{spectrum_output_dir}/{spec_index}_spectrum_(EDR_{start_index}-{end_index}).png', 
             bbox_inches='tight', dpi=300)
plt.close(fig1)
'''
# Add emission lines reference
emission_lines = {
    'He-II': [3202.15, 4685.74],
    'Ne-V': [3345.81, 3425.81],
    'O-II': [3726.03, 3728.73],
    'Ne-III': [3868.69, 3967.40],
    'H-ζ': 3889.05,
    'H-ε': 3970.07,
    'H-δ': 4101.73,
    'H-γ': 4340.46,
    'O-III': [4363.15, 4958.83, 5006.77],
    'Ar-IV': [4711.30, 4740.10],
    'H-β': 4861.32,
    'N-I': [5197.90, 5200.39],
    'He-I': 5875.60,
    'O-I': [6300.20, 6363.67],
    'N-II': [6547.96, 6583.34],
    'H-α': 6562.80,
    'S-II': [6716.31, 6730.68],
    'Ar-III': 7135.67
}

# Add emission line markers to the main plot and spectrum image
for line, wavelengths in emission_lines.items():
    if isinstance(wavelengths, list):
        for wav in wavelengths:
            ax1.axvline(x=wav, color='blue', linestyle='--', alpha=0.5)
            ax2.axvline(x=wav, color='blue', linestyle='--', alpha=0.5)
            ax3.axvline(x=wav, color='blue', linestyle='--', alpha=0.5)
    else:
        ax1.axvline(x=wavelengths, color='blue', linestyle='--', alpha=0.5)
        ax2.axvline(x=wavelengths, color='blue', linestyle='--', alpha=0.5)
        ax3.axvline(x=wavelengths, color='blue', linestyle='--', alpha=0.5)
'''
##########################
# lookback time (Gyr) plot
# ########################
fig2=plt.figure()
plt.xlim(0,14)
plt.xlabel('lookback time (Gyr)')
plt.ylabel('frequency')
plt.title(title_string)
plt.bar(10**(csp_age),csp_light,width=1,align='center',alpha=0.5)
plt.scatter(10**(csp_age),csp_light)
plt.annotate(params_text, xy=(0.02, 0.95), xycoords='axes fraction',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
            verticalalignment='top')
plt.savefig(f'{lookback_time_output_dir}/{spec_index}_lookback_time_(EDR_{start_index}-{end_index}).png', 
           bbox_inches='tight', dpi=300)
plt.close(fig2)

##################################
# [Z/H] Metallicity frequency plot
##################################
fig3=plt.figure()
plt.xlim(-2,2)
plt.xlabel('[Z/H] (dex)')
plt.ylabel('frequency')
plt.title(title_string)
plt.bar(csp_Z,csp_light,width=0.1,align='center',alpha=0.5)
plt.scatter(csp_Z,csp_light)
plt.annotate(params_text, xy=(0.02, 0.95), xycoords='axes fraction',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
            verticalalignment='top')
plt.savefig(f'{metallicity_output_dir}/{spec_index}_metallicity_(EDR_{start_index}-{end_index}).png', 
           bbox_inches='tight', dpi=300)
plt.close(fig3)

##############################
# Write/update master CSV file
##############################
def extract_parameters(hdul):
    """Extract all relevant parameters from a FIREFLY output file."""
    params = {
        # Basic parameters
        'spec_index': hdul[1].header.get('SPEC_IDX', -1),
        'id': hdul[0].header.get('ID', -1),
        'ra': hdul[0].header.get('RA', -1),
        'dec': hdul[0].header.get('DEC', -1),
        'redshift': hdul[1].header.get('redshift', -1),
        'snr': hdul[0].header.get('SNR', -1),
        'survey_type': hdul[0].header.get('SURVEY_TYPE', 'Unknown'),
        'model_type': hdul[1].header.get('MODEL', 'Unknown'),
        'imf': hdul[1].header.get('IMF', 'Unknown'),
        'converged': hdul[1].header.get('converged', False),
        
        # Light-weighted parameters
        'age_lightW': hdul[1].header.get('age_lightW', -1),
        'age_lightW_up_1sig': hdul[1].header.get('age_lightW_up_1sig', -1),
        'age_lightW_low_1sig': hdul[1].header.get('age_lightW_low_1sig', -1),
        'metallicity_lightW': hdul[1].header.get('metallicity_lightW', -1),
        'metallicity_lightW_up_1sig': hdul[1].header.get('metallicity_lightW_up_1sig', -1),
        'metallicity_lightW_low_1sig': hdul[1].header.get('metallicity_lightW_low_1sig', -1),
        
        # Mass-weighted parameters
        'age_massW': hdul[1].header.get('age_massW', -1),
        'age_massW_up_1sig': hdul[1].header.get('age_massW_up_1sig', -1),
        'age_massW_low_1sig': hdul[1].header.get('age_massW_low_1sig', -1),
        'metallicity_massW': hdul[1].header.get('metallicity_massW', -1),
        'metallicity_massW_up_1sig': hdul[1].header.get('metallicity_massW_up_1sig', -1),
        'metallicity_massW_low_1sig': hdul[1].header.get('metallicity_massW_low_1sig', -1),
        
        # Mass parameters
        'total_mass': hdul[1].header.get('total_mass', -1),
        'stellar_mass': hdul[1].header.get('stellar_mass', -1),
        'living_stars_mass': hdul[1].header.get('living_stars_mass', -1),
        'remnant_mass': hdul[1].header.get('remnant_mass', -1),
        'remnant_mass_in_whitedwarfs': hdul[1].header.get('remnant_mass_in_whitedwarfs', -1),
        'remnant_mass_in_neutronstars': hdul[1].header.get('remnant_mass_in_neutronstars', -1),
        'remnant_mass_blackholes': hdul[1].header.get('remnant_mass_blackholes', -1),
        'mass_of_ejecta': hdul[1].header.get('mass_of_ejecta', -1),
        
        # Dust parameter
        'EBV': hdul[1].header.get('EBV', -1),
        
        # SSP parameters
        'ssp_number': hdul[1].header.get('ssp_number', 0)
    }
    
    # Add SSP-specific parameters
    ssp_number = params['ssp_number']
    for i in range(ssp_number):
        ssp_params = {
            f'total_mass_ssp_{i}': hdul[1].header.get(f'total_mass_ssp_{i}', -1),
            f'stellar_mass_ssp_{i}': hdul[1].header.get(f'stellar_mass_ssp_{i}', -1),
            f'log_age_ssp_{i}': hdul[1].header.get(f'log_age_ssp_{i}', -1),
            f'metal_ssp_{i}': hdul[1].header.get(f'metal_ssp_{i}', -1),
            f'SFR_ssp_{i}': hdul[1].header.get(f'SFR_ssp_{i}', -1),
            f'weightMass_ssp_{i}': hdul[1].header.get(f'weightMass_ssp_{i}', -1),
            f'weightLight_ssp_{i}': hdul[1].header.get(f'weightLight_ssp_{i}', -1)
        }
        params.update(ssp_params)
    return params

def update_master_file(params, master_file):
    """Update or create master CSV file with FIREFLY results."""
    # Create directory if it doesn't exist
    os.makedirs(os.path.dirname(master_file), exist_ok=True)

    # Define column names for better readability
    column_names = {
        'spec_index': 'Spectrum_Index',
        'id': 'ID',
        'ra': 'RA',
        'dec': 'DEC',
        'redshift': 'Redshift',
        'snr': 'Signal_to_Noise_Ratio',
        'survey_type': 'Survey_Type',
        'model_type': 'Model_Type',
        'imf': 'Initial_Mass_Function',
        'converged': 'Converged',
        'age_lightW': 'Light_Weighted_Age_log_Gyr',
        'metallicity_lightW': 'Light_Weighted_Z_H',
        'stellar_mass': 'log_M_Msun',
        'EBV': 'E_B_V',
    }

    # Use file locking for thread safety
    with open(master_file + '.lock', 'w') as lockfile:
        try:
            # Acquire exclusive lock
            fcntl.flock(lockfile, fcntl.LOCK_EX)
            
            try:
                # Try to read existing file
                df = pd.read_csv(master_file)
                reverse_mapping = {v: k for k, v in column_names.items()}
                df = df.rename(columns=reverse_mapping)
                
                new_row = pd.DataFrame([params])
                if params['spec_index'] in df['spec_index'].values:
                    idx = df.index[df['spec_index'] == params['spec_index']].item()
                    df.loc[idx, new_row.columns] = new_row.iloc[0]
                else:
                    df = pd.concat([df, new_row], ignore_index=True)
            except (FileNotFoundError, pd.errors.EmptyDataError):
                df = pd.DataFrame([params])
            
            # Sort and save
            df = df.sort_values('spec_index')
            df = df.rename(columns=column_names)
            df.to_csv(master_file, index=False)
            
        finally:
            # Release lock
            fcntl.flock(lockfile, fcntl.LOCK_UN)

# After loading and processing the FITS file
params = extract_parameters(hdul)
df = update_master_file(params, master_file)

# Print summary of the processed data
print("\nProcessed parameters for spectrum #", params['spec_index'])
print(f"Redshift: {params['redshift']:.4f}")
print(f"SNR: {params['snr']:.2f}")
print(f"Age (light-weighted): {10**params['age_lightW']:.2f} Gyr")
print(f"[Z/H] (light-weighted): {params['metallicity_lightW']:.2f}")
print(f"log M/Msun: {params['stellar_mass']:.2f}")
print(f"E(B-V): {params['EBV']:.2f}")