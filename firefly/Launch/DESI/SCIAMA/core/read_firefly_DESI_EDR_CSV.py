"""
===================================================================================
|                         FIREFLY — Full Spectral Fitting                         |
===================================================================================
 <> Module: firefly_DESI_EDR_CSV.py

 <> Author:
	- Samuel Helps ~ <samuel.helps.sh__at__gmail.com>

 <> Purpose:
	- Read and update the master CSV for the DESI EDR spectra FIREFLY outputs.
___________________________________________________________________________________
|Institute of Cosmology and Gravitation - University of Portsmouth, Portsmouth, UK|
‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
"""
from astropy.io import fits
import numpy as np
import pandas as pd
import fcntl
import os
import sys

if len(sys.argv) != 4:
    print("Usage: python read_firefly_DESI_EDR_CSV.py <spectrum_index> <file_range> <temp_file>")
    sys.exit(1)

spec_index = int(sys.argv[1])
file_range = sys.argv[2]  # e.g., "0-10000"
temp_file = sys.argv[3]   # Full path to temp file

script_dir = os.path.dirname(os.path.abspath(__file__))  # .../Firefly/Launch/DESI/SCIAMA/core
repo_root = os.path.abspath(os.path.join(script_dir, "..", "..", "..", ".."))  # <repo>/Firefly

# Temp output file path for this spectrum
temp_file = os.path.join(
    repo_root, "Fitting_Engine", "temp_output",
    f"spFly-firefly_release/desidata_{file_range}",
    f"DESI_Firefly_EDR_GALAXY_{spec_index}.fits"
)
hdul = fits.open(temp_file)

# Define output directories
##################### ONLY CHANGE GALAXY INDEX ##############################
start_index, end_index = map(int, file_range.split('-'))
#############################################################################

remote_base_dir = os.path.join(repo_root, "Results", "DESI_EDR_VAC")
master_file = os.path.join(remote_base_dir, f'master_files/EDR_master_{file_range}.csv')

# Create master directory only
os.makedirs(os.path.dirname(master_file), exist_ok=True)

data = hdul[1].data
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

# Write to file
# ################################################################################################
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
        
        # SSP specific parameters
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
