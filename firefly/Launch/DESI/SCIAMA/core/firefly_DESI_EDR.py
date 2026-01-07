"""
===================================================================================
|                         FIREFLY — Full Spectral Fitting                         |
===================================================================================
 <> Module: firefly_DESI_EDR.py

 <> Author:
	- Samuel Helps ~ <samuel.helps.sh__at__gmail.com>

 <> Purpose:
	- Initialises Firefly for DESI EDR data, runs per spectrum in parallel jobs
      when submitted to slurm via sbatch.
___________________________________________________________________________________
|Institute of Cosmology and Gravitation - University of Portsmouth, Portsmouth, UK|
‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
"""

import sys, os
import subprocess

script_dir = os.path.dirname(os.path.abspath(__file__))                    # .../Launch/DESI/SCIAMA/core
repo_root = os.path.abspath(os.path.join(script_dir, "..", "..", "..", ".."))   # <repo>/firefly
fitting_engine_path = os.path.join(repo_root, "Fitting_Engine")
if fitting_engine_path not in sys.path:
    sys.path.append(fitting_engine_path)
os.environ["FF_DIR"] = repo_root
os.environ["STELLARPOPMODELS_DIR"] = os.path.join(fitting_engine_path, "stellar_population_models")
base_dir = repo_root

import numpy as np
from astropy.io import fits
import astropy.cosmology as co
import firefly_setup as fs
import firefly_models as fm
import time
import warnings
from astropy.io.fits.verify import VerifyWarning
warnings.filterwarnings('ignore', category=VerifyWarning, message='.*greater than 8 characters.*')

t0=time.time()
cosmo = co.Planck15

################################
#      CSV ONLY setting        # Will only output the master CSV files (useful for creating large VACs)
csv_only = False               # Set to True to run read_firefly_DESI_EDR_CSV.py
################################

######################################## ONLY CHANGE INPUT FILE #############################################
input_file = os.environ.get('DESI_INPUT_FILE')
if not input_file:
    print("DESI_INPUT_FILE environment variable not set")
    sys.exit(1)

input_file = os.path.join(base_dir, 'Data', 'DESI', 'DESI_EDR_data', input_file)
############################################################################################################

hdul = fits.open(input_file)
suffix = ""

# Parse arguments
if len(sys.argv) != 3:
    print("Usage: python firefly_DESI_EDR.py <start_index> <end_index>")
    sys.exit(1)

start_index = int(sys.argv[1])
end_index = int(sys.argv[2])

# Add check for actual number of spectra
actual_spectra = len(hdul[1].data)
if end_index > actual_spectra:
    print(f"Warning: Requested end index {end_index} exceeds available spectra ({actual_spectra})")
    end_index = actual_spectra
spectra_to_process = range(start_index, end_index)

# #######################################################################################
print(f"\nProcessing spectra: {list(spectra_to_process)}")

outputFolder = os.path.join(repo_root, "Fitting_Engine", "temp_output")
os.makedirs(outputFolder, exist_ok=True)
file_range = os.path.basename(input_file).split('_')[2].split('.')[0]  # Gets '0-10000'
tempFolder = os.path.join(outputFolder, f'spFly-firefly_release/desidata_{file_range}')
os.makedirs(tempFolder, exist_ok=True)
output_file = os.path.join(tempFolder, 'DESI_Firefly_EDR_GALAXY.fits')

failed_indices = []
for i in spectra_to_process:
    try:
        print(f"\nProcessing spectrum {i}")
        
        # Get spectrum data
        wavelength = hdul[1].data['WAVE'][i]
        flux = hdul[1].data['FLUX'][i]
        error = hdul[1].data['IVAR'][i]
        redshift = hdul[1].data['REDSHIFT'][i]
        ra = hdul[1].data['RA'][i]
        dec = hdul[1].data['DEC'][i]
        vdisp = hdul[1].data['VDISP'][i]
        snr = hdul[1].data['SNR'][i]
        id = hdul[1].data['ID'][i]
        survey_type = hdul[1].data['SURVEY_TYPE'][i]

        # Setup parameters
        r_instrument = np.zeros(len(wavelength))
        r_instrument.fill(5000)

        # Model parameters
        model_key = 'MaStar' #'m11'
        model_lib = ['gold'] #['MILES'] 
        imfs = ['kr']
        age_limits = [0, 'AoU'] 
        Z_limits = [-3., 3.]
        
        # Data parameters
        data_wave_medium = 'vacuum'
        fit_wave_medium = 'vacuum'
        flux_units = 10**(-17)
        
        # Dust parameters
        milky_way_reddening = True
        hpf_mode = 'on'
        dust_law = 'calzetti'
        max_ebv = 1.5
        num_dust_vals = 200
        dust_smoothing_length = 200
        max_iterations = 10
        pdf_sampling = 300

        # Emission lines to mask
        N_angstrom_masked = 0
        emlines = ['He-II', 'Ne-V', 'O-II', 'Ne-III', 'H-ζ', 'H-ε', 'H-δ', 'H-γ',
                  'O-III', 'Ar-IV', 'H-β', 'N-I', 'He-I', 'O-I', 'N-II', 'H-α',
                  'S-II', 'Ar-III']

        print('\nStarting firefly ...')

        # Set age limits
        age_min = age_limits[0]
        age_max = cosmo.age(redshift).value if age_limits[1] == 'AoU' else age_limits[1]
        Z_min, Z_max = Z_limits

        # Create primary HDU
        prihdr = fm.pyfits.Header()
        prihdr['FILE'] = os.path.basename(output_file)
        prihdr['MODELS'] = model_key
        prihdr['FITTER'] = "FIREFLY"
        prihdr['AGEMIN'] = str(age_min)
        prihdr['AGEMAX'] = str(age_max)
        prihdr['ZMIN'] = str(Z_min)
        prihdr['ZMAX'] = str(Z_max)
        prihdr['redshift'] = redshift
        prihdr['SNR'] = snr
        prihdr['ID'] = id
        prihdr['RA'] = ra
        prihdr['DEC'] = dec
        prihdr['SURVEY_TYPE'] = survey_type
        prihdr['HIERARCH age_universe'] = np.round(cosmo.age(redshift).value, 3)
        prihdu = fm.pyfits.PrimaryHDU(header=prihdr)
        tables = [prihdu]

        # Run FIREFLY
        spec = fs.firefly_setup(
            input_file, 
            milky_way_reddening=milky_way_reddening,
            N_angstrom_masked=N_angstrom_masked,
            hpf_mode=hpf_mode,
            data_wave_medium=data_wave_medium
        )
        
        spec.openSingleSpectrum(wavelength, flux, error, redshift, ra, dec, vdisp, emlines, r_instrument)

        try:
            model = fm.StellarPopulationModel(
                spec, output_file, cosmo,
                models=model_key,
                model_libs=model_lib,
                imfs=imfs,
                age_limits=[age_min, age_max],
                downgrade_models=True,
                data_wave_medium=data_wave_medium,
                fit_wave_medium=fit_wave_medium,
                Z_limits=Z_limits,
                suffix='',
                use_downgraded_models=False,
                dust_law=dust_law,
                max_ebv=max_ebv,
                num_dust_vals=num_dust_vals,
                dust_smoothing_length=dust_smoothing_length,
                max_iterations=max_iterations,
                pdf_sampling=pdf_sampling,
                flux_units=flux_units
            )

            model.fit_models_to_data()
            model.tbhdu.header['SPEC_IDX'] = i
            tables.append(model.tbhdu)

            # Save results and overwrite existing file
            complete_hdus = fm.pyfits.HDUList(tables)
            temp_file = os.path.join(tempFolder, f'DESI_Firefly_EDR_GALAXY_{i}.fits')
            complete_hdus.writeto(temp_file, overwrite=True)

            # Run LOOP_read_firefly_desi.py after each successful fit
            print("\nRunning analysis...")
            if csv_only == True:
                cmd = [sys.executable, 
                    os.path.join(repo_root, "Launch", "DESI", "SCIAMA", "core", "LOOP_read_firefly_desi_CSV.py"),
                    str(i),  # spectrum index
                    file_range,  # range like "0-10000"
                    temp_file]  # actual temp file path
            else:
                cmd = [sys.executable, 
                    os.path.join(repo_root, "Launch", "DESI", "SCIAMA", "core", "LOOP_read_firefly_desi.py"),
                    str(i),  # spectrum index
                    file_range,  # range like "0-10000"
                    temp_file]  # actual temp file path
            subprocess.run(cmd, check=True)
            
            # Clean up temp file
            if os.path.exists(temp_file):
                os.remove(temp_file)

            complete_hdus.close()

        except ValueError:
            print(f'Spectrum {i} did not converge')
            continue

    except Exception as e:
        print(f"Failed to process spectrum {i}: {str(e)}")
        failed_indices.append(i)
        continue

# Print summary
print("\nProcessing complete!")
print(f"Successfully processed: {len(spectra_to_process) - len(failed_indices)} spectra")
if failed_indices:
    print(f"Failed to process: {len(failed_indices)} spectra")
    print("Failed indices:", failed_indices)

print(f"\nTotal time: {int(time.time()-t0)} seconds")
hdul.close()
