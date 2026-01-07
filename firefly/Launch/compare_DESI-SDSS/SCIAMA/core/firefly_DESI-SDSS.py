"""
===================================================================================
|                         FIREFLY — Full Spectral Fitting                         |
===================================================================================

 <> Module: firefly_DESI-SDSS.py

 <> Author:
	- Samuel Helps ~ <samuel.helps.sh__at__gmail.com>

 <> Purpose:
	- Run FIREFLY on matched DESI and SDSS spectra for comparison (dual-fitting).

 <> Acknowledgements:
    - Spectra retrieved via SPARCL (sparcl.client). Please cite SPARCL developers
      and data providers when publishing.
___________________________________________________________________________________
|Institute of Cosmology and Gravitation - University of Portsmouth, Portsmouth, UK|
‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
"""
import sys, os, subprocess
import sys, os

script_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.abspath(os.path.join(script_dir, "..", "..", ".."))  # firefly\Launch
sys.path.append(os.path.join(os.path.dirname(repo_root), "Fitting_Engine"))  # firefly\Fitting_Engine
os.environ["FF_DIR"] = os.path.dirname(repo_root)  # firefly
os.environ["STELLARPOPMODELS_DIR"] = os.path.join(os.environ["FF_DIR"], "Fitting_Engine", "stellar_population_models")

import numpy as np
from astropy.io import fits
import astropy.cosmology as co
import firefly_setup as fs
import firefly_models as fm
import time
import warnings
from astropy.io.fits.verify import VerifyWarning

warnings.filterwarnings('ignore', category=VerifyWarning, message='.*greater than 8 characters.*')

t0 = time.time()
cosmo = co.Planck15

#########################################
#    Load environment and inputs       #
#########################################
desi_input = os.environ.get('DESI_INPUT_FILE')
sdss_input = os.environ.get('SDSS_INPUT_FILE')
if not desi_input or not sdss_input:
    print("DESI_INPUT_FILE or SDSS_INPUT_FILE environment variable not set")
    sys.exit(1)

base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Parse arguments for spectrum indices
if len(sys.argv) != 3:
    print("Usage: python firefly_DESI-SDSS.py <start_index> <end_index>")
    sys.exit(1)

start_index = int(sys.argv[1])
end_index = int(sys.argv[2])

# Open FITS files for DESI and SDSS
hdul_desi = fits.open(desi_input)
hdul_sdss = fits.open(sdss_input)

# Adjust end_index if beyond available spectra
actual_desi = len(hdul_desi[1].data)
actual_sdss = len(hdul_sdss[1].data)
actual_spectra = min(actual_desi, actual_sdss)
if end_index > actual_spectra:
    print(f"Warning: Requested end index {end_index} exceeds available spectra ({actual_spectra}). Adjusting.")
    end_index = actual_spectra
spectra_to_process = range(start_index, end_index)

print(f"\nProcessing spectra indices: {list(spectra_to_process)}")

# FIREFLY Fitting Engine temp_output
outputFolder = os.path.join(os.environ["FF_DIR"], "Fitting_Engine", "temp_output")
os.makedirs(outputFolder, exist_ok=True)

# Temporary directory per sample/file range
file_range = os.path.basename(desi_input).split('_')[2].split('.')[0]
temp_dir = os.path.join(outputFolder, f'compare_{file_range}')
os.makedirs(temp_dir, exist_ok=True)

failed_indices = []

for i in spectra_to_process:
    try:
        print(f"\n=== Processing spectrum index {i} ===")

        # --- SDSS Data ---
        wavelength_sdss = np.array(hdul_sdss[1].data['WAVE'][i])
        flux_sdss = np.array(hdul_sdss[1].data['FLUX'][i])
        ivar_sdss = np.array(hdul_sdss[1].data['IVAR'][i])

        # Compute error safely: where ivar>0 take ivar**-0.5 else 0.0
        error_sdss = np.where(ivar_sdss > 0, ivar_sdss ** -0.5, 0.0)

        # SDSS metadata 
        try:
            redshift_sdss = hdul_sdss[1].data['REDSHIFT'][i]
        except Exception:
            redshift_sdss = hdul_sdss[1].header.get('redshift', -1)
        ra_sdss = hdul_sdss[1].data['RA'][i] if 'RA' in hdul_sdss[1].columns.names else hdul_sdss[0].header.get('RA', np.nan)
        dec_sdss = hdul_sdss[1].data['DEC'][i] if 'DEC' in hdul_sdss[1].columns.names else hdul_sdss[0].header.get('DEC', np.nan)
        vdisp_sdss = hdul_sdss[1].data['VDISP'][i] if 'VDISP' in hdul_sdss[1].columns.names else 0.0
        survey_sdss = 'SDSS'

        # SDSS SNR: prefer per-row table value
        snr_sdss = hdul_sdss[0].header.get('SNR', -1)
        try:
            # Attempt to read table column 'SNR' if present
            if 'SNR' in hdul_sdss[1].columns.names:
                snr_val = hdul_sdss[1].data['SNR'][i]
                snr_sdss = float(np.atleast_1d(snr_val)[0])
            else:
                cols = [c.lower() for c in hdul_sdss[1].columns.names]
                if 'snr' in cols:
                    idx = cols.index('snr')
                    colname = hdul_sdss[1].columns.names[idx]
                    snr_val = hdul_sdss[1].data[colname][i]
                    snr_sdss = float(np.atleast_1d(snr_val)[0])
        except Exception:
            pass

        # --- DESI Data ---
        wavelength_desi = np.array(hdul_desi[1].data['WAVE'][i])
        flux_desi = np.array(hdul_desi[1].data['FLUX'][i])
        ivar_desi = np.array(hdul_desi[1].data['IVAR'][i])

        # Safe error for DESI: where ivar>0 use 1/sqrt(ivar), else 0.0
        error_desi = np.where(ivar_desi > 0, 1.0 / np.sqrt(ivar_desi), 0.0)

        redshift_desi = hdul_desi[1].data['REDSHIFT'][i]
        ra_desi = hdul_desi[1].data['RA'][i]
        dec_desi = hdul_desi[1].data['DEC'][i]
        vdisp_desi = hdul_desi[1].data['VDISP'][i]
        survey_desi = 'DESI'

        # DESI SNR: per-row table value preferred, fallback to header
        snr_desi = hdul_desi[0].header.get('SNR', -1)
        try:
            if 'SNR' in hdul_desi[1].columns.names:
                snr_val_d = hdul_desi[1].data['SNR'][i]
                snr_desi = float(np.atleast_1d(snr_val_d)[0])
            else:
                cols_d = [c.lower() for c in hdul_desi[1].columns.names]
                if 'snr' in cols_d:
                    idxd = cols_d.index('snr')
                    colname_d = hdul_desi[1].columns.names[idxd]
                    snr_val_d = hdul_desi[1].data[colname_d][i]
                    snr_desi = float(np.atleast_1d(snr_val_d)[0])
        except Exception:
            pass

        # --- Apply wavelength cut if requested in sbatch settings ---
        cut_mode = os.environ.get('CUT_MODE', '').lower()
        if cut_mode == 'cut':
            # Identify non-zero range in SDSS
            valid = (flux_sdss != 0) & ~np.isnan(flux_sdss)
            if np.any(valid):
                min_wave = np.min(wavelength_sdss[valid])
                max_wave = np.max(wavelength_sdss[valid])
                # Crop SDSS
                mask_sdss = (wavelength_sdss >= min_wave) & (wavelength_sdss <= max_wave)
                wavelength_sdss = wavelength_sdss[mask_sdss]
                flux_sdss = flux_sdss[mask_sdss]
                error_sdss = error_sdss[mask_sdss]
                # Crop DESI to same range (no resampling)
                mask_desi = (wavelength_desi >= min_wave) & (wavelength_desi <= max_wave)
                wavelength_desi = wavelength_desi[mask_desi]
                flux_desi = flux_desi[mask_desi]
                error_desi = error_desi[mask_desi]
                print(f"Applied wavelength cut: {min_wave:.1f} - {max_wave:.1f} Å")
            else:
                print("No valid SDSS data for cutting; skipping cut.")
        else:
            print("No cut applied; using full wavelength range.")

        # Check enough data points
        if len(wavelength_sdss) < 10 or len(wavelength_desi) < 10:
            print(f"Insufficient data after cut for index {i}. Skipping.")
            failed_indices.append(i)
            continue

        # Firefly parameters
        r_instrument = np.full_like(wavelength_desi, 5000.0)
        r_inst_sdss = np.full_like(wavelength_sdss, 2000.0)

        model_key = 'MaStar'
        model_lib = ['gold']
        imfs = ['kr']
        age_limits = [0, 'AoU']
        Z_limits = [-3.0, 3.0]
        data_wave_medium = 'vacuum'
        fit_wave_medium = 'vacuum'
        flux_unit = 10**(-17)
        milky_way_reddening = True
        hpf_mode = 'on'
        dust_law = 'calzetti'
        max_ebv = 1.5
        num_dust_vals = 200
        dust_smoothing_length = 200
        max_iterations = 10
        pdf_sampling = 300

        # Output FITS files for each survey
        temp_file_sdss = os.path.join(temp_dir, f'SDSS_Firefly_output_{i}.fits')
        temp_file_desi = os.path.join(temp_dir, f'DESI_Firefly_output_{i}.fits')

        failed = False
        # --- Run Firefly for SDSS ---
        print("Running FIREFLY on SDSS spectrum...")
        prihdr_sdss = fm.pyfits.Header()
        prihdr_sdss['FILE'] = os.path.basename(temp_file_sdss)
        prihdr_sdss['MODELS'] = model_key
        prihdr_sdss['FITTER'] = "FIREFLY"
        prihdr_sdss['redshift'] = redshift_sdss
        prihdr_sdss['ID'] = hdul_sdss[0].header.get('PLATE', -1)  # or other ID field
        prihdr_sdss['RA'] = ra_sdss
        prihdr_sdss['DEC'] = dec_sdss
        prihdr_sdss['SNR'] = snr_sdss  # write per-spectrum SNR into header
        prihdu_sdss = fm.pyfits.PrimaryHDU(header=prihdr_sdss)
        tables_sdss = [prihdu_sdss]

        spec_sdss = fs.firefly_setup(
            sdss_input,
            milky_way_reddening=milky_way_reddening,
            N_angstrom_masked=0,
            hpf_mode=hpf_mode,
            data_wave_medium=data_wave_medium
        )
        spec_sdss.openSingleSpectrum(wavelength_sdss, flux_sdss, error_sdss,
                                     redshift_sdss, ra_sdss, dec_sdss,
                                     vdisp_sdss, [], r_inst_sdss)
        try:
            model_sdss = fm.StellarPopulationModel(
                spec_sdss, temp_file_sdss, cosmo,
                models=model_key, model_libs=model_lib, imfs=imfs,
                age_limits=[age_limits[0], cosmo.age(redshift_sdss).value],
                data_wave_medium=data_wave_medium, fit_wave_medium=fit_wave_medium,
                Z_limits=Z_limits,
                suffix='', use_downgraded_models=False,
                dust_law=dust_law, max_ebv=max_ebv, num_dust_vals=num_dust_vals,
                dust_smoothing_length=dust_smoothing_length,
                max_iterations=max_iterations, pdf_sampling=pdf_sampling,
                flux_units=flux_unit
            )
            model_sdss.fit_models_to_data()
            model_sdss.tbhdu.header['SPEC_IDX'] = i
            tables_sdss.append(model_sdss.tbhdu)
            complete_hdus_sdss = fm.pyfits.HDUList(tables_sdss)
            complete_hdus_sdss.writeto(temp_file_sdss, overwrite=True)
            complete_hdus_sdss.close()
            print("SDSS FIREFLY fit complete.")
        except Exception as e:
            print(f"SDSS spectrum {i} did not converge: {e}")
            failed = True

        # --- Run Firefly for DESI ---
        print("Running FIREFLY on DESI spectrum...")
        prihdr_desi = fm.pyfits.Header()
        prihdr_desi['FILE'] = os.path.basename(temp_file_desi)
        prihdr_desi['MODELS'] = model_key
        prihdr_desi['FITTER'] = "FIREFLY"
        prihdr_desi['redshift'] = redshift_desi
        prihdr_desi['ID'] = hdul_desi[1].data['ID'][i] if 'ID' in hdul_desi[1].columns.names else -1
        prihdr_desi['RA'] = ra_desi
        prihdr_desi['DEC'] = dec_desi
        prihdr_desi['SNR'] = snr_desi  # write per-spectrum SNR into header
        prihdu_desi = fm.pyfits.PrimaryHDU(header=prihdr_desi)
        tables_desi = [prihdu_desi]

        spec_desi = fs.firefly_setup(
            desi_input,
            milky_way_reddening=milky_way_reddening,
            N_angstrom_masked=0,
            hpf_mode=hpf_mode,
            data_wave_medium=data_wave_medium
        )
        spec_desi.openSingleSpectrum(wavelength_desi, flux_desi, error_desi,
                                     redshift_desi, ra_desi, dec_desi,
                                     vdisp_desi, [], r_instrument)
        try:
            model_desi = fm.StellarPopulationModel(
                spec_desi, temp_file_desi, cosmo,
                models=model_key, model_libs=model_lib, imfs=imfs,
                age_limits=[age_limits[0], cosmo.age(redshift_desi).value],
                data_wave_medium=data_wave_medium, fit_wave_medium=fit_wave_medium,
                Z_limits=Z_limits,
                suffix='', use_downgraded_models=False,
                dust_law=dust_law, max_ebv=max_ebv, num_dust_vals=num_dust_vals,
                dust_smoothing_length=dust_smoothing_length,
                max_iterations=max_iterations, pdf_sampling=pdf_sampling,
                flux_units=flux_unit
            )
            model_desi.fit_models_to_data()
            model_desi.tbhdu.header['SPEC_IDX'] = i
            tables_desi.append(model_desi.tbhdu)
            complete_hdus_desi = fm.pyfits.HDUList(tables_desi)
            complete_hdus_desi.writeto(temp_file_desi, overwrite=True)
            complete_hdus_desi.close()
            print("DESI FIREFLY fit complete.")
        except Exception as e:
            print(f"DESI spectrum {i} did not converge: {e}")
            failed = True

        if failed:
            failed_indices.append(i)
            # Clean up any output files for failed galaxies
            for f in [temp_file_sdss, temp_file_desi]:
                if os.path.exists(f):
                    os.remove(f)
            continue

        # --- Run read/analysis script ---
        print("Running comparison analysis script...")
        read_script = os.path.join(os.environ["FF_DIR"], "Launch", "compare_DESI-SDSS", "SCIAMA", "core", "read_firefly_DESI-SDSS.py")
        cmd = [
            sys.executable,
            read_script,
            str(i),
            file_range,
            temp_file_desi,
            temp_file_sdss
        ]
        subprocess.run(cmd, check=True)

        # Remove any temp files after analysis
        for f in [temp_file_sdss, temp_file_desi]:
            if os.path.exists(f):
                os.remove(f)

    except Exception as e:
        print(f"Failed to process spectrum {i}: {e}")
        failed_indices.append(i)
        continue

# Summary
print("\nProcessing complete!")
print(f"Successfully processed: {len(spectra_to_process) - len(failed_indices)} spectra")
if failed_indices:
    print(f"Failed to process: {len(failed_indices)} spectra")
    print("Failed indices:", failed_indices)
print(f"Total time: {int(time.time()-t0)} seconds")

hdul_desi.close()
hdul_sdss.close()
