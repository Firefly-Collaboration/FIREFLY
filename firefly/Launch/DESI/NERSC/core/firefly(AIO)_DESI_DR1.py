"""
===================================================================================
|                         FIREFLY — Full Spectral Fitting                         |
===================================================================================
 <> Module: firefly(AIO)_DESI_DR1.py

 <> Author:
	- Kieran Graham ~ <kieran.graham__at__port.ac.uk>

 <> Contributors:
	- Samuel Helps ~ <samuel.helps.sh__at__gmail.com>

 <> Purpose:
	- Quickest Firefly pipeline for fitting DESI DR1 'Iron' data on NERSC HPC.
___________________________________________________________________________________
|Institute of Cosmology and Gravitation - University of Portsmouth, Portsmouth, UK|
‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
"""
# ==============================================================================
# 1. IMPORTS
# ==============================================================================
# Standard Python libraries
import os
import gc
import time
import sys
import warnings
import argparse # Added for command-line arguments
from collections import defaultdict
import traceback
from pathlib import Path

# Scientific computing libraries
import numpy as np
import matplotlib.pyplot as plt

# Astropy for data handling
from astropy.io import fits
from astropy.table import Table, join
import astropy.cosmology as co
import fitsio

# DESI-specific modules
from desimodel.footprint import radec2pix
import psutil

# Multiprocessing for parallel execution
from multiprocessing import Manager
from concurrent.futures import ProcessPoolExecutor, as_completed

# Matplotlib setup
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 12
from matplotlib.font_manager import FontProperties
font1 = FontProperties(family='sans-serif', size=14, weight='regular')


# ==============================================================================
# 1.5 HELPER FUNCTIONS
# ==============================================================================

def log_memory_usage(stage):
    """Logs the current and peak memory usage of the main process."""
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    rss_gb = mem_info.rss / (1024 ** 3)
    print(f"[Memory Log] {stage}: Current RSS = {rss_gb:.2f} GB")


# ==============================================================================
# 2. FIREFLY WRAPPER LOGIC
# ==============================================================================
# This section contains the refactored, importable version of the Firefly script.

# Set up Firefly environment (this must run before importing firefly modules)
script_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.abspath(os.path.join(script_dir, "..", "..", "..", ".."))
fe_path = os.path.join(repo_root, "Fitting_Engine")
if fe_path not in sys.path:
    sys.path.append(fe_path)
os.environ["FF_DIR"] = repo_root
os.environ["STELLARPOPMODELS_DIR"] = os.path.join(fe_path, "stellar_population_models")

# Now we can import the firefly modules
import firefly_setup as fs
import firefly_models as fm

def setup_firefly():
    """
    Gathers all static Firefly configuration into a dictionary.
    This is called once per worker process.
    """
    config = {
        'cosmo': co.Planck15,
        'emlines': [
            'He-II', 'Ne-V', 'O-II', 'Ne-III', 'H-ζ', 'H-ε', 'H-δ', 'H-γ',
            'O-III', 'Ar-IV', 'H-β', 'N-I', 'He-I', 'O-I', 'N-II', 'H-α',
            'S-II', 'Ar-III'
        ],
        'N_angstrom_masked': 0,
        'model_key': 'MaStar',
        'model_lib': ['gold'],
        'imfs': ['ss'],
        'age_limits': [0, 'AoU'],
        'Z_limits': [-3., 3.],
        'data_wave_medium': 'vacuum',
        'fit_wave_medium': 'vacuum',
        'flux_units': 1e-17,
        'milky_way_reddening': True,
        'hpf_mode': 'on',
        'dust_law': 'calzetti',
        'max_ebv': 1.5,
        'num_dust_vals': 200,
        'dust_smoothing_length': 200,
        'max_iterations': 10,
        'pdf_sampling': 300,
        # IMPORTANT: Set your final output folder here (adjust as needed per-user)
        'outputFolder': "/global/cfs/cdirs/desi/users/helpss/FIREFLY/firefly/Results/DESI_DR1(AIO)_output"
    }
    # Create output directory if it doesn't exist
    if not os.path.isdir(config['outputFolder']):
        os.makedirs(config['outputFolder'], exist_ok=True)
        
    return config

def run_firefly_on_spectrum(config, galaxy_info):
    """
    Runs the Firefly fitting routine for a single galaxy spectrum.
    """
    # Unpack galaxy data
    ID = galaxy_info['TARGETID']
    redshift = galaxy_info['Z']
    ra = galaxy_info['TARGET_RA']
    dec = galaxy_info['TARGET_DEC']
    vdisp = galaxy_info['VDISP']
    wavelength = galaxy_info['wavelength']
    flux = galaxy_info['flux']
    error = np.zeros_like(galaxy_info['ivar'])
    good_ivar = galaxy_info['ivar'] > 0
    error[good_ivar] = galaxy_info['ivar'][good_ivar]**(-0.5)

    output_file = os.path.join(config['outputFolder'], f"spFly-{ID}.fits")
    if os.path.exists(output_file):
        os.remove(output_file)

    age_min = config['age_limits'][0]
    age_max = config['cosmo'].age(redshift).value if config['age_limits'][1] == 'AoU' else config['age_limits'][1]
    
    spec = fs.firefly_setup("in-memory-object",
                              milky_way_reddening=config['milky_way_reddening'],
                              N_angstrom_masked=config['N_angstrom_masked'],
                              hpf_mode=config['hpf_mode'],
                              data_wave_medium=config['data_wave_medium'])
    
    r_instrument = np.full_like(wavelength, 2000.0)
    spec.openSingleSpectrum(wavelength, flux, error, redshift, ra, dec, vdisp, 
                              config['emlines'], r_instrument)

    model = fm.StellarPopulationModel(spec, output_file, config['cosmo'], 
                                        models=config['model_key'], model_libs=config['model_lib'], imfs=config['imfs'],
                                        age_limits=[age_min, age_max], downgrade_models=True,
                                        data_wave_medium=config['data_wave_medium'], fit_wave_medium=config['fit_wave_medium'], 
                                        Z_limits=config['Z_limits'], use_downgraded_models=False, dust_law=config['dust_law'], 
                                        max_ebv=config['max_ebv'], num_dust_vals=config['num_dust_vals'],
                                        dust_smoothing_length=config['dust_smoothing_length'], max_iterations=config['max_iterations'],
                                        pdf_sampling=config['pdf_sampling'], flux_units=config['flux_units'])

    # Initiate fit and save results
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', message='divide by zero encountered', category=RuntimeWarning)
            model.fit_models_to_data()
    except ValueError as e:
        if "empty sequence" in str(e):
            # This error is handled more gracefully inside the main worker loop now
            pass
        # Re-raise the exception so the main loop can handle it
        raise e

    
    prihdr = fits.Header()
    prihdr['FILE'] = os.path.basename(output_file)
    prihdr['MODELS'] = config['model_key']
    prihdr['AGEMIN'] = str(age_min)
    prihdr['AGEMAX'] = str(age_max)
    prihdr['ZMIN'] = str(config['Z_limits'][0])
    prihdr['ZMAX'] = str(config['Z_limits'][1])
    prihdr['redshift'] = redshift
    prihdu = fits.PrimaryHDU(header=prihdr)
    
    complete_hdus = fits.HDUList([prihdu, model.tbhdu])
    complete_hdus.writeto(output_file, overwrite=True)

# ==============================================================================
# 3. WORKER INITIALISER & DIRECT-CALL FUNCTION
# ==============================================================================

def init_worker():
    """
    Initialises a worker process by setting up the Firefly configuration.
    """
    global firefly_config
    firefly_config = setup_firefly()

def process_coadd_chunk(coadd_task_chunk, counter, completed_ids):
    """
    Worker function that processes an assigned CHUNK of targets from a single coadd file.
    It now returns a list of result dictionaries instead of printing directly.
    """
    import desispec.io
    from desispec import coaddition

    global firefly_config
    coadd_path = coadd_task_chunk['coadd_path']
    targets_in_chunk = coadd_task_chunk['targets_in_chunk']
    
    if not os.path.exists(coadd_path):
        results = []
        for target_info in targets_in_chunk:
            ID = target_info['TARGETID']
            results.append({'id': ID, 'status': 'Error: Coadd file not found', 'duration': 0})
            with counter['lock']:
                counter['errors'] += 1
            completed_ids[ID] = True
        return results

    try:
        target_ids_to_read = [t['TARGETID'] for t in targets_in_chunk]
        coadd_obj = desispec.io.read_spectra(coadd_path, targetids=target_ids_to_read, skip_hdus=['RESOLUTION'])
    except Exception as e:
        results = []
        for target_info in targets_in_chunk:
            ID = target_info['TARGETID']
            results.append({'id': ID, 'status': f'Error: Failed to read coadd file - {e}', 'duration': 0})
            with counter['lock']:
                counter['errors'] += 1
            completed_ids[ID] = True
        return results

    results = []
    for target_info in targets_in_chunk:
        start_time = time.time()
        ID = target_info['TARGETID']

        try:
            coadd_spec = coadd_obj.select(targets=[ID])
            if len(coadd_spec.flux) == 0:
                raise ValueError("Spectrum not found in coadd object.")

            spec_combined = coaddition.coadd_cameras(coadd_spec)
            ivar = spec_combined.ivar['brz'][0]
            target_info['wavelength'] = spec_combined.wave['brz']
            target_info['flux'] = spec_combined.flux['brz'][0]
            target_info['ivar'] = np.where(ivar == 0, 1e-8, ivar)
            
            run_firefly_on_spectrum(firefly_config, target_info)
            
            expected_output = os.path.join(firefly_config['outputFolder'], f"spFly-{ID}.fits")
            if os.path.exists(expected_output):
                with counter['lock']:
                    counter['completed'] += 1
                status = "Complete"
            else:
                raise RuntimeError("Firefly ran but output file was not created.")
        
        except Exception as e:
            with counter['lock']:
                counter['errors'] += 1
            if "wavelength coverage" in str(e):
                status = f"Error: Wavelength mismatch (z={target_info.get('Z', 'N/A'):.2f})"
            else:
                status = f"Error: {e}"

        duration = time.time() - start_time
        completed_ids[ID] = True
        
        results.append({'id': ID, 'status': status, 'duration': duration})
        
    return results


# ==============================================================================
# 4. MAIN SCRIPT LOGIC
# ==============================================================================
if __name__ == '__main__':

    # Setup command-line argument parsing
    parser = argparse.ArgumentParser(description="Run Firefly in parallel on DESI DR1 Iron data.")
    parser.add_argument('--start', type=int, required=True, help='The starting row index to process from the fastspecfit catalogue.')
    parser.add_argument('--end', type=int, required=True, help='The ending row index (exclusive) to process.')
    args = parser.parse_args()

    log_memory_usage("Script Start")

    # --- Load and filter initial catalogues (DR1 Iron Version) ---
    specprod = 'iron'
    specprod_dir = '/dvs_ro/cfs/cdirs/desi/public/dr1/spectro/redux/iron'
    healpix_dir = f'{specprod_dir}/healpix'
    
    # Use command-line arguments for the processing range
    start_index = args.start
    end_index = args.end
    print(f"Processing fastspecfit rows from {start_index} to {end_index-1}")

    print("Gathering slice of Fastspecfit Catalog...")
    fastspec_path = '/dvs_ro/cfs/cdirs/desi/spectro/fastspecfit/iron/v2.1/catalogs/fastspec-iron.fits'
    rows_to_read = np.arange(start_index, end_index)
    batch_spec_table = Table(fitsio.read(fastspec_path, columns=['TARGETID', 'VDISP'], rows=rows_to_read))
    print(f"Finished gathering {len(batch_spec_table)} rows from Fastspecfit Catalog.")

    print("Gathering Z Catalog to get coordinates and survey info...")
    zpix_cat = Table.read(f"{specprod_dir}/zcatalog/v1/zall-pix-{specprod}.fits", hdu="ZCATALOG")
    print("Finished Gathering Z Catalog")
    log_memory_usage("After Full Catalogue Load")

    # --- High-performance task preparation ---
    print("\nPreparing task list by joining and filtering...")
    prep_start_time = time.time()

    # Join the small slice of fastspecfit data with the large zcatalog
    joined_table_full = join(batch_spec_table, zpix_cat, keys='TARGETID', join_type='inner')
    
    # CRUCIAL: Immediately delete the large zcatalog to free memory
    del zpix_cat
    gc.collect()
    log_memory_usage("After Join and zcat Deletion")
    
    # Now, perform the final filtering on the much smaller, joined table
    galaxy_mask = joined_table_full['SPECTYPE'] == 'GALAXY'
    single_mask = joined_table_full['ZCAT_NSPEC'] == 1
    pos_mask = joined_table_full['TARGETID'] > 0
    final_table = joined_table_full[galaxy_mask & single_mask & pos_mask]

    # Handle missing VDISP values
    if 'VDISP' in final_table.colnames and hasattr(final_table['VDISP'], 'filled'):
        final_table['VDISP'] = final_table['VDISP'].filled(125.0)

    num_targets_in_batch = len(final_table)
    print(f"Selected {num_targets_in_batch} valid galaxy targets for processing.")
    
    # --- Group targets by coadd file and create "chunks" for balanced workload ---
    print("Grouping tasks by coadd file and creating chunks...")
    TARGETS_PER_CHUNK = 10  # Tune this parameter for optimal load balancing
    tasks_by_file = defaultdict(list)
    for row in final_table:
        hpx = row['HEALPIX']
        survey = row['SURVEY']
        program = row['PROGRAM']
        tgt_dir = f'{healpix_dir}/{survey}/{program}/{hpx//100}/{hpx}'
        coadd_path = f'{tgt_dir}/coadd-{survey}-{program}-{hpx}.fits'
        
        tasks_by_file[coadd_path].append({
            'TARGETID': row['TARGETID'], 'Z': row['Z'], 'TARGET_RA': row['TARGET_RA'],
            'TARGET_DEC': row['TARGET_DEC'], 'VDISP': row['VDISP']
        })

    coadd_task_chunks = []
    for path, targets in tasks_by_file.items():
        if len(targets) <= TARGETS_PER_CHUNK:
            coadd_task_chunks.append({'coadd_path': path, 'targets_in_chunk': targets})
        else:
            for i in range(0, len(targets), TARGETS_PER_CHUNK):
                chunk = targets[i:i + TARGETS_PER_CHUNK]
                coadd_task_chunks.append({'coadd_path': path, 'targets_in_chunk': chunk})
    
    prep_end_time = time.time()
    print(f"Finished preparing {len(coadd_task_chunks)} balanced task chunks in {prep_end_time - prep_start_time:.2f} seconds.")

    # --- Clean up large tables to free up memory before forking ---
    print("\nClearing large catalogues from memory...")
    del joined_table_full, final_table, tasks_by_file, batch_spec_table
    gc.collect()
    log_memory_usage("After Cleanup (Final Pre-Fork)")

    # ==========================================================================
    # 5. RESILIENT PARALLEL EXECUTION BLOCK
    # ==========================================================================
    MAX_FAILURES = 2
    total_chunks = len(coadd_task_chunks)
    
    process = psutil.Process(os.getpid())
    base_mem_gb = process.memory_info().rss / (1024 ** 3)
    
    start_time_workers = time.time()

    with Manager() as manager:
        counter = manager.dict({'completed': 0, 'errors': 0, 'lock': manager.Lock()})
        completed_ids = manager.dict()
        max_workers = 128

        with ProcessPoolExecutor(max_workers=max_workers, initializer=init_worker) as executor:
            future_to_chunk = {
                executor.submit(process_coadd_chunk, chunk, counter, completed_ids): chunk
                for chunk in coadd_task_chunks
            }

            print(f"\nPool started. Submitting {total_chunks} task chunks to {max_workers} workers.")
            
            for i, future in enumerate(as_completed(future_to_chunk)):
                chunk_results = future.result()
                
                print(f"--- Chunk {i+1}/{total_chunks} Finished ---")
                
                for res in chunk_results:
                    progress = len(completed_ids)
                    target_id = res['id']
                    status = res['status']
                    duration = res['duration']
                    
                    print(f"[{progress}/{num_targets_in_batch}] ID {target_id}: {status} | "
                          f"Done: {counter['completed']}, Errors: {counter['errors']} | Time: {duration:.2f}s")

        worker_duration = time.time() - start_time_workers
        final_target_count = len(completed_ids)

    print(f"\n\n--- Job Summary ---")
    print(f"Processing finished in {worker_duration:.2f} seconds.")
    print(f"Total targets marked as complete or skipped: {final_target_count}/{num_targets_in_batch}")