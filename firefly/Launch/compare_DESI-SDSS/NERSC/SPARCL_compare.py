"""
===================================================================================
|                         FIREFLY — Full Spectral Fitting                         |
===================================================================================

 <> Module: SPARCL_compare.py

 <> Author:
	- Samuel Helps ~ <samuel.helps.sh__at__gmail.com>

 <> Purpose:
	- Query SPARCL, match DESI and SDSS spectra, resample/trim them, and write
      Firefly-ready FITS files.

 <> Acknowledgements:
    - Spectra retrieved via SPARCL (sparcl.client). Please cite SPARCL developers
      and data providers when publishing.
___________________________________________________________________________________
|Institute of Cosmology and Gravitation - University of Portsmouth, Portsmouth, UK|
‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
"""
import os
import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord, search_around_sky
import astropy.units as u
from astropy.io import fits
from scipy.interpolate import interp1d

from sparcl.client import SparclClient

# -------------------- SETTINGS --------------------
RA_MIN = 208.4
RA_MAX = 210.2
DEC_MIN = 4.8
DEC_MAX = 6.4
Z_MIN = 0.1
Z_MAX = 0.15

MATCH_RADIUS_ARCSEC = 0.5  # Set max difference in distance for matching objects
LIMIT = 5000
OUT_PREFIX = "greatwall_sample"

FASTSPEC_PATH = "/global/cfs/cdirs/desi/public/edr/vac/edr/fastspecfit/fuji/v3.2/catalogs/fastspec-fuji.fits"

# Output directory
OUT_DIR = os.path.join("firefly", "Data", "DESI-SDSS_samples")
os.makedirs(OUT_DIR, exist_ok=True)

# May need to resample if wavelength arrays differ slightly
FORCE_RESAMPLE = False

# Minimum number of good pixels after trimming to keep a spectrum:
MIN_GOOD_PIXELS = 10
# --------------------------------------------------

def load_fastspec_vdisp(path):
    if not path:
        return {}
    try:
        t = Table.read(path, hdu=1)
        if 'TARGETID' in t.colnames and 'VDISP' in t.colnames:
            return {int(tid): float(v) for tid, v in zip(t['TARGETID'], t['VDISP']) if np.isfinite(v)}
    except Exception as e:
        print("Warning: cannot read fastspec file:", e)
    return {}

def compute_snr_from_flux_ivar(flux, ivar):
    f = np.asarray(flux)
    iv = np.asarray(ivar)
    good = (iv > 0) & np.isfinite(f)
    if np.any(good):
        return float(np.median(np.abs(f[good] * np.sqrt(iv[good]))))
    return np.nan

def estimate_vdisp_from_resolution(lam, R):
    lam = np.asarray(lam)
    if R is None:
        return 125.0
    R = np.asarray(R)
    if R.size == 1:
        Rarr = R * np.ones_like(lam)
    else:
        if len(R) == len(lam):
            Rarr = R
        else:
            try:
                xold = np.linspace(lam.min(), lam.max(), len(R))
                Rarr = np.interp(lam, xold, R, left=np.nanmedian(R), right=np.nanmedian(R))
            except Exception:
                Rarr = np.nanmedian(R) * np.ones_like(lam)
    sigma_lambda = (lam / Rarr) / 2.355
    with np.errstate(invalid='ignore'):
        v_disp = 3e5 * np.nanmedian(sigma_lambda / lam)
    if not np.isfinite(v_disp) or v_disp < 20:
        return 125.0
    return float(v_disp)

def all_waves_equal(wave_list, tol=1e-8):
    base = wave_list[0]
    for w in wave_list[1:]:
        if len(w) != len(base):
            return False
        if not np.allclose(w, base, atol=tol, rtol=tol):
            return False
    return True

def per_spectrum_trim(base_wave, flux_res, ivar_res):
    """
    Given resampled arrays on base_wave, find first/last good pixel for this spectrum
    Good pixel defined as: finite(flux) AND ivar>0 AND flux != 0
    Return: trimmed_wave, trimmed_flux, trimmed_ivar
    If no good pixels, return empty arrays.
    """
    good = np.isfinite(flux_res) & (ivar_res > 0) & (flux_res != 0.0)
    if not np.any(good):
        return np.array([], dtype=np.float64), np.array([], dtype=np.float64), np.array([], dtype=np.float64)
    first = np.argmax(good)  # index of first True
    last = len(good) - 1 - np.argmax(good[::-1])
    return base_wave[first:last+1].astype(np.float64), flux_res[first:last+1].astype(np.float64), ivar_res[first:last+1].astype(np.float64)

def build_varlen_hdu(spec_records, label, waves_list, identical_flag, fastspec_map):
    N = len(spec_records.records)

    # Choose a base wave grid if necessary (resample onto it)
    if identical_flag:
        base_wave = np.array(waves_list[0])
    else:
        mins = np.array([np.nanmin(w) for w in waves_list])
        maxs = np.array([np.nanmax(w) for w in waves_list])
        deltas = np.array([np.median(np.diff(w)) for w in waves_list if len(w) > 2])
        step = float(np.nanmedian(deltas)) if deltas.size > 0 else 1.0
        wave_min = float(np.nanmin(mins))
        wave_max = float(np.nanmax(maxs))
        print(f"WARNING: {label} WAVE arrays differ; resampling to linear grid {wave_min:.1f}-{wave_max:.1f} Å step {step:.3f}")
        base_wave = np.arange(wave_min, wave_max + step, step)

    # Lists for variable-length columns
    wave_list = []
    flux_list = []
    ivar_list = []
    orig_flux_list = []
    orig_ivar_list = []

    id_list = []
    red_list = []
    ra_list = []
    dec_list = []
    vdisp_list = []
    snr_list = []
    s_type_list = []

    dropped_indices = []

    for i, rec in enumerate(spec_records.records):
        # Resample spectrum to base_wave
        try:
            wave_i = np.array(rec.wavelength)
            flux_i = np.array(rec.flux)
            ivar_i = np.array(rec.ivar)
        except Exception:
            wave_i = np.array([])
            flux_i = np.array([])
            ivar_i = np.array([])

        if len(wave_i) == 0:
            # Treat as empty
            f_res = np.full(len(base_wave), np.nan, dtype=np.float64)
            iv_res = np.zeros(len(base_wave), dtype=np.float64)
        else:
            try:
                f_int = interp1d(wave_i, flux_i, bounds_error=False, fill_value=np.nan)
                iv_int = interp1d(wave_i, ivar_i, bounds_error=False, fill_value=0.0)
                f_res = f_int(base_wave)
                iv_res = iv_int(base_wave)
            except Exception:
                f_res = np.full(len(base_wave), np.nan, dtype=np.float64)
                iv_res = np.zeros(len(base_wave), dtype=np.float64)

        # Per-spectrum trimming to remove leading/trailing zeros/nans (seen in SDSS)
        w_trim, f_trim, iv_trim = per_spectrum_trim(base_wave, f_res, iv_res)

        # Decide keep/drop
        if f_trim.size < MIN_GOOD_PIXELS:
            # Drop
            dropped_indices.append(i)
            continue

        # Append trimmed arrays (variable-length)
        wave_list.append(w_trim)
        flux_list.append(f_trim)
        ivar_list.append(iv_trim)
        orig_flux_list.append(f_trim.copy())
        orig_ivar_list.append(iv_trim.copy())

        # Scalars
        sid = rec.get('specid') or rec.get('sparcl_id') or ''
        try:
            id_val = int(sid)
        except Exception:
            try:
                id_val = int(float(rec.get('sparcl_id', -1)))
            except Exception:
                id_val = -1
        id_list.append(id_val)
        red_list.append(float(rec.get('redshift') or np.nan))
        ra_list.append(float(rec.get('ra') or np.nan))
        dec_list.append(float(rec.get('dec') or np.nan))

        vdisp = np.nan
        tkey = rec.get('targetid') or rec.get('sparcl_id') or rec.get('specid')
        try:
            if tkey is not None:
                tid = int(tkey)
                if tid in fastspec_map:
                    vdisp = float(fastspec_map[tid])
        except Exception:
            pass
        if not np.isfinite(vdisp):
            R = rec.get('wave_sigma', None)
            try:
                vdisp = estimate_vdisp_from_resolution(base_wave, R)
            except Exception:
                vdisp = 125.0
        vdisp_list.append(float(vdisp))
        snr_list.append(compute_snr_from_flux_ivar(f_trim, iv_trim))
        s_type_list.append(str(rec.get('datasetgroup') or rec.get('data_release') or label))

    if len(id_list) == 0:
        raise RuntimeError(f"No {label} spectra left after trimming/dropping (MIN_GOOD_PIXELS={MIN_GOOD_PIXELS}).")

    # Build variable-length FITS columns using 'PD()' format (variable-length doubles)
    cols = []
    cols.append(fits.Column(name='ID', format='K', array=np.array(id_list, dtype=np.int64)))
    cols.append(fits.Column(name='WAVE', format='PD()', array=wave_list))
    cols.append(fits.Column(name='FLUX', format='PD()', array=flux_list))
    cols.append(fits.Column(name='IVAR', format='PD()', array=ivar_list))
    # original_data/original_ivar now identical to FLUX/IVAR trimmed arrays:
    cols.append(fits.Column(name='original_data', format='PD()', array=orig_flux_list))
    cols.append(fits.Column(name='original_ivar', format='PD()', array=orig_ivar_list))
    cols.append(fits.Column(name='REDSHIFT', format='D', array=np.array(red_list, dtype=np.float64)))
    cols.append(fits.Column(name='RA', format='D', array=np.array(ra_list, dtype=np.float64)))
    cols.append(fits.Column(name='DEC', format='D', array=np.array(dec_list, dtype=np.float64)))
    cols.append(fits.Column(name='VDISP', format='D', array=np.array(vdisp_list, dtype=np.float64)))
    cols.append(fits.Column(name='SNR', format='D', array=np.array(snr_list, dtype=np.float64)))
    cols.append(fits.Column(name='SURVEY_TYPE', format='15A', array=np.array([s.encode() for s in s_type_list])))

    hdu = fits.BinTableHDU.from_columns(cols, name=f'META_{label}')
    return hdu, len(id_list), dropped_indices

def main():
    print("Starting SPARCL retrieval with settings:")
    print(f" RA {RA_MIN} - {RA_MAX}   Dec {DEC_MIN} - {DEC_MAX}   z {Z_MIN} - {Z_MAX}")
    print(" Match radius (arcsec):", MATCH_RADIUS_ARCSEC, "  limit:", LIMIT)
    client = SparclClient()

    sdss_cons = {'spectype':['GALAXY'], 'ra':[RA_MIN,RA_MAX], 'dec':[DEC_MIN,DEC_MAX],
                 'specprimary':[True], 'datasetgroup':['SDSS_BOSS'], 'redshift':[Z_MIN,Z_MAX]}
    desi_cons = {'spectype':['GALAXY'], 'ra':[RA_MIN,RA_MAX], 'dec':[DEC_MIN,DEC_MAX],
                 'specprimary':[True], 'data_release':['DESI-DR1'], 'redshift':[Z_MIN,Z_MAX]}

    include_fields = ['wavelength','flux','ivar','sparcl_id','specid','ra','dec','redshift',
                      'datasetgroup','data_release','wave_sigma','wavemin','wavemax','targetid']

    print("Query SPARCL for SDSS...")
    rec_sdss = client.find(outfields=include_fields, constraints=sdss_cons, limit=LIMIT)
    print("Query SPARCL for DESI...")
    rec_desi = client.find(outfields=include_fields, constraints=desi_cons, limit=LIMIT)

    ra_sdss = np.array([r['ra'] for r in rec_sdss.records])
    dec_sdss = np.array([r['dec'] for r in rec_sdss.records])
    ra_desi = np.array([r['ra'] for r in rec_desi.records])
    dec_desi = np.array([r['dec'] for r in rec_desi.records])

    cat_desi = SkyCoord(ra=ra_desi*u.deg, dec=dec_desi*u.deg)
    cat_sdss = SkyCoord(ra=ra_sdss*u.deg, dec=dec_sdss*u.deg)

    print(f"Searching matches within {MATCH_RADIUS_ARCSEC} arcsec ...")
    ii_desi, ii_sdss, d2d, _ = search_around_sky(cat_desi, cat_sdss, MATCH_RADIUS_ARCSEC*u.arcsec)
    if len(ii_desi) == 0:
        print("No matches found. Exiting.")
        return

    pairs = sorted(list(zip(ii_desi, ii_sdss, d2d.arcsec)), key=lambda x:(x[1], x[2]))
    picked = {}
    for d,s,sep in pairs:
        if s not in picked:
            picked[s] = (d,s,sep)
    chosen = list(picked.values())
    Nmatch = len(chosen)
    print("Matched (unique SDSS -> DESI) count:", Nmatch)

    sdss_uuids = [str(rec_sdss.records[s]['sparcl_id']) for (_, s, _) in chosen]
    desi_uuids = [str(rec_desi.records[d]['sparcl_id']) for (d, _, _) in chosen]

    # Retrieve SDSS and DESI spectra
    try:
        spec_sdss = client.retrieve(uuid_list=sdss_uuids, include=include_fields, dataset_list=['BOSS-DR16','SDSS-DR16'])
    except Exception as e:
        print("Retrieval warning:", str(e))
        spec_sdss = client.retrieve(uuid_list=sdss_uuids, include=include_fields, dataset_list=None)

    spec_desi = client.retrieve(uuid_list=desi_uuids, include=include_fields, dataset_list=['DESI-DR1','DESI-EDR'])

    spec_sdss = spec_sdss.reorder(sdss_uuids)
    spec_desi = spec_desi.reorder(desi_uuids)

    fastspec_map = load_fastspec_vdisp(FASTSPEC_PATH) if FASTSPEC_PATH else {}

    sdss_waves = [np.array(r.wavelength) for r in spec_sdss.records]
    desi_waves = [np.array(r.wavelength) for r in spec_desi.records]

    sdss_identical = all_waves_equal(sdss_waves)
    desi_identical = all_waves_equal(desi_waves)
    if FORCE_RESAMPLE:
        sdss_identical = False
        desi_identical = False

    print("SDSS waves identical across sample?", sdss_identical)
    print("DESI waves identical across sample?", desi_identical)

    sdss_hdu, sdss_kept, sdss_dropped_list = build_varlen_hdu(spec_sdss, 'SDSS', sdss_waves, sdss_identical, fastspec_map)
    desi_hdu, desi_kept, desi_dropped_list = build_varlen_hdu(spec_desi, 'DESI', desi_waves, desi_identical, fastspec_map)

    out_sdss = os.path.join(OUT_DIR, f"{OUT_PREFIX}_sdss_dr16.fits")
    out_desi = os.path.join(OUT_DIR, f"{OUT_PREFIX}_desi_dr1.fits")

    fits.HDUList([fits.PrimaryHDU(), sdss_hdu]).writeto(out_sdss, overwrite=True)
    fits.HDUList([fits.PrimaryHDU(), desi_hdu]).writeto(out_desi, overwrite=True)

    print(f"Wrote: {out_sdss}  (kept {sdss_kept}, dropped {len(sdss_dropped_list)})")
    if sdss_dropped_list:
        print("Dropped SDSS indices (original ordering):", sdss_dropped_list)
    print(f"Wrote: {out_desi}  (kept {desi_kept}, dropped {len(desi_dropped_list)})")
    if desi_dropped_list:
        print("Dropped DESI indices (original ordering):", desi_dropped_list)
    print("Total matched objects initially:", Nmatch)

    # Quick verification
    for fn in (out_sdss, out_desi):
        with fits.open(fn, mode='readonly') as hd:
            print(f"\nInspecting: {fn}")
            hd.info()
            data = hd[1].data
            print("Columns:", list(data.names))
            try:
                w0 = data['WAVE'][0]
                print("First WAVE length:", len(w0))
                print("First FLUX length:", len(data['FLUX'][0]))
                print("Has original_data column:", 'original_data' in data.names)
                od = np.asarray(data['original_data'][0])
                fl = np.asarray(data['FLUX'][0])
                iv = np.asarray(data['IVAR'][0])
                print("original_data equals FLUX (first row)?", np.array_equal(od, fl))
                print("original_ivar equals IVAR (first row)?", np.array_equal(np.asarray(data['original_ivar'][0]), iv))
            except Exception as e:
                print("Verification error for file:", fn, e)

if __name__ == '__main__':
    main()