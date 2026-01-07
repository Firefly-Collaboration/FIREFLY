<div align="center">

```
███████╗██╗██████╗ ███████╗███████╗██╗  ██╗   ██╗
██╔════╝██║██╔══██╗██╔════╝██╔════╝██║  ╚██╗ ██╔╝
█████╗  ██║██████╔╝█████╗  █████╗  ██║   ╚████╔╝ 
██╔══╝  ██║██╔══██╗██╔══╝  ██╔══╝  ██║    ╚██╔╝  
██║     ██║██║  ██║███████╗██║     ███████╗██║   
╚═╝     ╚═╝╚═╝  ╚═╝╚══════╝╚═╝     ╚══════╝╚═╝   
```

# FIREFLY - Full Spectral Fitting

### **F**itting **I**te**R**ativ**E**ly **F**or **L**ikelihood Anal**Y**sis

**[📄 Cite](CITATION.cff)** • **[🌐 Website](https://firefly-collaboration.github.io/FIREFLY/)** • **[🔗 Old Release](https://github.com/FireflySpectra/firefly_release)** • **[🚀 Get Started](https://github.com/Firefly-Collaboration/FIREFLY/tree/main/firefly/Launch)**

---

</div>

<div align="center">

![Spectral example](docs/images/EDR_2814R_bigpng.png)

</div>

---

<div align="center">

## 📚 Documentation • Tutorials • Examples

[![🌐 Visit Our Website](https://img.shields.io/badge/🌐_Visit_Our_Website-Full_Documentation-orange?style=for-the-badge)](https://firefly-collaboration.github.io/FIREFLY/)

**For the full styled experience with interactive features, visit our website above**

</div>

---

<div align="center">

![Spectral example](docs/images/BGS_329.png)

![Spectral example](docs/images/BGS_4748.png)

</div>

---

## 📖 Table of Contents

<table width="100%">
<tr>
<th width="50%">Core Documentation</th>
<th width="50%">Advanced Topics</th>
</tr>
<tr>
<td>

• [Overview](#-overview)  
• [Key Features](#-key-features)  
• [Installation](#-installation)  
• [Quick Start](#-quick-start)  
• [Usage](#-usage)  
• [Repository Structure](#-repository-structure)  
• [Stellar Population Models](#-stellar-population-models)

</td>
<td>

• [Output Products](#-output-products)  
• [Value Added Catalogs](#-value-added-catalogs)  
• [Citation](#-citation)  
• [Contributing](#-contributing)  
• [FIREFLY Collaboration](#-firefly-collaboration)  
• [License](#-license)  
• [Original Release](#-original-release)

</td>
</tr>
</table>

---

<div align="center">

## 🔭 Overview

</div>

> **FIREFLY is a chi-squared minimization fitting code designed for deriving the stellar population properties of stellar systems from spectroscopic data. Whether analysing the very lastest observed galaxy spectra, samples from huge surveys, or model spectra from simulations, FIREFLY provides a prior-free fitting approach.**

<table>
<tr>
<td width="50%" valign="top">

### ⚙️ **How It Works**

- **FIREFLY fits combinations of single-burst stellar population (SSP) models** to spectroscopic data following an iterative best-fitting process controlled by the Bayesian Information Criterion (BIC).

- **No imposed priors:** All solutions within a statistical cut are retained with their weights.

- **No additive or multiplicative polynomials:** Spectral shape is not artificially adjusted.

- **No regularisation:** Maximum fitting freedom to map intrinsic SED degeneracies, such as age, metallicity, dust reddening on stellar population properties.

- **High-Pass Filter (HPF) dust treatment:** Novel procedure for continuum rectification and dust attenuation. The returned attenuation array is then matched to known analytical approximations to return an E(B-V) value. This procedure allows for removal of large scale modes of the spectrum associated with dust and/or poor flux calibration.

</td>
<td width="50%" valign="top">

### 🎯 **Performance**

- **Comprehensive VACs produced** from a variety of telescopes and different data sources.

- **Extensivley tested** on real galaxy spectra from the Sloan Digital Sky Survey (SDSS).

- **Applied to spectra** from SDSS-IV/MaNGA integral field spectroscopy to analyse millions of galaxies.

- **Upgraded to derive** the stellar population properties of galaxies from the latest Dark Energy Spectroscopic Instrument (DESI) observations.

- **Tested on data** from the DEEP2 survey and Milky Way globular clusters.

- **Robust recovery** down to S/N ~ 5 for moderately dusty systems.

</td>
</tr>
</table>

---

<div align="center">

![FIREFLY software workflow](docs/images/firefly_workflow.png)

</div>

---

<div align="center">

## ✨ Key Features

</div>

- **Flexible Input:** FIREFLY can be tailored to fit spectra from a variety of different data formats, spectral resolutions and wavelength ranges.
- **Comprehensive Output:** Provides both light and mass weighted stellar population properties including age, metallicity, dust attenuation E(B-V), stellar mass and remnant mass partition (white dwarfs, neutron stars, black holes), star formation rates and histories, and SSP-specific component weights.
- **Multiple Model Libraries:** Support for MaStar (Maraston et al. 2020) and M11 (Maraston & Strömback 2011) stellar population models.
- **Survey Integration:** Dedicated pipelines for DESI, SDSS, and MaNGA data sources.
- **Emission Line Masking:** Configurable masking settings for accurate continuum fitting.
- **Multiple Viable Solutions:** Using the Bayesian Information Criterion (BIC), FIREFLY retains several fits rather than a single best-fit, improving convergence with models and allowing for realistic error estimates on derived parameters.

---

<div align="center">

## 🚀 Installation

</div>

### **Requirements**

```yaml
Python: 3.6+
Core Packages: numpy, astropy, matplotlib
Optional Packages: argparse, multiprocessing, concurrent.futures, subprocess, fcnt, pandas
```

---

### **Installation Steps**

#### **1️⃣ Clone the repository:**

```bash
git clone https://github.com/Firefly-Collaboration/FIREFLY.git
cd firefly
```

#### **2️⃣ Install Python dependencies:**

```bash
pip install -r requirements.txt
```

#### **3️⃣ Environment Variables (Optional):**

<table>
<tr>
<td width="50%">

**For bash/zsh (.bashrc or .bash_profile):**

```bash
export FF_DIR='/path/to/FIREFLY'
export PYTHONPATH="${FF_DIR}/firefly/Fitting_Engine:$PYTHONPATH"
export STELLARPOPMODELS_DIR="${FF_DIR}/firefly/Fitting_Engine/stellar_population_models"
```

</td>
<td width="50%">

</td>
</tr>
</table>

---

### **Stellar Population Models**

> Stellar population model templates are included in the repository: M11 Models (Maraston & Strömback 2011): MILES, STELIB, ELODIE, MARCS libraries. MaStar Models: High-resolution empirical stellar library. IMF Options: Kroupa and Salpeter initial mass functions. Models are located in: `firefly/Fitting_Engine/stellar_population_models/` and new models can be added to the module providing they are formatted for FIREFLY compatability.

---

<div align="center">

## ⚡ Quick Start

</div>

```bash
cd firefly/Launch/generic
python firefly.py
```

> **This will fit the example spectrum located in `firefly/Data/example_data/`.**

---

<div align="center">

## 📊 Usage

</div>

<table>
<tr>
<td width="50%" valign="top">

### 📄 **Generic ASCII Input**

1. Navigate: `cd firefly/Launch/generic`
2. Edit `firefly.py` to set input_file, redshift, and parameters.
3. Run: `python firefly.py`
4. Read output: `python read_firefly.py path/to/output.fits`

**Input Format:**
> ASCII with wavelength, flux, and error columns.

</td>
<td width="50%" valign="top">

### 🌌 **SDSS Spectra**

```bash
cd firefly/Launch/SDSS
python firefly_SDSS.py
```

> Edit the script to point to new SDSS spec files; redshift and metadata are read from FITS headers.

</td>
</tr>
<tr>
<td width="50%" valign="top">

### 🔬 **MaNGA Data Cubes**

```bash
cd firefly/Launch/MANGA
python firefly_MANGA.py
```

> Processes MaNGA data cubes and fits Voronoi binned spectra. Configure paths to logcube and DAP files inside the script.

</td>
<td width="50%" valign="top">

### 🌠 **DESI Spectra**

**NERSC (DESI-DR1):**

The quickest version of FIREFLY, firefly(AIO), can be run on NERSC in the DESI-DR1 example pipleine.

```bash
cd FIREFLY/Launch/DESI/NERSC/run_scripts
sbatch SBATCH_Iron.sh
```

**SCIAMA HPC (DESI-EDR):**

For users who don't have access to NERSC a DESI-EDR example pipeline uses the SCIAMA HPC (Insitute of Cosmology and Gravitation, Portsmouth, UK) but can be easily updated for use on another HPC.

```bash
cd FIREFLY/Launch/DESI/SCIAMA/run_scripts/
sbatch SBATCH_Fuji.sh DESI_EDR_10000-20000.fits
```

</td>
</tr>
</table>

---

### ⚙️ **Configuration Options**

<table width="100%">
<tr>
<th width="40%">Parameter</th>
<th width="60%">Options</th>
</tr>
<tr>
<td><strong>Stellar Population Models</strong></td>
<td><code>'m11'</code> or <code>'MaStar'</code></td>
</tr>
<tr>
<td><strong>Model Library</strong></td>
<td><code>'MILES'</code>, <code>'STELIB'</code>, <code>'ELODIE'</code>, <code>'MARCS'</code>, <code>'gold'</code></td>
</tr>
<tr>
<td><strong>Initial Mass Function</strong></td>
<td><code>'kr'</code> (Kroupa) or <code>'ss'</code> (Salpeter)</td>
</tr>
<tr>
<td><strong>Wavelength Range</strong></td>
<td>Customise fitting limits</td>
</tr>
<tr>
<td><strong>Emission Line Masking</strong></td>
<td>Enable/disable and specify lines to mask</td>
</tr>
<tr>
<td><strong>Dust Treatment</strong></td>
<td>Configure High-Pass Filter parameters</td>
</tr>
</table>

> **💡 See the example scripts in `firefly/Launch/` for detailed configuration examples.**

---

<div align="center">

## 📁 Repository Structure

</div>

![FIREFLY Repository Structure](docs/images/RepoStructure.jpg)

---

<div align="center">

## 🌟 Stellar Population Models

</div>

### **MaStar Models**

- **v0.2:** Initial release
- **v1.1:** Updated version with improved calibration
- **gold:** Highest quality stellar templates

> MaStar: Maraston et al. 2020 — [ADS](https://ui.adsabs.harvard.edu/abs/2020MNRAS.496.2962M) | [BibTeX](https://ui.adsabs.harvard.edu/abs/2020MNRAS.496.2962M/exportcitation)

---

### **M11 Models**

- **MILES:** Medium resolution (FWHM ~ 2.5Å)
- **STELIB:** Empirical stellar library
- **ELODIE:** High-resolution stellar library
- **MARCS:** Theoretical stellar atmospheres (Kroupa IMF only)

> M11: Maraston & Strömbäck 2011 — [ADS](https://ui.adsabs.harvard.edu/abs/2011MNRAS.418.2785M) | [BibTeX](https://ui.adsabs.harvard.edu/abs/2011MNRAS.418.2785M/exportcitation)

---

### **Available Initial Mass Functions**

- **Kroupa IMF (`'kr'`):** Broken power law IMF (Kroupa 2001)
- **Salpeter IMF (`'ss'`):** Single power law IMF (Salpeter 1955)

---

<div align="center">

## 📤 Output Products

</div>

### **Header Information**

| **Parameter** | **Description** |
|:-------------|:---------------|
| `stellar_mass` | log(M*/M☉) |
| `age_lightW` | Light-weighted log(age/yr) |
| `age_massW` | Mass-weighted log(age/yr) |
| `metallicity_lightW` | Light-weighted [Z/H] |
| `metallicity_massW` | Light-weighted [Z/H] |
| `EBV` | Dust attenuation E(B-V) in magnitudes |
| `ssp_number` | Number of SSP components in best fit |
| Individual SSP properties | `log_age_ssp_X`, `metal_ssp_X`, `weightLight_ssp_X`, `weightMass_ssp_X` |

---

### **Reading Output**

#### **Load FIREFLY output:**

```python
from astropy.io import fits
import numpy as np

hdul = fits.open('output_file.fits')
data = hdul[1].data
```

---

#### **Extract best-fit parameters:**

```python
wave  = data['wavelength']
flux  = data['original_data']
model = data['firefly_model']

stellar_mass = hdul[1].header['stellar_mass']
age_lightW = hdul[1].header['age_lightW']
metallicity = hdul[1].header['metallicity_lightW']
ebv = hdul[1].header['EBV']
```

---

#### **Display results:**

```python
hdul.info()
hdul.close()

print('Age: ' + str(np.around(10**age_lightW, decimals=2)) + ' Gyr')
print('[Z/H]: ' + str(np.around(metallicity, decimals=2)) + ' dex')
print('log M/M☉: ' + str(np.around(stellar_mass, decimals=2)))
print('E(B-V): ' + str(np.around(ebv, decimals=2)) + ' mag')
```

---

#### **Plot spectrum and best-fit model:**

```python
import matplotlib.pyplot as plt

plt.plot(wave, flux, label='Observed Spectrum')
plt.plot(wave, model, label='FIREFLY Model')
plt.xlabel('Wavelength')
plt.ylabel('Flux')
plt.legend()
plt.show()
```

---

<div align="center">

## 📊 Value Added Catalogs

</div>

> The `Galaxy_VACs/` directory is designed to store store just some Value Added Catalogs (VACs) created with FIREFLY for large galaxy surveys. Presented in various formats, these catalogs contain stellar population properties of millions of galaxies derived from a variety of spectroscopic surveys.

> **📝 Note:** The VACs are typically very large files and may be hosted separately.

---

<div align="center">

## 📝 Citation

</div>

> **If you use FIREFLY or its resources for work/research presented in a publication we ask that you please cite the following papers:**

<table width="100%">
<tr>
<th width="50%">Publication</th>
<th width="50%">Links</th>
</tr>
<tr>
<td><strong>FIREFLY:</strong> Wilkinson et al. 2017</td>
<td><a href="https://ui.adsabs.harvard.edu/abs/2017MNRAS.472.4297W">ADS</a> • <a href="https://ui.adsabs.harvard.edu/abs/2017MNRAS.472.4297W/exportcitation">BibTeX</a></td>
</tr>
<tr>
<td><strong>FIREFLY (ASCL):</strong> Wilkinson et al. 2021</td>
<td><a href="https://ui.adsabs.harvard.edu/abs/2021ascl.soft08010W">ADS</a> • <a href="https://ui.adsabs.harvard.edu/abs/2021ascl.soft08010W/exportcitation">BibTeX</a></td>
</tr>
<tr>
<td><strong>FIREFLY:</strong> Neumann et al. 2022</td>
<td><a href="https://ui.adsabs.harvard.edu/abs/2022MNRAS.513.5988N">ADS</a> • <a href="https://ui.adsabs.harvard.edu/abs/2022MNRAS.513.5988N/exportcitation">BibTeX</a></td>
</tr>
<tr>
<td><strong>MaStar:</strong> Maraston et al. 2020</td>
<td><a href="https://ui.adsabs.harvard.edu/abs/2020MNRAS.496.2962M">ADS</a> • <a href="https://ui.adsabs.harvard.edu/abs/2020MNRAS.496.2962M/exportcitation">BibTeX</a></td>
</tr>
<tr>
<td><strong>M11:</strong> Maraston & Strömbäck 2011</td>
<td><a href="https://ui.adsabs.harvard.edu/abs/2011MNRAS.418.2785M">ADS</a> • <a href="https://ui.adsabs.harvard.edu/abs/2011MNRAS.418.2785M/exportcitation">BibTeX</a></td>
</tr>
</table>

**📄 A BibTeX file is provided in [CITATION.cff](CITATION.cff) for your convenience.**

---

<div align="center">

## 🤝 Contributing

</div>

> **We welcome contributions to FIREFLY! Whether you're fixing bugs, adding new features, or improving documentation, your help is appreciated. If you fit spectra from a new survey with FIREFLY, adding the new pipeline to the project with help future users.**

### **How to Contribute**

```bash
# 1. Fork the repository

# 2. Create a feature branch
git checkout -b feature/YourFeature

# 3. Commit your changes
git commit -m 'Add YourFeature'

# 4. Push to the branch
git push origin feature/YourFeature

# 5. Open a Pull Request
```

> Please ensure your code follows the existing style and includes appropriate documentation.

---

<div align="center">

## 👥 FIREFLY Collaboration

</div>

<table width="100%">
<tr>
<td colspan="2" align="center">

### **Original Module Author**
**David M. Wilkinson** - Core fitting engine development

</td>
</tr>
<tr>
<td width="50%" align="center">

### **Principal Investigators**
**Daniel Thomas**  
**Claudia Maraston**

</td>
<td width="50%" align="center">

### **Website & Repository Developer**
**Samuel Helps**

</td>
</tr>
</table>

### **Module Contributors**

<table width="100%">
<tr>
<th width="40%">Contributor</th>
<th width="60%">Role</th>
</tr>
<tr>
<td><strong>Daniel Thomas</strong></td>
<td>Core scripting and SDSS/generic pipelines</td>
</tr>
<tr>
<td><strong>Johan Comparat</strong></td>
<td>Spectral setup utilities and SDSS pipeline</td>
</tr>
<tr>
<td><strong>Justus Neumann</strong></td>
<td>MaNGA pipeline and fitting engine support</td>
</tr>
<tr>
<td><strong>Violeta Gonzalez-Perez</strong></td>
<td>Utilities, dust and models support</td>
</tr>
<tr>
<td><strong>Daniel Goddard</strong></td>
<td>Instrument and MaNGA pipeline support</td>
</tr>
<tr>
<td><strong>Sofia Meneses-Goytia</strong></td>
<td>Estimations, setup and models support</td>
</tr>
<tr>
<td><strong>Samuel Helps</strong></td>
<td>DESI-EDR pipeline and repository development</td>
</tr>
<tr>
<td><strong>Kieran Graham</strong></td>
<td>DESI DR1 (AIO) pipeline</td>
</tr>
<tr>
<td><strong>Harry Hicks</strong></td>
<td>Models module support</td>
</tr>
<tr>
<td><strong>Kyle Westfall</strong></td>
<td>Utilities and MaNGA DAP (constants) integration</td>
</tr>
</table>

### **Institution**

> **Institute of Cosmology and Gravitation**  
> University of Portsmouth  
> Portsmouth, United Kingdom

---

<div align="center">

## 📜 License

</div>

> This project is licensed under the [**MIT License**](https://opensource.org/licenses/MIT) - see the [LICENSE.md](LICENSE.md) file for details.
>
> Unless otherwise specified, media content is licensed under [**CC BY 4.0**](https://creativecommons.org/licenses/by/4.0/).

---

<div align="center">

## 🔗 Original Release

</div>

> This new release of FIREFLY is based on the original repository. To explore the old github, visit: 
>
> **https://github.com/FireflySpectra/firefly_release**

---

<div align="center">

### ⚠️ **Disclaimer**

> **FIREFLY is provided as-is for academic purposes. The FIREFLY Collaboration assumes no liability for misuse of this software, website, media content, or any outputs generated by this code. While extensively tested, users are responsible for validating results for their specific applications.**

---

</div>