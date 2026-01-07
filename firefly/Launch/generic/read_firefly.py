"""
===================================================================================
|                         FIREFLY — Full Spectral Fitting                         |
===================================================================================
 <> Module: read_firefly.py

 <> Author:
	- Daniel Thomas ~ <d.g.thomas__at__leeds.ac.uk>

 <> Contributors:
 	- Samuel Helps ~ <samuel.helps.sh__at__gmail.com>

 <> Purpose:
	- Read in and plot FIREFLY output files (most basic version).
___________________________________________________________________________________
|Institute of Cosmology and Gravitation - University of Portsmouth, Portsmouth, UK|
‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
"""
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import sys

# Read in FIREFLY output file (alter as needed using sys.argv or hard-coded paths)
hdul = fits.open(sys.argv[1])  # hdul = fits.open("spFly-spec-0266-51602-0001.fits")
data=hdul[1].data
wave = data['wavelength']
flux = data['original_data']
model = data['firefly_model']

# Print out some of the header information
hdul.close()
hdul.info()

# Extract parameters of the best-fit
csp_age=np.ndarray(hdul[1].header['ssp_number'])
csp_Z=np.ndarray(hdul[1].header['ssp_number'])
csp_light=np.ndarray(hdul[1].header['ssp_number'])
csp_mass=np.ndarray(hdul[1].header['ssp_number'])
for i in range(len(csp_age)):
	csp_age[i]=hdul[1].header['log_age_ssp_'+str(i)]
	csp_Z[i]=hdul[1].header['metal_ssp_'+str(i)]
	csp_light[i]=hdul[1].header['weightLight_ssp_'+str(i)]
	csp_mass[i]=hdul[1].header['weightMass_ssp_'+str(i)]

# Print out some results
print()
print(hdul[0].header)
print(hdul[1].header)
print()
print('age: '+str(np.around(10**hdul[1].header['age_lightW'],decimals=2))+' Gyr')
print('[Z/H]: '+str(np.around(hdul[1].header['metallicity_lightW'],decimals=2))+' dex')
print('log M/Msun: '+str(np.around(hdul[1].header['stellar_mass'],decimals=2)))
print('E(B-V): '+str(np.around(hdul[1].header['EBV'],decimals=2))+' mag')

# Plot the orginal spectrum and the best-fit model
plt.plot(wave,flux)
plt.plot(wave,model)
plt.show()

# Plot the star formation history
fig1=plt.figure()
plt.xlim(0,15)
plt.xlabel('lookback time (Gyr)')
plt.ylabel('frequency')
#plt.bar(10**(csp_age),csp_light,width=1,align='center',edgecolor='k',linewidth=2)
plt.bar(10**(csp_age),csp_light,width=1,align='center',alpha=0.5)
plt.scatter(10**(csp_age),csp_light)
plt.show()

# Plot the metallicity distribution
fig2=plt.figure()
plt.xlim(-2,0.5)
plt.xlabel('[Z/H] (dex)')
plt.ylabel('frequency')
#plt.bar(10**(csp_age),csp_light,width=1,align='center',edgecolor='k',linewidth=2)
plt.bar(csp_Z,csp_light,width=0.1,align='center',alpha=0.5)
plt.scatter(csp_Z,csp_light)
plt.show()