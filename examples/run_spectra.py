# This example focuses on creating individual spectra from different models:

# Top level imports
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Qt5Agg')
import numpy as np

# Local imports
from src.EUVpy.tools import toolbox
from src.EUVpy.NEUVAC import neuvac
from src.EUVpy.empiricalModels.models.EUVAC import euvac
from src.EUVpy.empiricalModels.models.HEUVAC import heuvac
from src.EUVpy.empiricalModels.models.SOLOMON import solomon

# Global Plotting Settings:
import matplotlib.pylab as pylab
params = {'legend.fontsize': 'large',
          'figure.figsize': (16, 8),
         'axes.labelsize': 'large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

# Sample values for F10.7, F10.7A, and F10.7B
f107 = 120
f107a = 85
f107b = 87

# In this example, we'll only plot the wavelength ranges, not the individual lines. In order to do that, we need to get
# the boundaries of the wavelength ranges to make a nice stair plot (think Fig. 8 from Nishimoto, et al:
# https://link.springer.com/article/10.1186/s40623-021-01402-7).

euvacTable = euvac.euvacTable
leftsides = euvacTable[:, 1]
rightsides = euvacTable[:, 2]
band_indices, band_boundaries = toolbox.band_info(leftsides, rightsides)

# In the EUVAC bins:
neuvacIrr, _, _, _ = neuvac.neuvacEUV(f107, f107b, bands='EUVAC')
euvacFlux, euvacIrr, _, _, _ = euvac.euvac(f107, f107a)
heuvac_wav, heuvacFlux, heuvacIrr, _, _, _ = heuvac.heuvac(f107, f107a, torr=True)

fig, ax = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True)
ax.stairs(values=neuvacIrr[0, band_indices], edges=band_boundaries, label='NEUVAC-37', lw=3, color='tab:orange')
ax.stairs(values=euvacIrr[0, band_indices], edges=band_boundaries, label='EUVAC-37', lw=3, color='tab:green')
ax.stairs(values=heuvacIrr[0, band_indices], edges=band_boundaries, label='HEUVAC-37', lw=3, color='tab:red')
ax.set_yscale('log')
ax.legend(loc='best')
ax.grid()
ax.set_xlabel('Wavelength ($\mathrm{\AA}$)')
ax.set_ylabel('Irradiance (W/m$^2$)')
ax.set_title('Individual Solar Spectra in EUVAC Bins for (F10.7, F10.7A, F10.7B) = ('+str(f107)+', '+str(f107a)+', '+str(f107b)+')')
plt.show()

# In the SOLOMON bins:
solomonTable = solomon.solomonBands
leftsides_solomon = solomonTable['short']
rightsides_solomon = solomonTable['long']
band_indices_solomon, band_boundaries_solomon = toolbox.band_info(leftsides_solomon, rightsides_solomon, solomon=True)

neuvac_tableFile_Stan_Bands = '../NEUVAC/neuvac_table_stan_bands.pkl'
neuvacIrrSolomon, _, _, _ = neuvac.neuvacEUV(f107, f107b, bands='SOLOMON')
solomonFluxHFG, solomonIrrHFG = solomon.solomon(f107, f107a, model='HFG')
solomonFluxEUVAC, solomonIrrEUVAC = solomon.solomon(f107, f107a, model='EUVAC')

neuvacIrrSolomon_for_plotting = toolbox.solomonRebin(neuvacIrrSolomon)
hfgIrrSolomon_for_plotting = toolbox.solomonRebin(solomonIrrHFG)
euvacIrrSolomon_for_plotting = toolbox.solomonRebin(solomonIrrEUVAC)

fig2, ax = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True)
ax.stairs(values=neuvacIrrSolomon_for_plotting[0, band_indices_solomon], edges=band_boundaries_solomon, label='NEUVAC-22', lw=3, color='tab:orange')
ax.stairs(values=hfgIrrSolomon_for_plotting[0, band_indices_solomon], edges=band_boundaries_solomon, label='HFG', lw=3, color='purple')
ax.stairs(values=euvacIrrSolomon_for_plotting[0, band_indices_solomon], edges=band_boundaries_solomon, label='EUVAC-22', lw=3, color='tab:green')
ax.set_yscale('log')
ax.legend(loc='best')
ax.grid()
ax.set_xlabel('Wavelength ($\mathrm{\AA}$)')
ax.set_ylabel('Irradiance (W/m$^2$)')
ax.set_title('Individual Solar Spectra in SOLOMON Bins for (F10.7, F10.7A, F10.7B) = ('+str(f107)+', '+str(f107a)+', '+str(f107b)+')')
plt.show()
