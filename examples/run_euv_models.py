# The script contains an example of how to run ALL the models contained with the solarEUV package and plot them
# for comparisons to each other.

# Top level imports
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Qt5Agg')

# Local imports
from src.EUVpy.tools import processIndices
from src.EUVpy.NEUVAC import neuvac
from src.EUVpy.empiricalModels.models.EUVAC import euvac
from src.EUVpy.empiricalModels.models.HEUVAC import heuvac
from src.EUVpy.empiricalModels.models.SOLOMON import solomon

# Get some F10.7 for the entirety of 2018
f107times, f107, f107a, f107b = processIndices.getCLSF107('2018-01-01', '2018-12-31', truncate=False)

# Call the models in the EUVAC bins!
neuvacIrr, _, _, _ = neuvac.neuvacEUV(f107, f107b, bands='EUVAC')
euvacFlux, euvacIrr, _, _, _ = euvac.euvac(f107, f107a)
heuvac_wav, heuvacFlux, heuvacIrr, _, _, _ = heuvac.heuvac(f107, f107a, torr=True)

# Call the models in the SOLOMON bins!
neuvac_tableFile_Stan_Bands = '../NEUVAC/neuvac_table_stan_bands.pkl'
neuvacIrrSolomon, _, _, _ = neuvac.neuvacEUV(f107, f107b, bands='SOLOMON')
solomonFluxHFG, solomonIrrHFG = solomon.solomon(f107, f107a, model='HFG')
solomonFluxEUVAC, solomonIrrEUVAC = solomon.solomon(f107, f107a, model='EUVAC')

# Plotting...
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
# Left plot will be total irradiance in the EUVAC bins
ax[0].plot(f107times, np.sum(neuvacIrr, axis=-1), label='NEUVAC-37', color='orange')
ax[0].plot(f107times, np.sum(euvacIrr, axis=-1), label='EUVAC-37', color='green')
ax[0].plot(f107times, np.sum(heuvacIrr, axis=-1), label='HEUVAC-37', color='red')
ax[0].set_ylabel('Irradiance (W/m$^2$)')
ax[0].set_xlabel('Date')
ax[0].set_title('Total Irradiance (EUVAC Bins)')
ax[0].legend(loc='best')
# Right subplot will be total irradiance in the SOLOMON bins
solomonTable = solomon.solomonTable
mids = 0.5*(solomonTable[:, 1] + solomonTable[:, 2])*10
ind2 = 7
ax[1].plot(f107times, np.sum(neuvacIrrSolomon, axis=-1), label='NEUVAC-22', color='orange')
ax[1].plot(f107times, np.sum(solomonIrrHFG, axis=-1), label='HFG', color='purple')
ax[1].plot(f107times, np.sum(solomonIrrEUVAC, axis=-1), label='HEUVAC-22', color='red')
ax[1].set_ylabel('Irradiance (W/m$^2$)')
ax[1].set_xlabel('Date')
ax[1].set_title('Total Irradiance (SOLOMON bins)')
ax[1].legend(loc='best')
fig.tight_layout()
plt.show()