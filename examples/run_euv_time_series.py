# This example focuses on creating simple irradiance time series for NEUVAC.

# Top level imports
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Qt5Agg')
from tools import processIndices
# import numpy as np

# Local imports
from NEUVAC import neuvac
from empiricalModels.models.EUVAC import euvac
from empiricalModels.models.HEUVAC import heuvac

# Global Plotting Settings:
import matplotlib.pylab as pylab
params = {'legend.fontsize': 'large',
          'figure.figsize': (16, 8),
         'axes.labelsize': 'large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

# Get some F10.7 data:
# f107 = np.random.uniform(low=60, high=200, size=(100,))
# f107a = f107
f107times, f107, f107a, f107b = processIndices.getCLSF107('2008-12-01', '2019-12-31', truncate=False)
# Above, note that 'truncate' should be set to True only when you are looking at more than 3 solar rotations of F10.7
# data, and that data INCLUDES data up until the present time.

# neuvac_tableFile = '../NEUVAC/neuvac_table.txt'
neuvacIrradiance, _, _, _ = neuvac.neuvacEUV(f107, f107a)#, tableFile=neuvac_tableFile)
euvacFlux, euvacIrr, _, _, _ = euvac.euvac(f107, f107a)
# heuvac_wav, heuvacFlux, heuvacIrr, _, _, _ = heuvac.heuvac(f107, f107a, torr=True)

# Just plot the irradiance in the 10th wavelength band:
ind = 11
mids = 0.5*(euvac.euvacTable[:, 1] + euvac.euvacTable[:, 2] )
plt.figure(figsize=(12, 8))
plt.plot(f107times, neuvacIrradiance[:, ind], label='NEUVAC', color='tab:orange')
plt.plot(f107times, euvacIrr[:, ind], label='EUVAC', color='tab:green')
# plt.plot(f107times, heuvacIrr[:, ind], label='HEUVAC', color='tab:red')
plt.legend(loc='best')
plt.xlabel('Date')
plt.ylabel('Irradiance W/m$^2$')
plt.title('Solar Irradiance at '+str(mids[ind])+' Angstroms (SC24)')
plt.tight_layout()
plt.show()

