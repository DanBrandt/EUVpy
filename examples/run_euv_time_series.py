# This example focuses on creating simple irradiance time series for NEUVAC.

# Top level imports
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Qt5Agg')
from tools import processIndices
# import numpy as np

# Local imports
from NEUVAC import neuvac

# Get some F10.7 data:
# f107 = np.random.uniform(low=60, high=200, size=(100,))
# f107a = f107
f107times, f107, f107a, f107b = processIndices.getCLSF107('2008-12-01', '2019-12-31', truncate=False)
# Above, note that 'truncate' should be set to True only when you are looking at more than 3 solar rotations of F10.7
# data, and that data INCLUDES data up until the present time.

# neuvac_tableFile = '../NEUVAC/neuvac_table.txt'
neuvacIrradiance, _, _, _ = neuvac.neuvacEUV(f107, f107a)#, tableFile=neuvac_tableFile)

# Just plot the irradiance in the 10th wavelength band:
plt.figure(figsize=(12, 8))
plt.plot(f107times, neuvacIrradiance[:, 9])
plt.xlabel('Date')
plt.ylabel('Irradiance W/m$^2$')
plt.tight_layout()
plt.show()

