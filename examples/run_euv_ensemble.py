# This script contains an example of how to run NEUVAC to generate an ensemble of irradiances and visualize the result.

# Top-level imports
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Qt5Agg')

# Local imports
from src.EUVpy.tools import processIndices
from src.EUVpy.tools import ensemble

# Global Plotting Settings:
import matplotlib.pylab as pylab
params = {'legend.fontsize': 'large',
          'figure.figsize': (16, 8),
         'axes.labelsize': 'large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

# Get some F10.7
f107times, f107, f107a, f107b = processIndices.getCLSF107('2014-08-01', '2014-10-01', truncate=False)

# Run the NEUVAC ensemble - do so with the version of NEUVAC in the EUVAC bins...
iterations = 50 # Number of ensemble members to generate
ensemble_NeuvacIrr, ensemble_average_NeuvacIrr, ensemble_stddev_NeuvacIrr = ensemble.irradiance_ensemble(f107, f107b,
                                                                                                iterations=iterations,
                                                                                                model='NEUVAC-E')

# Plot the ensemble, with the confidence bands corresponding to the ensemble spread (the 8th band this time)
ind = 7
plt.figure(figsize=(12,8))
plt.fill_between(f107times, (ensemble_average_NeuvacIrr-ensemble_stddev_NeuvacIrr)[:, ind],
                              (ensemble_average_NeuvacIrr+ensemble_stddev_NeuvacIrr)[:, ind],
                     color='orange', alpha=0.6)
plt.plot(f107times, ensemble_average_NeuvacIrr[:, ind], label='NEUVAC-37 (n='+str(iterations)+')', color='tab:orange', linewidth=5)
plt.xlabel('Date')
plt.ylabel('Irradiance W/m$^2$')
plt.legend(loc='best')
plt.tight_layout()
plt.show()