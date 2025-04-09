# Contains draft files for creating unit tests.
from EUVpy.NEUVAC import neuvac
from EUVpy.tools import processIndices
import sys

if sys.version_info < (3, 9):
    from importlib_resources import files
else:
    from importlib.resources import files

# Testing different parts of EUVpy

if __name__ == '__main__':
    # Get some F10.7 for the entirety of 2018:
    f107times, f107, f107a, f107b = processIndices.getCLSF107('2018-01-01', '2018-12-31', truncate=False)

    # Call the model in the EUVAC bins:
    neuvacIrr, _, _, _ = neuvac.neuvacEUV(f107, f107b, bands='EUVAC')
    # Integrals:
    from scipy.integrate import simpson
    all_ints = []
    for i in range(neuvacIrr.shape[1]):
        all_ints.append(simpson(neuvacIrr[:, i]))

    # Call the model in the GITM bins and compute integrals:
    neuvacIrrG, _, _, _ = neuvac.neuvacEUV(f107, f107b)
    all_ints_g = []
    for j in range(neuvacIrrG.shape[1]):
        all_ints_g.append(simpson(neuvacIrrG[:, j]))

    # Call the model in the SOLOMON bands and compute integrals:
    neuvacIrrS, _, _, _ = neuvac.neuvacEUV(f107, f107b, bands='SOLOMON')
    all_ints_s = []
    for k in range(neuvacIrrS.shape[1]):
        all_ints_s.append(simpson(neuvacIrrS[:, k]))

    sys.exit(0)