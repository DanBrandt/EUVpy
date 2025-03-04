from EUVpy.NEUVAC import neuvac
from EUVpy.tools import processIndices
import sys

if sys.version_info < (3, 9):
    from importlib_resources import files
else:
    from importlib.resources import files


# Testing different parts of EUVpy

if __name__ == '__main__':
    # Get some F10.7 for the entirety of 2018
    f107times, f107, f107a, f107b = processIndices.getCLSF107('2018-01-01', '2018-12-31', truncate=False)

    # Call the models in the EUVAC bins!
    neuvacIrr, _, _, _ = neuvac.neuvacEUV(f107, f107b, bands='EUVAC')


    sys.exit(0)