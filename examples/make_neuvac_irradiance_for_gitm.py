# This script gives an example of how to generate a set of NEUVAC irradiances for ingestion directly into the
# Global Ionosphere Thermosphere Model

# Local imports
from src.EUVpy.tools import processIndices
from src.EUVpy.NEUVAC import neuvac

# Get some F10.7 for the entirety of 2022
f107times, f107, f107a, f107b = processIndices.getCLSF107('2022-01-01', '2023-01-01', truncate=False)

# Generate the NEUVAC irradiance but output them to a file that can be used by GITM
out = neuvac.gitmNEUVAC(f107times, f107, f107b)

# It's as simple as that!
