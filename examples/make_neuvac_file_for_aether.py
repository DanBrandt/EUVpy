# This script gives an example of how to convert the NEUVAC coefficients in the table file to a format for direct use
# by the Aether model.

# Local imports
from src.EUVpy.NEUVAC import neuvac

# Convert the table file:
out = neuvac.aetherFile()