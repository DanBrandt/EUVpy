# This script gives an example of how to convert the NEUVAC coefficients in the table file to a format for direct use
# by the Aether model.

# Local imports
from NEUVAC import neuvac

# Convert the table file:
tableFilename = '../NEUVAC/neuvac_table.txt'
out = neuvac.aetherFile(tableFilename)