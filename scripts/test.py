# Contains draft files for creating unit tests.
from EUVpy.NEUVAC import neuvac
from EUVpy.tools import processIndices
from EUVpy.empiricalModels.models.EUVAC import euvac
from EUVpy.empiricalModels.models.SOLOMON import solomon
from EUVpy.empiricalModels.models.HEUVAC import heuvac
import sys
import pathlib

#-----------------------------------------------------------------------------------------------------------------------
# Directory management:
here = pathlib.Path(__file__).parent.resolve()
#-----------------------------------------------------------------------------------------------------------------------

if sys.version_info < (3, 9):
    from importlib_resources import files
else:
    from importlib.resources import files

# Testing different parts of EUVpy

if __name__ == '__main__':
    # Get some F10.7 for the entirety of 2018:
    f107times, f107, f107a, f107b = processIndices.getCLSF107('2018-01-01', '2018-12-31', truncate=False)

    # # Call the model in the EUVAC bins:
    # neuvacIrr, _, _, _ = neuvac.neuvacEUV(f107, f107b, bands='EUVAC')
    # # Integrals:
    # from scipy.integrate import simpson
    # all_ints = []
    # for i in range(neuvacIrr.shape[1]):
    #     all_ints.append(simpson(neuvacIrr[:, i]))
    #
    # # Call the model in the GITM bins and compute integrals:
    # neuvacIrrG, _, _, _ = neuvac.neuvacEUV(f107, f107b)
    # all_ints_g = []
    # for j in range(neuvacIrrG.shape[1]):
    #     all_ints_g.append(simpson(neuvacIrrG[:, j]))
    #
    # # Call the model in the SOLOMON bands and compute integrals:
    # neuvacIrrS, _, _, _ = neuvac.neuvacEUV(f107, f107b, bands='SOLOMON')
    # all_ints_s = []
    # for k in range(neuvacIrrS.shape[1]):
    #     all_ints_s.append(simpson(neuvacIrrS[:, k]))
    #
    # # Call the EUVAC model:
    # euvacFlux, euvacIrr, _, _, _ = euvac.euvac(f107, f107a)
    # all_ints_e = []
    # for l in range(euvacIrr.shape[1]):
    #     all_ints_e.append(simpson(euvacIrr[:, l]))
    #
    # # Call the HFG model:
    # solomonFluxHFG, solomonIrrHFG = solomon.solomon(f107, f107a, model='HFG')
    # all_ints_hfg = []
    # for m in range(solomonIrrHFG.shape[1]):
    #     all_ints_hfg.append(simpson(solomonIrrHFG[:, m]))
    #
    # # Call the EUVAC Model in the SOLOMON bands:
    # solomonFluxEUVAC, solomonIrrEUVAC = solomon.solomon(f107, f107a, model='EUVAC')
    # all_ints_euvac_solomon = []
    # for n in range(solomonIrrEUVAC.shape[1]):
    #     all_ints_euvac_solomon.append(simpson(solomonIrrEUVAC[:, n]))
    #
    # # Call the HEUVAC Model:
    # heuvac_wav, heuvacFlux, heuvacIrr, _, _, _ = heuvac.heuvac(f107, f107a, torr=True)
    # all_ints_heuvac = []
    # for o in range(heuvacIrr.shape[1]):
    #     all_ints_heuvac.append(simpson(heuvacIrr[:, o]))

    tableFilename = '../src/EUVpy/data/neuvac_table.txt'
    fileLoc = here.joinpath(tableFilename)
    out = neuvac.aetherFile(fileLoc)

    sys.exit(0)