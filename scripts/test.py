# Contains drafts for creating unit tests.
import numpy as np
from EUVpy.NEUVAC import neuvac
from EUVpy.tools import processIndices
from EUVpy.empiricalModels.models.EUVAC import euvac
from EUVpy.empiricalModels.models.SOLOMON import solomon
from EUVpy.empiricalModels.models.HEUVAC import heuvac
from EUVpy.tools import ensemble
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

#-----------------------------------------------------------------------------------------------------------------------
# Helper functions:
def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    # https://stackoverflow.com/questions/5595425/how-to-compare-floats-for-almost-equality-in-python
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def all_same(items):
    return all(x == items[0] for x in items)
#-----------------------------------------------------------------------------------------------------------------------

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
    #
    # # Put the NEUVAC coefficients into a format that Aether can use:
    # tableFilename = '../src/EUVpy/data/neuvac_table.txt'
    # fileLoc = here.joinpath(tableFilename)
    # out = neuvac.aetherFile(fileLoc)
    #
    # Generate NEUVAC Irradiances but put them in a form that GITM can directly use:
    # Get some F10.7 for the entirety of 2022
    # f107times, f107, f107a, f107b = processIndices.getCLSF107('2022-01-01', '2023-01-01', truncate=False)
    #
    # # Generate the NEUVAC irradiance but output them to a file that can be used by GITM
    # out = neuvac.gitmNEUVAC(f107times, f107, f107b)
    #
    # # Verify that the data in the generated file is reasonable:
    # with open(out, 'r') as outfile:
    #     contents = outfile.readlines()
    #
    # Run the irradiance ensemble command, with NEUVAC as the model:
    # f107times, f107, f107a, f107b = processIndices.getCLSF107('2014-08-01', '2014-10-01', truncate=False)
    #
    # # Run the NEUVAC ensemble - do so with the version of NEUVAC in the EUVAC bins...
    # iterations = 50  # Number of ensemble members to generate
    # ensemble_NeuvacIrr, ensemble_average_NeuvacIrr, ensemble_stddev_NeuvacIrr = ensemble.irradiance_ensemble(f107,
    #                                                                                                          f107b,
    #                                                                                                          iterations=iterations,
    #                                                                                                          model='NEUVAC-E')
    #
    # sample_spectrum = np.array([3.45237838e-04, 1.08617287e-04, 6.02517592e-04, 4.29327521e-04, 7.41131797e-05,
    #                             7.82372944e-05, 1.66828616e-04, 5.44773995e-04, 1.31203587e-05, 1.84494178e-04,
    #                             5.25233401e-05, 1.12263426e-04, 4.88559021e-05, 1.62028541e-05, 6.36173582e-05,
    #                             5.69629170e-05, 2.46465099e-05, 5.77028724e-05, 2.54878162e-05, 3.16616388e-05,
    #                             5.57739347e-05, 2.29024036e-05, 1.76708493e-05, 1.09410397e-05, 2.06660115e-05,
    #                             1.20924642e-05, 1.38322114e-05, 1.47788316e-05, 5.70492198e-05, 1.12358584e-04,
    #                             1.93896513e-04, 1.56549582e-04, 1.49034033e-04, 5.55388980e-05, 9.73276207e-05,
    #                             6.74656691e-05, 7.97451541e-05])
    # new_spectrum = ensemble_average_NeuvacIrr[30, :]
    #
    # results = [isclose(new_spectrum[int(idx)], sample_spectrum[int(idx)], rel_tol=1e1) for idx in np.linspace(0, len(new_spectrum)-1, len(new_spectrum))]

    sys.exit(0)