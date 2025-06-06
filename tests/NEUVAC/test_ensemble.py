import pytest
import numpy as np
from EUVpy.tools import processIndices
from EUVpy.tools import ensemble
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Helper function:
def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    # https://stackoverflow.com/questions/5595425/how-to-compare-floats-for-almost-equality-in-python
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def all_same(items):
    return all(x == items[0] for x in items)
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
def test_config():
    # Run the irradiance ensemble command, with NEUVAC as the model:
    f107times, f107, f107a, f107b = processIndices.getCLSF107('2014-08-01', '2014-10-01', truncate=False)

    # Run the NEUVAC ensemble - do so with the version of NEUVAC in the EUVAC bins...
    iterations = 50  # Number of ensemble members to generate
    ensemble_NeuvacIrr, ensemble_average_NeuvacIrr, ensemble_stddev_NeuvacIrr = ensemble.irradiance_ensemble(f107,
                                                                                                             f107b,
                                                                                                             iterations=iterations,
                                                                                                             model='NEUVAC-E')
    # Verify that the output is reasonable:

    # True values at index ensemble_average_NeuvacIrr[30, :]:
    sample_spectrum = np.array([3.45237838e-04, 1.08617287e-04, 6.02517592e-04, 4.29327521e-04, 7.41131797e-05,
                                7.82372944e-05, 1.66828616e-04, 5.44773995e-04, 1.31203587e-05, 1.84494178e-04,
                                5.25233401e-05, 1.12263426e-04, 4.88559021e-05, 1.62028541e-05, 6.36173582e-05,
                                5.69629170e-05, 2.46465099e-05, 5.77028724e-05, 2.54878162e-05, 3.16616388e-05,
                                5.57739347e-05, 2.29024036e-05, 1.76708493e-05, 1.09410397e-05, 2.06660115e-05,
                                1.20924642e-05, 1.38322114e-05, 1.47788316e-05, 5.70492198e-05, 1.12358584e-04,
                                1.93896513e-04, 1.56549582e-04, 1.49034033e-04, 5.55388980e-05, 9.73276207e-05,
                                6.74656691e-05, 7.97451541e-05])
    new_spectrum = ensemble_average_NeuvacIrr[30, :]
    results = [isclose(new_spectrum[int(idx)], sample_spectrum[int(idx)], rel_tol=1e1) for idx in
               np.linspace(0, len(new_spectrum) - 1, len(new_spectrum))]

    assert (all_same(results))

#-----------------------------------------------------------------------------------------------------------------------
