import pytest
import numpy as np
from EUVpy.NEUVAC import neuvac
from EUVpy.tools import processIndices
from scipy.integrate import simpson
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Helper function:
def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    # https://stackoverflow.com/questions/5595425/how-to-compare-floats-for-almost-equality-in-python
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
def test_config():
    # Get some F10.7 for the entirety of a given month in 2018:
    f107times, f107, f107a, f107b = processIndices.getCLSF107('2018-01-01', '2018-12-31',
                                                              truncate=False)

    # Generate irradiances in the EUVAC bands (validate them):
    neuvacIrr, _, _, _ = neuvac.neuvacEUV(f107, f107b, bands='EUVAC')

    integrals = [0.05246648908572817, 0.023497427718739458, 0.13806354412432875, 0.07631886793016032,
                 0.016723701365305498, 0.0051957608462539235, 0.028292161459003737, 0.1577763106130986,
                 0.0009930146094422707, 0.037566784591576495, 0.01443216526079433, 0.02049586310530997,
                 0.01269908762484041, 0.005469908401913766, 0.0139109311371798, 0.013797460183869877,
                 0.008512049943612458, 0.01654861207416694, 0.006648913079692914, 0.007203248069337905,
                 0.018754505395893087, 0.005027181893561718, 0.005594639014209394, 0.0038743017631303685,
                 0.006275032969400698, 0.0041669621899827845, 0.004683215512056784, 0.004835992522933198,
                 0.01782572798660837, 0.03361503371757845, 0.05419734192300687, 0.04450112093103549,
                 0.04615363618215064, 0.017159618520822418, 0.026531640122067855, 0.019221480713113928,
                 0.022672198761951815]

    all_ints = np.zeros_like(integrals)
    for i in range(neuvacIrr.shape[1]):
        all_ints[i] = simpson(neuvacIrr[:, i])
        assert(isclose(all_ints[i], integrals[i], rel_tol=1e-15))

    # Do the same as the above with the other wavelength bands:

    assert True

