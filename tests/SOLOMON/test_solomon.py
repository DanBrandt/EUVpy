import pytest
from EUVpy.empiricalModels.models.SOLOMON import solomon
from EUVpy.tools import processIndices
from scipy.integrate import simpson
import numpy as np
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

    # Generate some irradiances for the HFG model:
    solomonFluxHFG, solomonIrrHFG = solomon.solomon(f107, f107a, model='HFG')
    integrals = [6.55704969816892e-07, 1.680484239373774e-05, 0.0012430751892519498, 0.002330916819952905,
                 0.023751001626232005, 0.02587941084129497, 0.0803080693615398, 0.04979410405554018,
                 0.16099968222979127, 0.04704180853120706, 0.061030137904925105, 0.015582566288276693,
                 0.012575201807334744, 0.025315524584710522, 0.0391976374718764, 0.008575833635471517,
                 0.004148143049799358, 0.01130302041502086, 0.006585123131982026, 0.04444690518941401,
                 0.0397885641751608, 0.04370576417351968]
    all_ints = np.zeros_like(integrals)
    for i in range(solomonIrrHFG.shape[1]):
        all_ints[i] = simpson(solomonIrrHFG[:, i])
        assert (isclose(all_ints[i], integrals[i], rel_tol=1e-10))

    # Generate some irradiances for the EUVAC Model in the SOLOMON bands:
    solomonFluxEUVAC, solomonIrrEUVAC = solomon.solomon(f107, f107a, model='EUVAC')
    integrals_e = [8.438500737538878e-07, 3.2749536340391615e-05, 0.0011343034700311068, 0.0030751116787642355,
                   0.06517048080504681, 0.07391306749833923, 0.18288559534521728, 0.09039490962196206,
                   0.18410722304448676, 0.04096834056458914, 0.05509968188041358, 0.013698503025717509,
                   0.010954565359970336, 0.018865470573571762, 0.02918318525373133, 0.006732711490980054,
                   0.0027727818850081447, 0.007462593570200892, 0.005223969328462236, 0.031640904662614704,
                   0.030051136515095066, 0.027984039519858353]
    all_ints_e = np.zeros_like(integrals_e)
    for j in range(solomonIrrEUVAC.shape[1]):
        all_ints_e[j] = simpson(solomonIrrEUVAC[:, j])
        assert (isclose(all_ints_e[j], integrals_e[j], rel_tol=1e-10))