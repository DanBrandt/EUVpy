# Code for computing solar irradiance according to the EUVAC model.
# Reference: Richard, P. G., Fennelly, J. A., and Torr, D. G., EUVAC: A solar EUV flux model for aeronomic calculations,
# Journal of Geophysical Research, 99, A5, 8981-8991, 1994

#-----------------------------------------------------------------------------------------------------------------------
# Top-level imports:
import numpy as np
from tqdm import tqdm
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Local Imports:
from EUVpy.tools import spectralAnalysis
from EUVpy import tools

#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Constants
h = 6.62607015e-34 # Planck's constant in SI units of J s
c = 299792458 # Speed of light in m s^-1
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Global variables:
euvacTable = np.array([
    [1, 50, 100, 1.200, 1.0017e-2],
    [2, 100, 150, 0.450, 7.1250e-3],
    [3, 150, 200, 4.800, 1.3375e-2],
    [4, 200, 250, 3.100, 1.9450e-2],
    [5, 256.32, 256.32, 0.460, 2.7750e-3],
    [6, 284.15, 284.15, 0.210, 1.3768e-1],
    [7, 250, 300, 1.679, 2.6467e-2],
    [8, 303.31, 303.31, 0.800, 2.5000e-2],
    [9, 303.78, 303.78, 6.900, 3.3333e-3],
    [10, 300, 350, 0.965, 2.2450e-2],
    [11, 368.07, 368.07, 0.650, 6.5917e-3],
    [12, 350, 400, 0.314, 3.6542e-2],
    [13, 400, 450, 0.383, 7.4083e-3],
    [14, 465.22, 465.22, 0.290, 7.4917e-3],
    [15, 450, 500, 0.285, 2.0225e-2],
    [16, 500, 550, 0.452, 8.7583e-3],
    [17, 554.37, 554.37, 0.720, 3.2667e-3],
    [18, 584.33, 584.33, 1.270, 5.1583e-3],
    [19, 550, 600, 0.357, 3.6583e-3],
    [20, 609.76, 609.76, 0.530, 1.6175e-2],
    [21, 629.73, 629.73, 1.590, 3.3250e-3],
    [22, 600, 650, 0.342, 1.1800e-2],
    [23, 650, 700, 0.230, 4.2667e-3],
    [24, 703.36, 703.36, 0.360, 3.0417e-3],
    [25, 700, 750, 0.141, 4.7500e-3],
    [26, 765.15, 765.15, 0.170, 3.8500e-3],
    [27, 770.41, 770.41, 0.260, 1.2808e-2],
    [28, 789.36, 789.36, 0.702, 3.2750e-3],
    [29, 750, 800, 0.758, 4.7667e-3],
    [30, 800, 850, 1.625, 4.8167e-3],
    [31, 850, 900, 3.537, 5.6750e-3],
    [32, 900, 950, 3.000, 4.9833e-3],
    [33, 977.02, 977.02, 4.400, 3.9417e-3],
    [34, 950, 1000, 1.475, 4.4167e-3],
    [35, 1025.72, 1025.72, 3.500, 5.1833e-3],
    [36, 1031.91, 1031.91, 2.100, 5.2833e-3],
    [37, 1000, 1050, 2.467, 4.3750e-3]
    ])
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Functions:
def refSpec(i):
    """
    Helper function for the EUVAC Model: returns the standard solar flux in 37 bands from the F74113 Spectrum
    (pp. 584-585 in Schunk and Nagy).

    For additional resources, see the following:
    - Richard, P. G., Fennelly, J. A., and Torr, D. G., EUVAC: A solar EUV flux model for aeronomic calculations,
    Journal of Geophysical Research, 99, A5, 8981-8992, 1994.
    - Heroux, L. and Hinteregger, H. E., Aeronomical Reference Spectrum for Solar UV Below 2000 A, Journal of
    Geophysical Research, 83, A11, 1978.

    Parameters
    ----------
    i : int
        The index for the wavelength. Must be between 0 and 37.

    Returns
    -------
    F74113_i : float
        The reference solar flux in units of photons/m^2/s^.
    A_i : float
        The scaling factor for the wavelength interval.
    """
    lookUpIdx = np.where(euvacTable[:, 0] == i)[0]
    F74113_i = euvacTable[lookUpIdx, 3][0]*(1e13) # Multiply by 1e13 to obtain photons m^-2 s^-1.
    A_i = euvacTable[lookUpIdx, 4][0]
    return F74113_i, A_i

def euvac(F107, F107A, statsFiles=None):
    """
    Compute the solar flux from F10.7, according to the EUVAC model. Return the solar flux across 37 wavelength
    bands in units of photons m^-2 s^-1.

    Parameters
    ----------
    F107 : numpy.ndarray
        Values of the F10.7 solar flux.
    F107A : numpy.ndarray
        Values of the 81-day averaged solar flux, centered on the present day.
    statsFiles : list
        A 2 element list where the first element is a file containing the 37x37 correlation matrix and the second
        element is a file containing the 1x37 standard deviation values for EUVAC. Optional.

    Returns
    -------
    euvacFlux: numpy.ndarray
        Values of the solar radiant flux in 37 distinct wavelength bands. In photons/m^2/s.
    euvacIrr: numpy.ndarray
        Values of the solar spectral irradiance in 37 distinct wavelength bands. In W/m^2.
    """
    P = (F107A + F107)/2.0
    if type(F107) == np.ndarray:
        euvacFlux = np.zeros((len(F107), 37)) # Columns represent each wavelength band 37 (59).
        euvacIrr = np.zeros((len(F107), 37))
    else:
        euvacFlux = np.zeros((1, 37))
        euvacIrr = np.zeros((1, 37))
        P = np.array([P])
    if not statsFiles:
        # Loop across the F10.7 values
        for i in tqdm(range(euvacIrr.shape[0])):
            k = 0
            # Loop across all the bands:
            for j in range(37):
                # Compute the flux:
                wav = 0.5 * (euvacTable[j, 2] + euvacTable[j, 1])
                F74113_i, A_i = refSpec(j + 1)
                fluxFactor = (1. + A_i * (P[i] - 80.))
                photonFlux = (F74113_i) * fluxFactor
                try:
                    photonFlux[photonFlux < 0] = 0
                except:
                    if photonFlux < 0:
                        photonFlux = 0
                euvacFlux[i, k] = photonFlux
                irrRes = spectralAnalysis.spectralIrradiance(photonFlux, wavelength=wav)
                euvacIrr[i, k] = irrRes
                k += 1
        return euvacFlux, euvacIrr, None, None, None
    else:
        perturbedEuvIrradiance = np.zeros_like(euvacIrr)
        savedPerts = np.zeros_like(euvacIrr)
        # Include statistical data for calculating uncertainties via perturbations:
        corMatFile = statsFiles[0]  # '../../../experiments/corMatEUVAC.pkl'
        corMatEUVAC = src.EUVpy.tools.toolbox.loadPickle(corMatFile)
        sigmaFileEUVAC = statsFiles[1]  # '../../../experiments/sigma_EUVAC.pkl'
        STDEuvacResids = src.EUVpy.tools.toolbox.loadPickle(sigmaFileEUVAC)
        # Loop over the F10.7 values
        for i in tqdm(range(euvacIrr.shape[0])):
            k = 0
            P_n = []
            # Loop across all the bands:
            for j in range(37):
                # Percentage perturbation:
                P_j = np.random.normal(0, 1.0)
                P_n.append(P_j)
                P_1 = P_n[0]
                # Normalized Correlated Perturbation:
                C_j1 = corMatEUVAC[0, j]  # Only consider correlation with the first wavelength bin
                N_j = C_j1 * P_1 + (1.0 - C_j1) * P_j
                # Actual Normalized Correlated Perturbation:
                A_j = STDEuvacResids[j] * N_j
                # Compute the flux:
                wav = 0.5*(euvacTable[j, 2] + euvacTable[j, 1])
                # dWav = euvacTable[j, 2] - euvacTable[j, 1]
                # if dWav == 0:
                #     dWav = None
                F74113_i, A_i = refSpec(j+1)
                fluxFactor = (1. + A_i*(P[i]-80.))
                # if fluxFactor < 0.8:
                    # fluxFactor = 0.8
                photonFlux = (F74113_i)*fluxFactor
                # If P-80 is negative, set the flux to ZERO.
                try:
                    photonFlux[photonFlux < 0] = 0
                except:
                    if photonFlux < 0:
                        photonFlux = 0
                euvacFlux[i, k] = photonFlux
                irrRes = spectralAnalysis.spectralIrradiance(photonFlux, wavelength=wav)
                euvacIrr[i, k] = irrRes
                perturbedEuvIrradiance[i, k] = irrRes + A_j
                savedPerts[i, j] = A_j
                k += 1

        # Generate a correlation matrix of the perturbations:
        cc2 = np.zeros((37, 37))
        for iW1 in range(37):
            for iW2 in range(37):
                cc = src.EUVpy.tools.toolbox.get_cc(savedPerts[:, iW1], savedPerts[:, iW2])
                cc2[iW1, iW2] = cc

        return euvacFlux, euvacIrr, perturbedEuvIrradiance, savedPerts, cc2
#-----------------------------------------------------------------------------------------------------------------------

