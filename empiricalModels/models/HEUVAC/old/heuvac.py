# Code for computing solar irradiance according to the HEUVAC model.
# Reference: Richard, P. G., Woods, T. N., and Peterson, W. K., HEUVAC: A new high resolution solar EUV proxy model,
# Advances in Space Research, 37, 2, 315-322, 2006.

#-----------------------------------------------------------------------------------------------------------------------
# Top-level imports:
import numpy as np
from random import randrange
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Local Imports:
from tools.spectralAnalysis import spectralIrradiance
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Global Variable:
heuvacTable = np.array([
        [1, 50, 100, 3.106, 2.58],
        [2, 100, 150, 2.1, 1.53],
        [3, 150, 200, 3.5, 1.5],
        [4, 200, 250, 1.888, 2.48],
        [5, 256.32, 256.32, 4.4, 1.38],
        [6, 284.15, 284.15, 4.433, 1.85],
        [7, 250, 300, 4.734, 1.93],
        [8, 303.31, 303.31, 2.42, 1.72],
        [9, 303.78, 303.78, 0.731, 1.29],
        [10, 300, 350, 0.7, 1.35],
        [11, 368.07, 368.07, 0.26, 2.24],
        [12, 350, 400, 0.1832, 1.27],
        [13, 400, 450, 0.353, 1.78],
        [14, 465.22, 465.22, 0.36, 1.3],
        [15, 450, 500, 0.5454, 1.65],
        [16, 500, 550, 0.712, 2.15],
        [17, 554.37, 554.37, 0.795, 1.51],
        [18, 584.33, 584.33, 0.265, 3.02],
        [19, 550, 600, 0.349, 1.03],
        [20, 609.76, 609.76, 0.791, 2.15],
        [21, 629.73, 629.73, 0.5, 1.4],
        [22, 600, 650, 0.775, 2.19],
        [23, 650, 700, 0.817, 1.9],
        [24, 703.36, 703.36, 0.29, 1.73],
        [25, 700, 750, 0.649, 1.96],
        [26, 765.15, 765.15, 1.051, 2.94],
        [27, 770.41, 770.41, 0.65, 1.64],
        [28, 789.36, 789.36, 0.4824, 8.02],
        [29, 750, 800, 6.3736, 1.43],
        [30, 800, 850, 0.739, 3.71],
        [31, 850, 900, 4.51, 1.84],
        [32, 900, 950, 0.21, 14.39],
        [33, 977.02, 977.02, 0.46, 1.57],
        [34, 950, 1000, 3.558, 1.78],
        [35, 1025.72, 1025.72, 6.26, 1.51],
        [36, 1031.91, 1031.91, 0.28, 1.5],
        [37, 1000, 1050, 0.4721, 1.82]
    ])
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Functions:
def refSpecHEUVAC(i):
    """
    Return the standard solar flux in 37 bands using HEUVAC.
    :param: i: int
        The index for the wavelength. Must be between 0 and 37.
    :return: flux68: float
        The HEUVAC reference solar flux in units of photons m^-2 s^-1.
    :return: SEEFAC: float
        The scaling factor for the wavelength interval.
    """
    lookUpIdx = np.where(heuvacTable[:, 0] == i)[0]
    flux68 = heuvacTable[:, 3][lookUpIdx][0]
    SEEFAC = np.flip(heuvacTable[:, 4])[lookUpIdx][0]
    return flux68, SEEFAC

def heuvac(F107, F107A):
    """
    Compute the solar flux from F10.7, according to the HEUVAC model. Return the solar flux across 37 wavelength
    bands in units of photons m^-2 s^-1.
    Originally written by P. Richards in January 2009 in Fortran.
    :param F107: ndarray
        Values of the F10.7 solar flux.
    :param F107A: ndarray
        Values of the 81-day averaged solar flux, centered on the present day.
    :return heuvacFlux: ndarray
        Values of the solar radiant flux in 37 distinct wavelength bands.
    :return heuvacIrr: ndarray
        Values of the solar EUV irradiance in 37 distinct wavelength bands.
    """
    P = (F107A + F107)/2.0
    if type(F107) == np.ndarray:
        heuvacFlux = np.zeros((len(F107), 37)) # Columns represent each wavelength band 37
        heuvacIrr = np.zeros((len(F107), 37))
    else:
        heuvacFlux = np.zeros((1, 37))
        heuvacIrr = np.zeros((1, 37))
    # Loop over each wavelength band and compute the flux and irradiance:
    for i in range(37):
        wav = 0.5 * (heuvacTable[i, 2] + heuvacTable[i, 1])
        dWav = heuvacTable[i, 2] - heuvacTable[i, 1]
        if dWav == 0:
            dWav = None
        flux68, SEEFAC = refSpecHEUVAC(i+1)
        A = (SEEFAC-1.) / 128.0
        B = 1 - A*68.0
        UVFAC = A*P + B
        try:
            UVFAC[UVFAC < 0.8] = 0.8
        except:
            if UVFAC < 0.8:
                UVFAC = 0.8
        photonFlux = flux68*1.0e9*UVFAC
        # If photonFlux is negative, set the flux to ZERO.
        try:
            photonFlux[photonFlux < 0] = 0
        except:
            if photonFlux < 0:
                photonFlux = 0
        photonFlux = photonFlux / (1e-4) # Divide by a factor of 1e-4 to convert from cm^-2 s^-1 to m^-2 s^-1
        heuvacFlux[:, i] = photonFlux
        heuvacIrr[:, i] = spectralIrradiance(photonFlux, wavelength=wav, dWavelength=dWav)
    # Flip the results before they are returned:
    # heuvacFlux = np.flip(heuvacFlux, 1)
    # heuvacIrr = np.flip(heuvacIrr, 1)
    return heuvacFlux, heuvacIrr

# TODO: UPDATE THIS TO USE PYTHON-WRAPPED VERSION OF THE FORTRAN CODE WRITTEN BY PHIL RICHARDS


#-----------------------------------------------------------------------------------------------------------------------
