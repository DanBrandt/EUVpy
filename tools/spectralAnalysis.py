# Code for computing converting flux into irradiance and the reverse, along with miscelleneous helper functions.

#-----------------------------------------------------------------------------------------------------------------------
# Top-level imports:
import numpy as np
from tqdm import tqdm
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Global Variables
neuvac_tableFile = '../NEUVAC/neuvac_table.txt'
neuvac_tableFile_Stan_Bands = '../NEUVAC/neuvac_table_stan_bands.txt'
neuvacStatsFiles = ['../experiments/corMat.pkl', '../experiments/sigma_NEUVAC.pkl']
neuvacStatsFiles_Stan_Bands = ['../experiments/corMatStanBands.pkl', '../experiments/sigma_NEUVAC_StanBands.pkl']
euvacStatsFiles = ['../experiments/corMatEUVAC.pkl', '../experiments/sigma_EUVAC.pkl']
heuvacStatsFiles = ['../experiments/corMatHEUVAC.pkl', '../experiments/sigma_HEUVAC.pkl']
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Local Imports
from NEUVAC import neuvac
from empiricalModels.models.EUVAC import euvac
from empiricalModels.models.HEUVAC import heuvac
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Constants
h = 6.62607015e-34 # Planck's constant in SI units of J s
c = 299792458 # Speed of light in m s^-1
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Functions:
def spectralIrradiance(photonFlux, wavelength):
    """
    Convert the photon flux to the corresponding spectral irradiance, given a specific wavelength.
    :param: photonFlux: ndarray, float, or int
        Photon flux in units of photons s^-1 m^-2. For a singular wavelength, units are photons m^-2
    :param: wavelength: float
        A specific wavelength in Angstroms.
    :return: irradiance: ndarray or float
        The corresponding spectral irradiance in units of W/m^2/nm.
    """
    photonEnergy = (h*c) / (wavelength*1e-10) # Convert the wavelength in the denominator to meters.
    irradiance = photonFlux * photonEnergy
    # if dWavelength != None:
    #     irradiance = photonFlux * photonEnergy * (1./(dWavelength*0.1)) # Multiply the denominator by 0.1 in order to convert from an Angstrom interval to a nanometer interval.
    # else:
    #     irradiance = photonFlux * photonEnergy / wavelength
    return irradiance

def spectralFlux(irradiance, wavelength): #, dWavelength=10):
    """
    Convert the spectral irradiance to spectral flux, given a specific wavelength.
    :param: irradiance: ndarray, float, or int
        [Solar] spectral irradiance for a specific wavelength, in W/m^2/nm.
    :param: wavelength: float
        A specific wavelength in Angstroms.
    :param: dWavelength: float or int
        Wavelength bin width in Angstroms. Default is 1.
    :return: photonFlux: ndarray or float
        The corresponding spectral flux in units of Watts.
    """
    photonEnergy = (h * c) / (wavelength * 1e-10)
    photonFlux = irradiance / photonEnergy # (irradiance * dWavelength * 0.1) / photonEnergy
    return photonFlux

def irradiance_ensemble(F107, F107A, iterations=100, model='NEUVAC'):
    """
    Given F10.7 and F10.7A, run an ensemble of modeled EUV irradiances. Return the ensemble average, all of the members,
    and the standard deviations of the ensemble.
    :param F107: float or ndarray
        Solar flux at 10.7 cm.
    :param F107A: float or ndarray
        81 day-averaged solar flux at 10.7, centered on the current day.
    :param iterations: int
        The number of ensemble members.
    :param model: str
        The model with which to run the ensemble. May be 'NEUVAC', 'NEUVAC-E', 'EUVAC', or 'HEUVAC'. If the model is in
        the STAN BANDS, valid arguments include 'NEUVAC-S', 'EUVAC-S', or 'HFG'. 'NEUVAC' refers to the 59-band base
        model of NEUVAC, while 'NEUVAC-E' refers to the 37-band base model of NEUVAC.
    :return ensemble: ndarray
        A 3D array where each 2D element is an ensemble.
    :return ensemble_avg: ndarray
        A 2D array of the ensemble average of irradiances.
    :return ensemble_stddev: ndarray
        A 2D array of the ensemble standard deviations.
    """
    # Get the shape of the third dimension
    if model[-1] == 'S':
        lastDim = 22
    elif model == 'NEUVAC':
        lastDim = 59
    else:
        lastDim = 37
    # Instantiate the ensemble:
    ensemble = np.zeros((iterations, len(F107), lastDim))
    # Fill the ensemble:
    for i in tqdm(range(iterations)):
        if model=='NEUVAC-E':
            _, perturbedEuvIrradiance, _, _ = neuvac.neuvacEUV(F107, F107A, bands='EUVAC', tableFile=neuvac_tableFile, statsFiles=neuvacStatsFiles)
        elif model=='NEUVAC':
            _, perturbedEuvIrradiance, _, _ = neuvac.neuvacEUV(F107, F107A, bands=None, tableFile=neuvac_tableFile, statsFiles=neuvacStatsFiles)
        elif model=='EUVAC':
            _, _, perturbedEuvIrradiance, _, _ = euvac.euvac(F107, F107A, statsFiles=euvacStatsFiles)
        elif model=='HEUVAC':
            _, _, _, perturbedEuvIrradiance, _, _ = heuvac.heuvac(F107, F107A, statsFiles=heuvacStatsFiles)
        elif model=='NEUVAC-S' or model=='EUVAC-S' or model=='HFG':
            if model == 'NEUVAC-S':
                _, perturbedEuvIrradiance, _, _ = neuvac.neuvacEUV(F107, F107A, bands='SOLOMON',
                                                                   tableFile=neuvac_tableFile_Stan_Bands,
                                                                   statsFiles=neuvacStatsFiles)
            elif model == 'EUVAC-S':
                pass
            elif model == 'HFG':
                pass
            else:
                raise ValueError
        else:
            raise ValueError('The chosen model must either be: NEUVAC-E, NEUVAC, EUVAC, HEUVAC, or the SOLOMON version'
                             'of the aforementioned (NEUVAC-S, EUVAC-S, or HFG)!')
        ensemble[i, :, :] = perturbedEuvIrradiance
    # Compute the ensemble average and ensemble standard deviations:
    ensemble_average = np.nanmean(ensemble, axis=0)
    ensemble_stddev = np.nanstd(ensemble, axis=0)
    return ensemble, ensemble_average, ensemble_stddev
#-----------------------------------------------------------------------------------------------------------------------
