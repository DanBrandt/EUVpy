# Code for computing the irradiance ensemble.

#-----------------------------------------------------------------------------------------------------------------------
# Top-level imports:
import numpy as np
from tqdm import tqdm
import pathlib
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Local Imports
from EUVpy.NEUVAC import neuvac
from EUVpy.empiricalModels.models.EUVAC import euvac
from EUVpy.empiricalModels.models.HEUVAC import heuvac
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Global Variables
here = pathlib.Path(__file__).parent.resolve()
neuvac_tableFile = here.parent.joinpath('data/neuvac_table.txt')
neuvac_tableFile_Stan_Bands = here.parent.joinpath('data/neuvac_table_stan_bands.txt')
neuvacStatsFiles = [here.parent.joinpath('experiments/corMat.pkl'), here.parent.joinpath('experiments/sigma_NEUVAC.pkl')]
neuvacStatsFiles_Stan_Bands = [here.parent.joinpath('experiments/corMatStanBands.pkl'), here.parent.joinpath('experiments/sigma_NEUVAC_StanBands.pkl')]
euvacStatsFiles = [here.parent.joinpath('experiments/corMatEUVAC.pkl'), here.parent.joinpath('experiments/sigma_EUVAC.pkl')]
heuvacStatsFiles = [here.parent.joinpath('experiments/corMatHEUVAC.pkl'), here.parent.joinpath('experiments/sigma_HEUVAC.pkl')]
#-----------------------------------------------------------------------------------------------------------------------

def irradiance_ensemble(F107, F107A, iterations=100, model='NEUVAC'):
    """
    Given F10.7 and F10.7A, run an ensemble of modeled EUV irradiances. Return the ensemble average, all the individual
    members, and the standard deviations of the ensemble.

    Parameters
    ----------
    F107 : float or numpy.ndarray
        Solar flux at 10.7 cm.
    F107A : float or numpy.ndarray
        81 day-averaged solar flux at 10.7, centered on the current day.
    iterations : int
        The number of ensemble members.
    model : str
        The model with which to run the ensemble. May be 'NEUVAC', 'NEUVAC-E', 'EUVAC', or 'HEUVAC'. If the model is in
        the STAN BANDS, valid arguments include 'NEUVAC-S', 'EUVAC-S', or 'HFG'. 'NEUVAC' refers to the 59-band base
        model of NEUVAC, while 'NEUVAC-E' refers to the 37-band base model of NEUVAC.

    Returns
    -------
    ensemble : numpy.ndarray
        A 3D array where each 2D element is an ensemble.
    ensemble_avg : numpy.ndarray
        A 2D array of the ensemble average of irradiances.
    ensemble_stddev : numpy.ndarray
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
