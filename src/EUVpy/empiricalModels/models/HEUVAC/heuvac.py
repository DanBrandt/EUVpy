# Run HEUVAC
# HEUVAC was composed by Dr. Phil Chamberlin in 2006:
# Richards, Philip G., Thomas N. Woods, and William K. Peterson. "HEUVAC: A new high resolution solar EUV proxy model."
# Advances in Space Research 37.2 (2006): 315-322.

#-----------------------------------------------------------------------------------------------------------------------
# Top-level Imports:
import numpy as np
import os
from tqdm import tqdm
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Unchangeable filenames (lines 167-169 in HEUVAC-Driver.for), where HEUVAC outputs are stored:
topDir = os.getcwd()
directory = '../empiricalModels/models/HEUVAC/'
torrFluxFile = 'Torr-37-bins.txt'
userFluxFile = 'flux-User-bins-10A.txt'
userIonizationFile = 'XS-User-bins-10A.txt'
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Local Imports:
from src.EUVpy import tools
import src.EUVpy.tools.spectralAnalysis
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Helper Functions:
def writeInputFile(F107, F107A):
    """
    Write the formatted input file for running HEUVAC.
    :param F107: float
        A single value of F10.7.
    :param F107A: float
        A single value for the 81-day averaged F10.7, centered on the current day.
    :return:
        Returns nothing.
    """
    filename = 'HEUVAC-scratch.TXT'
    with open(filename, 'w') as heuvacFile:
        heuvacFile.write(str(F107)+'\n')
        heuvacFile.write(str(F107A)+'\n')
        heuvacFile.write(str(10)) # We manually keep the 10 Angstrom bin width (the highest resolution of HEUVAC)
    print(filename)

def getTorr(fluxFile):
    """
    Read the output file from HEUVAC in the Torr bins.
    :param fluxFile: str
        The filename where the HEUVAC fluxes in the Torr bins have been output.
    :return wav: ndarray
        The wavelength bin centers for the Torr bins (Angstroms).
    :return flux: ndarray
        The HEUVAC flux in the Torr bins (W/m2/nm).
    """
    wavs = np.zeros(37)
    fluxes = np.zeros(37)
    irrs = np.zeros(37)
    with open(fluxFile) as myFile:
        fileData = myFile.readlines()
        i = 0
        j = 0
        for line in fileData:
            if i > 0:
                wavs[j] = float(line.split()[1])
                fluxes[j] = float(line.split()[-1])
                irrs[j] = src.EUVpy.tools.spectralAnalysis.spectralIrradiance(fluxes[j], wavs[j])
                j += 1
            i += 1
    # The arrays will have values from largest to smallest - they should be flipped before being returned:
    wav = np.flip(wavs)
    flux = np.flip(fluxes)
    irr = np.flip(irrs)
    return wav, flux, irr

def getFlux(userFile):
    """
    Read the output file from HEUVAC in the Torr bins.
    :param userFile: str
        The filename where the HEUVAC fluxes in the user-defined bins have been output.
    :return wav: ndarray
        The wavelength bin centers for the Torr bins (Angstroms).
    :return flux: ndarray
        The HEUVAC flux in the user-defined bins (ph/cm2/s).
    :return: irr
        The HEUVAC irradiance in the user-defined bins (W/m2/nm).
    """
    wav = np.zeros(106)
    flux = np.zeros(106)
    irr = np.zeros(106)
    with open(userFile) as myFile:
        fileData = myFile.readlines()
        i = 0
        j = 0
        for line in fileData:
            if i >= 2:
                wav[j] = float(line.split()[3])
                flux[j] = float(line.split()[-2])
                irr[j] = float(line.split()[-1])
                j += 1
            i += 1
    return wav, flux, irr

def heuvac(F107, F107A, torr=True, statsFiles=None):
    """
    Call the HEUVAC Fortran code for each F10.7, F10.7A pair.
    :param F107: arraylike
        The values of [daily] F10.7.
    :param F107A: arraylike
        The values of 81-day averaged F10.7, centered on the current day.
    :param torr: bool
        Controls whether or not the binned data returned is in the 37 standard Torr et al bins or in the high-resolution
        10 Angstrom-wide bins (the standard high resolution of HEUVAC). Default is True.
    :return heuvacWav: ndarray
        The bin center wavelengths for the HEUVAC data.
    :return heuvacFlux: ndarray
        The solar EUV flux in different wavelength bins returned from HEUVAC.
    :return heuvacIrr: ndarray
        The solar EUV irradiance in different wavelength bins returend from HEUVAC.
    """
    if torr==True:
        if type(F107) == np.ndarray:
            heuvacFlux = np.zeros((len(F107), 37))  # Columns represent each wavelength band 37 (59).
            heuvacIrr = np.zeros((len(F107), 37))
        else:
            heuvacFlux = np.zeros((1, 37))
            heuvacIrr = np.zeros((1, 37))
            F107 = np.array([F107])
            F107A = np.array([F107A])
    else:
        if type(F107) == np.ndarray:
            heuvacFlux = np.zeros((len(F107), 106))  # Columns represent each wavelength band 37 (59).
            heuvacIrr = np.zeros((len(F107), 106))
        else:
            heuvacFlux = np.zeros((1, 106))
            heuvacIrr = np.zeros((1, 106))
            F107 = np.array([F107])
            F107A = np.array([F107A])
    if not statsFiles:
        # Loop across all the F10.7 values:
        os.chdir(directory)
        for i in tqdm(range(heuvacIrr.shape[0])):
            # Write the input file and run HEUVAC:
            writeInputFile(F107[i], F107A[i])
            # Run base HEUVAC:
            os.system('./HEUVAC.exe')
            # Read in the fluxes in the Torr bins (37 bins) and the user-specified bins:
            torrWav, torrFlux, torrIrr = getTorr(torrFluxFile)
            userWav, userFlux, userIrr = getFlux(userFluxFile)
            if torr == True:
                heuvacFlux[i, :] = torrFlux
                irrRes = torrIrr * (1e4)  # Convert to W/m^2
                heuvacIrr[i, :] = irrRes
            else:
                heuvacFlux[i, :] = userFlux
                heuvacIrr[i, :] = userIrr * (1e4)  # Convert to W/m^2
        os.chdir(topDir)
        if torr == True:
            heuvacWav = torrWav
        else:
            heuvacWav = userWav
        return heuvacWav, heuvacFlux, heuvacIrr, None, None, None
    else:
        perturbedEuvIrradiance = np.zeros_like(heuvacIrr)
        savedPerts = np.zeros_like(heuvacIrr)
        # Include statistical data for calculating uncertainties via perturbations:
        corMatFile = statsFiles[0]  # '../../../experiments/corMatEUVAC.pkl'
        corMatHEUVAC = src.EUVpy.tools.toolbox.loadPickle(corMatFile)
        sigmaFileHEUVAC = statsFiles[1]  # '../../../experiments/sigma_EUVAC.pkl'
        STDHeuvacResids = src.EUVpy.tools.toolbox.loadPickle(sigmaFileHEUVAC)
        os.chdir(directory)
        # Loop across all the F10.7 values:
        for i in tqdm(range(len(F107))):
            # Write the input file and run HEUVAC:
            writeInputFile(F107[i], F107A[i])
            P_n = []
            A_j_vals = []
            # Loop across all the bands to obtain perturbations:
            for j in range(37):
                # Percentage perturbation:
                P_j = np.random.normal(0, 1.0)
                P_n.append(P_j)
                P_1 = P_n[0]
                # Normalized Correlated Perturbation:
                C_j1 = corMatHEUVAC[0, j]
                N_j = C_j1 * P_1 + (1.0 - C_j1) * P_j
                # Actual Normalized Correlated Perturbation:
                A_j = STDHeuvacResids[j] * N_j
                A_j_vals.append(A_j)
                savedPerts[i, j] = A_j

            # Run base HEUVAC:
            os.system('./HEUVAC.exe')

            # Read in the fluxes in the Torr bins (37 bins) and the user-specified bins:
            torrWav, torrFlux, torrIrr = getTorr(torrFluxFile)
            userWav, userFlux, userIrr = getFlux(userFluxFile)

            # Collect the flux and irradiance into their respective arrays:
            if torr==True:
                heuvacFlux[i, :] = torrFlux
                irrRes = torrIrr*(1e4) # Convert to W/m^2
                heuvacIrr[i, :] = irrRes
                perturbedIrr = irrRes + np.asarray(A_j_vals)
                for k in range(perturbedIrr.shape[0]):
                    if perturbedIrr[k] < 0:
                        perturbedIrr[k] = 0
                perturbedEuvIrradiance[i, :] = perturbedIrr
            else:
                heuvacFlux[i, :] = userFlux
                heuvacIrr[i, :] = userIrr*(1e4) # Convert to W/m^2
                # TODO: Add perturbation functionality for variable bin widths
        os.chdir(topDir)
        if torr==True:
            heuvacWav = torrWav
        else:
            heuvacWav = userWav

        # Generate a correlation matrix of the perturbations:
        cc2 = np.zeros((37, 37))
        for iW1 in range(37):
            for iW2 in range(37):
                cc = src.EUVpy.tools.toolbox.get_cc(savedPerts[:, iW1], savedPerts[:, iW2])
                cc2[iW1, iW2] = cc

    return heuvacWav, heuvacFlux, heuvacIrr, perturbedEuvIrradiance, savedPerts, cc2
#-----------------------------------------------------------------------------------------------------------------------


