# This script performs uncertainty quantification on the NEUVAC model.

#-----------------------------------------------------------------------------------------------------------------------
# Top-level Imports:
import numpy as np
import matplotlib, sys
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from scipy.optimize import curve_fit
from scipy.stats import pearsonr
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Local Imports:
from NEUVAC import neuvac
from tools.EUV.fism2_process import read_euv_csv_file
from tools.processIrradiances import obtainFism2
from tools import toolbox
from empiricalModels.models.EUVAC import euvac

#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Directory Management
euv_folder = '../tools/EUV/'
neuvac_tableFile = '../NEUVAC/neuvac_table.txt'
neuvac_tableFile_Solomon = '../NEUVAC/neuvac_table_stan_bands.txt'
figures_folder = 'Uncertainty'
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Global Plotting Settings:
import matplotlib.pylab as pylab
params = {'legend.fontsize': 'large',
         'axes.labelsize': 'large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Execution
if __name__=="__main__":
    # Load in F10.7 data (OMNIWeb):
    omniTimesData = '../solarIndices/F107/OMNIWeb/OMNIF107times.pkl'
    omniF107Data = '../solarIndices/F107/OMNIWeb/OMNIF107vals.pkl'
    omniF107AveData = '../solarIndices/F107/OMNIWeb/OMNIF107averageVals.pkl'
    omniTimes = toolbox.loadPickle(omniTimesData)
    omniF107 = toolbox.loadPickle(omniF107Data)
    omniF107B = toolbox.loadPickle(omniF107AveData)
    # F10.7 data extends between 1963-11-28; 12:00 to 2023-09-27; 12:00.
    times = omniTimes
    F107 = omniF107
    F107A = toolbox.rollingAverage(omniF107, 81, True, True)
    F107B = omniF107B
    # Compute F10.7 uncertainties.
    rollingStdF107 = toolbox.rollingStd(F107, 2)
    rollingStdF107B = toolbox.rollingStd(F107B, 54)
    # View the results:
    plt.figure()
    plt.fill_between(np.linspace(0, len(F107) - 1, len(F107)), np.subtract(F107, rollingStdF107),
                     np.add(F107, rollingStdF107), color='cyan', linestyle='-')
    plt.plot(np.linspace(0, len(F107) - 1, len(F107)), F107, 'b-')
    plt.figure()
    plt.fill_between(np.linspace(0, len(F107B) - 1, len(F107B)), np.subtract(F107B, rollingStdF107B),
                     np.add(F107B, rollingStdF107B), color='limegreen', linestyle='-')
    plt.plot(np.linspace(0, len(F107B) - 1, len(F107B)), F107B, 'g-')

    # Generate NEUVAC data:
    neuvacIrr, _, _, _ = neuvac.neuvacEUV(F107, F107B, bands=None, tableFile=neuvac_tableFile) #, statsFiles=['corMat.pkl', 'sigma_NEUVAC.pkl']) # perturbedNeuvacIrr, savedPerts, cc2

    # Load in FISM2 data:
    euv_data_59 = read_euv_csv_file(euv_folder + 'euv_59.csv', band=False)
    mids = 0.5 * (euv_data_59['long'] + euv_data_59['short'])
    # FISM2 Results:
    fism2file = '../empiricalModels/irradiances/FISM2/daily_data_1947-2023.nc'
    myIrrTimesFISM2, wavelengthsFISM2, myIrrDataAllFISM2, myIrrUncAllFISM2 = obtainFism2(fism2file)
    # Rebin the data:
    myIrrDataWavelengthsFISM2, rebinnedIrrDataFISM2 = toolbox.newbins(wavelengthsFISM2, myIrrDataAllFISM2, euv_data_59,
                                                                    zero=False) # toolbox.rebin
    rebinnedIrrUncFISM2 = np.zeros_like(rebinnedIrrDataFISM2)
    for column in range(rebinnedIrrDataFISM2.shape[1]):
        rebinnedIrrUncFISM2[:, column] = toolbox.rollingStd(rebinnedIrrDataFISM2[:, column], 2)

    # Replace bad values with NaNs:
    myIrrDataAllFISM2Fixed = rebinnedIrrDataFISM2.copy()
    myIrrDataAllFISM2Fixed[myIrrDataAllFISM2Fixed <= 0] = np.nan # For plotting

    # Make a training set from data up to Solar Cycle 25:
    times = np.array([element + +timedelta(hours=12) for element in times])
    trainIndsOMNI = np.where(times < datetime(2019, 12, 31))[0]
    trainIndsFISM2 = np.where((myIrrTimesFISM2 >= times[0]) & (myIrrTimesFISM2 < datetime(2019, 12, 31)))[0]
    trainTimesOMNI = times[trainIndsOMNI]
    trainTimesFISM2 = myIrrTimesFISM2[trainIndsFISM2]
    trainF107 = F107[trainIndsOMNI]
    trainF107B = F107B[trainIndsOMNI]
    trainNEUVAC = neuvacIrr[trainIndsOMNI, :]
    trainFISM2 = myIrrDataAllFISM2Fixed[trainIndsFISM2, :]
    trainUncFISM2 = rebinnedIrrUncFISM2[trainIndsFISM2, :]

    # Make a test set out of Solar Cycle 25:
    testIndsOMNI = np.where(times >= datetime(2019, 12, 31))[0]
    testIndsFISM2 = np.where((myIrrTimesFISM2 > datetime(2019, 12, 31)) & (myIrrTimesFISM2 <= times[testIndsOMNI][-1]))[0]
    testTimesOMNI = times[testIndsOMNI]
    testTimesFISM2 = myIrrTimesFISM2[testIndsFISM2]
    testF107 = F107[testIndsOMNI]
    testF107B = F107B[testIndsOMNI]
    testNEUVAC = neuvacIrr[testIndsOMNI, :]
    testFISM2 = myIrrDataAllFISM2Fixed[testIndsFISM2, :]
    testUncFISM2 = rebinnedIrrUncFISM2[testIndsFISM2, :]
    # ------------------------------------------------------------------------------------------------------------------
    # UNCERTAINTY ANALYSIS

    # Harmonize the times for NEUVAC and FISM2:
    correspondingIndsFISM2 = np.where((myIrrTimesFISM2 >= times[0]) & (myIrrTimesFISM2 <= times[-1]))[0]
    correspondingIrrTimesFISM2 = myIrrTimesFISM2[correspondingIndsFISM2]
    correspondingIrrFISM2 = rebinnedIrrDataFISM2[correspondingIndsFISM2, :]

    # ------------------------------------------------------------------------------------------------------------------
    # 1: Compute the normalized cross-correlation matrix between residuals in different bins.
    residualsArray = np.subtract(neuvacIrr, correspondingIrrFISM2)
    toolbox.savePickle(residualsArray, 'residualsArray.pkl')

    corMat = toolbox.mycorrelate2d(residualsArray, normalized=True)
    toolbox.savePickle(corMat, 'corMat.pkl')
    # ------------------------------------------------------------------------------------------------------------------
    # 2: Compute the normalized standard deviation of NEUVAC irradiance residuals (in each band):
    STDNeuvacResids = np.zeros(neuvacIrr.shape[1])
    for i in range(STDNeuvacResids.shape[0]):
        STDNeuvacResids[i] = np.nanstd(residualsArray[:, i])
    # Save these values to be used later for running ensembles:
    toolbox.savePickle(STDNeuvacResids, 'sigma_NEUVAC.pkl')

    # ------------------------------------------------------------------------------------------------------------------
    # 3: View the correlation matrix for the residuals of the perturbed NEUVAC irradiances alongside the base NEUVAC irradiances:
    # fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(11, 6))
    # pos=axs[0].imshow(corMat, aspect='auto', cmap='bwr', vmin=-1.0, vmax=1.0, interpolation='none')
    # axs[0].set_xlabel('Wavelength Band')
    # axs[0].set_ylabel('Wavelength Band')
    # axs[0].set_title('Original Correlation Matrix (NEUVAC - FISM2)')
    # pos2=axs[1].imshow(cc2, aspect='auto', cmap='bwr', vmin=-1.0, vmax=1.0, interpolation='none')
    # axs[1].set_xlabel('Wavelength Band')
    # axs[1].set_ylabel('Wavelength Band')
    # axs[1].set_title('Perturbation Correlation Matrix (NEUVAC_P - NEUVAC)')
    # fig.colorbar(pos, ax=axs[0])
    # fig.colorbar(pos2, ax=axs[1])
    fig = plt.figure()
    pos = plt.imshow(corMat, aspect='auto', cmap='bwr', vmin=-1.0, vmax=1.0, interpolation='none') # 0.995
    # pos = plt.imshow(originalCorMat, aspect='auto', cmap='bwr', vmin=-1.0, vmax=1.0, interpolation='none')  # 0.995
    plt.xlabel('Wavelength Band')
    plt.ylabel('Wavelength Band')
    plt.title('Covariance Matrix (NEUVAC-59 - FISM2)')
    fig.colorbar(pos)
    plt.savefig('Uncertainty/corMat.png', dpi=300)

    # ------------------------------------------------------------------------------------------------------------------
    # 4: Look at the Correlation between FISM2 and NEUVAC in each band, and that of (NEUVAC-FISM2)^2 and NEUVAC in each band:
    pearsonVals = []
    for i in range(neuvacIrr.shape[1]):
        fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(11,6))
        # FISM2 Irradiance vs NEUVAC Irradiance
        sortInds = np.argsort(neuvacIrr[:, i])
        sortedNEUVAC = neuvacIrr[:, i][sortInds]
        sortedFISM2 = correspondingIrrFISM2[:, i][sortInds]
        popt, pcov = curve_fit(toolbox.linear, sortedNEUVAC, sortedFISM2)
        pearsonR = pearsonr(sortedNEUVAC, sortedFISM2)
        axs[0].scatter(sortedNEUVAC, sortedFISM2, color='b')
        xlims = axs[0].get_xlim()
        sample = np.linspace(xlims[0], xlims[-1], 250)
        axs[0].plot(sample, toolbox.linear(sample, *popt), 'r-', label='R='+str(np.round(pearsonR[0], 2)))
        axs[0].set_xlim(xlims)
        axs[0].set_xlabel('NEUVAC (W/m$^2$)')
        axs[0].set_ylabel('FISM2 (W/m$^2$)')
        axs[0].legend(loc='best')
        axs[0].set_xticklabels(axs[0].get_xticklabels(), rotation=45, ha='right')
        # Look at (Irradiance Predicted - Irradiance FISM2)^2 vs. Irradiance Predicted
        squareDiffs = toolbox.squareDiff(sortedNEUVAC, sortedFISM2)
        popt2, pcov2 = curve_fit(toolbox.linear, sortedNEUVAC, squareDiffs)
        pearsonR2 = pearsonr(sortedNEUVAC, squareDiffs)
        axs[1].scatter(sortedNEUVAC, squareDiffs, color='b')
        xlims2 = axs[0].get_xlim()
        sample2 = np.linspace(xlims2[0], xlims2[-1], 250)
        axs[1].plot(sample2, toolbox.linear(sample2, *popt2), 'r-', label='R='+str(np.round(pearsonR2[0], 2)))
        axs[1].set_xlim(xlims2)
        axs[1].set_yscale('log')
        axs[1].set_xlabel('NEUVAC (W/m$^2$)')
        axs[1].set_ylabel('(NEUVAC - FISM2)$^2$ (W/m$^2$)')
        axs[1].legend(loc='best')
        axs[1].set_xticklabels(axs[0].get_xticklabels(), rotation=45, ha='right')
        fig.suptitle('FISM2 and NEUVAC Correlation and Squared Differences ('+str(mids[i])+' $\mathrm{\AA}$)')
        plt.tight_layout()
        plt.savefig('Uncertainty/correlation_sqdf_band_'+str(i+1)+'.png', dpi=300)
        pearsonVals.append([pearsonR[0], pearsonR2[0]])
    # ------------------------------------------------------------------------------------------------------------------
    # 5: Do all of the above, but for the SOLOMON version of NEUVAC:

    # Generate NEUVAC data:
    neuvacIrrSolomon, _, _, _ = neuvac.neuvacEUV(F107, F107A, bands='SOLOMON', tableFile=neuvac_tableFile_Solomon) #, statsFiles=['corMatStanBands.pkl', 'sigma_NEUVAC_StanBands.pkl'])

    # Load in FISM2 STAN BAND data:
    # FISM2 Stan Band Results:
    fism2file = '../empiricalModels/irradiances/FISM2/daily_bands_1947-2024.nc'
    myIrrTimesFISM2Bands, wavelengthsFISM2Bands, myDataAllFISM2Bands, _ = obtainFism2(fism2file, bands=True)
    myFluxDataAllFISM2Bands, myIrrDataAllFISM2Bands = myDataAllFISM2Bands

    # Replace bad values with NaNs:
    myIrrDataAllFISM2BandsFixed = myIrrDataAllFISM2Bands.copy()
    myIrrDataAllFISM2BandsFixed[myIrrDataAllFISM2BandsFixed <= 0] = np.nan

    # Harmonize the times for NEUVAC and FISM2:
    correspondingIndsFISM2StanBands = np.where((myIrrTimesFISM2Bands >= times[0]) & (myIrrTimesFISM2Bands <= times[-1]))[0]
    correspondingIrrTimesFISM2StanBands = myIrrTimesFISM2Bands[correspondingIndsFISM2StanBands]
    correspondingIrrFISM2StanBands = myIrrDataAllFISM2BandsFixed[correspondingIndsFISM2StanBands, :]

    # Mask the NaNs:
    cleanCorrespondingIrrFISM2StanBands = np.ma.masked_invalid(correspondingIrrFISM2StanBands)

    # Compute the normalized cross-correlation matrix between residuals in different bins:
    residualsArrayStanBands = np.subtract(neuvacIrrSolomon, cleanCorrespondingIrrFISM2StanBands[:, :-1])
    toolbox.savePickle(residualsArrayStanBands, 'residualsArrayStanBands.pkl')
    corMatStanBands = toolbox.mycorrelate2d(residualsArrayStanBands, normalized=True)
    toolbox.savePickle(corMatStanBands, 'corMatStanBands.pkl')

    # Compute the normalized standard deviation of NEUVAC irradiance residuals (in each band):
    STDNeuvacResidsStanBands = np.zeros(neuvacIrrSolomon.shape[1])
    for i in range(STDNeuvacResidsStanBands.shape[0]):
        STDNeuvacResidsStanBands[i] = np.nanstd(residualsArrayStanBands[:, i])
    # Save these values to be used later for running ensembles:
    toolbox.savePickle(STDNeuvacResidsStanBands, 'sigma_NEUVAC_StanBands.pkl')

    # View the correlation matrix for the residuals of the perturbed NEUVAC irradiances alongside the base NEUVAC irradiances:
    fig = plt.figure()
    pos = plt.imshow(corMatStanBands, aspect='auto', cmap='bwr', vmin=-1.0, vmax=1.0, interpolation='none')
    plt.xlabel('Wavelength Band')
    plt.ylabel('Wavelength Band')
    plt.title('Covariance Matrix (NEUVAC-22 - FISM2)')
    fig.colorbar(pos)
    plt.savefig('Uncertainty/corMatStanBands.png', dpi=300)

    # ------------------------------------------------------------------------------------------------------------------
    # 6: Compute the normalized cross-correlation matrices and normalized standard deviations for ALL OTHER MODELS
    # 6a: EUVAC
    euvacFlux, euvacIrr, _, _, _ = euvac.euvac(F107, F107A)

    residualsArrayEUVAC = np.subtract(euvacIrr, correspondingIrrFISM2[:, 7:44])
    toolbox.savePickle(residualsArrayEUVAC, 'residualsArrayEUVAC.pkl')
    corMatEUVAC = toolbox.mycorrelate2d(residualsArrayEUVAC, normalized=True)
    toolbox.savePickle(corMatEUVAC, 'corMatEUVAC.pkl')
    STDNeuvacResidsEUVAC = np.zeros(euvacIrr.shape[1])
    for i in range(STDNeuvacResidsEUVAC.shape[0]):
        STDNeuvacResidsEUVAC[i] = np.nanstd(residualsArrayEUVAC[:, i])
    toolbox.savePickle(STDNeuvacResidsEUVAC, 'sigma_EUVAC.pkl')
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    pos = ax.imshow(corMatEUVAC, aspect='auto', cmap='bwr', vmin=0.9999, vmax=1., interpolation='none')
    ax.set_xlabel('Wavelength Band', fontsize=16)
    ax.set_ylabel('Wavelength Band', fontsize=16)
    fig.suptitle('Correlation Matrix (EUVAC - FISM2)', fontsize=18)
    plt.colorbar(pos)
    fig.subplots_adjust(top=0.9)
    plt.savefig('Uncertainty/corMatsEUVAC.png', dpi=300)

    # 6b: HEUVAC
    # heuvac_wav, heuvacFlux, heuvacIrr, _, _, _ = heuvac.heuvac(F107, F107A, torr=True, statsFiles=['corMatHEUVAC.pkl',
    #                                                                                        'sigma_HEUVAC.pkl'])
    # residualsArrayHEUVAC = np.subtract(heuvacIrr, correspondingIrrFISM2[:, 7:44])
    # toolbox.savePickle(residualsArrayHEUVAC, 'residualsArrayHEUVAC.pkl')
    # corMatHEUVAC = toolbox.mycorrelate2d(residualsArrayHEUVAC, normalized=True)
    # toolbox.savePickle(corMatHEUVAC, 'corMatHEUVAC.pkl')
    # STDNeuvacResidsHEUVAC = np.zeros(heuvacIrr.shape[1])
    # for i in range(STDNeuvacResidsHEUVAC.shape[0]):
    #     STDNeuvacResidsHEUVAC[i] = np.nanstd(residualsArrayHEUVAC[:, i])
    # toolbox.savePickle(STDNeuvacResidsHEUVAC, 'sigma_HEUVAC.pkl')
    # fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    # pos = ax.imshow(corMatHEUVAC, aspect='auto', cmap='bwr', vmin=0.9999, vmax=1., interpolation='none')
    # ax.set_xlabel('Wavelength Band', fontsize=16)
    # ax.set_ylabel('Wavelength Band', fontsize=16)
    # fig.suptitle('Correlation Matrix (HEUVAC - FISM2)', fontsize=18)
    # plt.colorbar(pos)
    # fig.subplots_adjust(top=0.9)
    # plt.savefig('Uncertainty/corMatsHEUVAC.png', dpi=300)

    # 6c: TODO: SOLOMON (HFG and EUVAC) -- Optional

    # ------------------------------------------------------------------------------------------------------------------
    # Exit with a zero error code:
    sys.exit(0)