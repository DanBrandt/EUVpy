# This script performs the model fits for NEUVAC.

#-----------------------------------------------------------------------------------------------------------------------
# Top-level Imports:
import numpy as np
import matplotlib, sys
matplotlib.use('Qt5Agg')
from datetime import datetime
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Local Imports:
from src.EUVpy.NEUVAC import neuvac
from src.EUVpy.empiricalModels.models.EUVAC import euvac
from src.EUVpy.tools.EUV.fism2_process import read_euv_csv_file
from src.EUVpy.tools.processIrradiances import obtainFism2
from src.EUVpy.tools import toolbox, processIndices
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Directory Management
neuvac_directory = '../NEUVAC/'
neuvac_tableFile = '../NEUVAC/neuvac_table.txt'
neuvac_tableFile_Stan_Bands = '../NEUVAC/neuvac_table_stan_bands.txt'
figures_directory = 'Figures/'
results_directory = 'Results/'
fism1_spectra_folder = '../empiricalModels/irradiances/FISM1/'
fism2_spectra_folder = '../empiricalModels/irradiances/FISM2/'
euv_folder = '../tools/EUV/'
preparedDataFolder = '../experiments/preparedData'
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Constants
euvacTable = euvac.euvacTable
from src.EUVpy.empiricalModels.models.SOLOMON import solomon

solomonTable = solomon.solomonBands
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Execution
if __name__=="__main__":
    # Load in F10.7 data (Penticton, CA):
    pentictonTimesData = '../solarIndices/F107/Penticton/F107times.pkl'
    pentictonF107Data = '../solarIndices/F107/Penticton/F107vals.pkl'
    pentictonF107AveData = '../solarIndices/F107/Penticton/F107averageVals.pkl'
    pentictonTimes = toolbox.loadPickle(pentictonTimesData)
    pentictonF107 = toolbox.loadPickle(pentictonF107Data)
    pentictonF107A = toolbox.loadPickle(pentictonF107AveData)
    # F10.7 data extends between 1947-02-14; 12:00 to 2008-02-03; 12:00.
    # Load in F10.7 data (OMNIWeb):
    omniTimesData = '../solarIndices/F107/OMNIWeb/OMNIF107times.pkl'
    omniF107Data = '../solarIndices/F107/OMNIWeb/OMNIF107vals.pkl'
    omniF107AveData = '../solarIndices/F107/OMNIWeb/OMNIF107averageVals.pkl'
    omniTimes = toolbox.loadPickle(omniTimesData)
    omniF107 = toolbox.loadPickle(omniF107Data)
    omniF107B = toolbox.loadPickle(omniF107AveData)
    # F10.7 data extends between 1963-11-28; 12:00 to 2023-09-27; 12:00.
    # Load in CLS F10.7 data:
    CLS_times, CLS_F107, CLS_F107A, CLS_F107B = processIndices.getCLSF107('1951-11-01', '2024-08-15')
    # CLS data extends back to 1951-11-01; 12:00.

    times = CLS_times # omniTimes
    F107 = CLS_F107 # omniF107
    F107B = CLS_F107B # omniF107B

    euv_data_59 = read_euv_csv_file(euv_folder + 'euv_59.csv', band=False)
    mids = 0.5 * (euv_data_59['long'] + euv_data_59['short'])

    refit = [True, True] # Control whether not refitting is done for the GITM bins and the Solomon bins, respectively.
    # ==================================================================================================================
    # GITM BINS
    # ==================================================================================================================

    # FISM2 Results:
    fism2file = '../empiricalModels/irradiances/FISM2/daily_data_1947-2023.nc'
    myIrrTimesFISM2, wavelengthsFISM2, myIrrDataAllFISM2, myIrrUncAllFISM2 = obtainFism2(fism2file)
    # Rebin the data:
    myIrrDataWavelengthsFISM2, rebinnedIrrDataFISM2 = toolbox.newbins(wavelengthsFISM2, myIrrDataAllFISM2, euv_data_59, zero=True)

    # Replace bad values with NaNs:
    myIrrDataAllFISM2Fixed = rebinnedIrrDataFISM2.copy()
    myIrrDataAllFISM2Fixed[myIrrDataAllFISM2Fixed <= 0 ] = np.nan
    # FISM2 data extends between 1947-02-14; 00:00 and 2023-08-29; 00:00.

    if refit[0] == True:
        print('Fitting the NEUVAC base model...')
        # Perform a non-linear fit between F10.7, F10.7B, and FISM2:
        neuvacTable = neuvac.neuvacFit([times, F107, F107B], myIrrTimesFISM2, myIrrDataAllFISM2Fixed, wavelengths=mids, label='FISM2')

        # Print the coefficients to a file:
        with open(neuvac_directory+'neuvac_table.txt', 'w') as output:
            # Write the header information:
            output.write('This file contains coefficients for the current iteration of NEUVAC.\n'
                         'NOTE THAT THIS FILE CONTAINS COEFFICIENTS FOR THE GITM BINS.\n'
                         'This file was created on '+datetime.strftime(datetime.now(), '%Y-%m-%dT%H:%M:%S')+'\n'
                         'File Authors: Brandt, Daniel A. and Ridley, Aaron J.\n'
                         'This version of NEUVAC was created by fitting nonlinear models between F10.7 and preceding\n'
                         '54-day averaged F10.7 and FISM2 decomposed into 59 wavelength bands conventionally used in\n'
                         'the GITM and Aether models.\n'
                         'The file is formatted as follows:\n'
                         ' - First Column: Lower limit of the given wavelength bin in Angstroms.\n'
                         ' - Second Column: Upper limit of the given wavelength bin in Angstroms.\n'
                         ' - Third through Eighth colummns: Coefficients for the model.\n'
                         'The functional form of the model is given by:\n'
                         'Irr_i(t) = A_i * (F107(t) ** B_i) + C_i * (F107B(t) ** D_i) + E_i * (F107B(t) - F107(t)) + F_i\n'
                         'where the irradiance in bin i (Irr_i) is a function of time t, and A_i through F_i are \n'
                         'coefficients for bin i, and F107(t) and F107B(t) represent values of the F10.7 and 54-day\n'
                         'averaged F10.7 computed with a backwards-looking window, respectively.\n'
                         '-----------------------------------------------------------------------------------------------\n'
                         'WAVES WAVEL A_i B_i C_i D_i E_i F_i\n')
            for i in range(neuvacTable.shape[0]):
                output.writelines(str(euv_data_59['short'][i]) +' ' + str(euv_data_59['long'][i]) +' ' + toolbox.stringList(neuvacTable[i, :]) + '\n')
        print('Fitting of NEUVAC complete!')
        neuvacIrr, _, _, _ = neuvac.neuvacEUV(F107, F107B, bands=None, tableFile=neuvac_tableFile)

        # View the result of the model fits, as a sanity check:
        # for i in range(neuvacIrr.shape[1]):
        #     plt.figure()
        #     plt.plot(myIrrTimesFISM2, myIrrDataAllFISM2Fixed[:, i], label='FISM2')
        #     plt.plot(times, neuvacIrr[:, i], label='NEUVAC')
        #     plt.legend(loc='best')
        #     plt.title('Irradiance Time Series at :'+str(np.round(myIrrDataWavelengthsFISM2[i],2))+' Angstroms')

    # ==================================================================================================================
    # STANDARD BANDS
    # ==================================================================================================================

    # Load in FISM2 STAN BAND data:
    # FISM2 Stan Band Results:
    fism2file = '../empiricalModels/irradiances/FISM2/daily_bands_1947-2024.nc'
    myIrrTimesFISM2Bands, wavelengthsFISM2Bands, myDataAllFISM2Bands, _ = obtainFism2(fism2file, bands=True)
    myFluxDataAllFISM2Bands, myIrrDataAllFISM2Bands = myDataAllFISM2Bands

    # Replace bad values with NaNs:
    myIrrDataAllFISM2BandsFixed = myIrrDataAllFISM2Bands.copy()
    myIrrDataAllFISM2BandsFixed[myIrrDataAllFISM2BandsFixed <= 0] = np.nan

    if refit[1] == True:
        print('Fitting the NEUVAC-22 base model...')
        # Perform a non-linear fit between F10.7, F10.7B, and FISM2:
        midsS = (wavelengthsFISM2Bands * 10)[:-1]
        neuvacTableS = neuvac.neuvacFit([times, F107, F107B], myIrrTimesFISM2Bands, myIrrDataAllFISM2BandsFixed[:, :-1], wavelengths=midsS,
                                        label='FISM2S')

        # Print the coefficients to a file:
        with open(neuvac_directory+'neuvac_table_stan_bands.txt', 'w') as output:
            # Write the header information:
            output.write('This file contains coefficients for the current iteration of NEUVAC.\n'
                         'NOTE THAT THIS FILE CONTAINS COEFFICIENTS FOR THE STAN BANDS.\n'
                         'This file was created on '+datetime.strftime(datetime.now(), '%Y-%m-%dT%H:%M:%S')+'\n'
                         'File Authors: Brandt, Daniel A. and Ridley, Aaron J.\n'
                         'This version of NEUVAC was created by fitting nonlinear models between F10.7 and preceding\n'
                         '54-day averaged F10.7 and FISM2 decomposed into 23 wavelength bands conventionally used in\n'
                         'numerous atmospheric models, such as SAMI3.\n'
                         'The file is formatted as follows:\n'
                         ' - First Column: Lower limit of the given wavelength bin in Angstroms.\n'
                         ' - Second Column: Upper limit of the given wavelength bin in Angstroms.\n'
                         ' - Third through Eighth colummns: Coefficients for the model.\n'
                         'The functional form of the model is given by:\n'
                         'Irr_i(t) = A_i * (F107(t) ** B_i) + C_i * (F107B(t) ** D_i) + E_i * (F107B(t) - F107(t)) + F_i\n'
                         'where the irradiance in bin i (Irr_i) is a function of time t, and A_i through F_i are \n'
                         'coefficients for bin i, and F107(t) and F107B(t) represent values of the F10.7 and 54-day\n'
                         'averaged F10.7 computed with a backwards-looking window, respectively.\n'
                         '-----------------------------------------------------------------------------------------------\n'
                         'WAVES WAVEL A_i B_i C_i D_i E_i F_i\n')
            for i in range(neuvacTableS.shape[0]):
                output.writelines(str(solomonTable['short'][i]) +' ' + str(solomonTable['long'][i]) +' ' + toolbox.stringList(neuvacTableS[i, :]) + '\n')

        neuvacIrrS, _, _, _ = neuvac.neuvacEUV(F107, F107B, bands='SOLOMON', tableFile=neuvac_tableFile_Stan_Bands)
        # View the result of the model fits, as a sanity check:
        # for i in range(neuvacIrrS.shape[1]):
        #     plt.figure()
        #     plt.plot(myIrrTimesFISM2Bands, myIrrDataAllFISM2BandsFixed[:, i], label='FISM2S')
        #     plt.plot(times, neuvacIrrS[:, i], label='NEUVAC')
        #     plt.legend(loc='best')
        #     plt.title('Irradiance Time Series at :'+str(np.round(midsS[i],2))+' Angstroms')
        print('Fitting of NEUVAC-22 complete!')

    sys.exit(0)
#-----------------------------------------------------------------------------------------------------------------------
