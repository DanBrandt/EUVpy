# This module contains functions used for processing/loading in solar irradiances.

#-----------------------------------------------------------------------------------------------------------------------
# Top-level Imports
import os, sys

import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np
import pickle
from datetime import datetime, timedelta
import matplotlib
matplotlib.use('Qt5Agg')
from netCDF4 import Dataset
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Local Imports
from tools import toolbox
from tools import spectralAnalysis
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Directory Management
fism1_spectra_folder = '../empiricalModels/irradiances/FISM1/'
fism2_spectra_folder = '../empiricalModels/irradiances/FISM2/'
TIMED_spectra_folder = '../measurements/TIMED_SEE_Level_3/'
euv_folder = '../tools/EUV/'
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Global Variable:
SEEBands = np.array([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5,
         9.5,  10.5,  11.5,  12.5,  13.5,  14.5,  15.5,  16.5,  17.5,
        18.5,  19.5,  20.5,  21.5,  22.5,  23.5,  24.5,  25.5,  26.5,
        27.5,  28.5,  29.5,  30.5,  31.5,  32.5,  33.5,  34.5,  35.5,
        36.5,  37.5,  38.5,  39.5,  40.5,  41.5,  42.5,  43.5,  44.5,
        45.5,  46.5,  47.5,  48.5,  49.5,  50.5,  51.5,  52.5,  53.5,
        54.5,  55.5,  56.5,  57.5,  58.5,  59.5,  60.5,  61.5,  62.5,
        63.5,  64.5,  65.5,  66.5,  67.5,  68.5,  69.5,  70.5,  71.5,
        72.5,  73.5,  74.5,  75.5,  76.5,  77.5,  78.5,  79.5,  80.5,
        81.5,  82.5,  83.5,  84.5,  85.5,  86.5,  87.5,  88.5,  89.5,
        90.5,  91.5,  92.5,  93.5,  94.5,  95.5,  96.5,  97.5,  98.5,
        99.5, 100.5, 101.5, 102.5, 103.5, 104.5, 105.5, 106.5, 107.5,
       108.5, 109.5, 110.5, 111.5, 112.5, 113.5, 114.5, 115.5, 116.5,
       117.5, 118.5, 119.5, 120.5, 121.5, 122.5, 123.5, 124.5, 125.5,
       126.5, 127.5, 128.5, 129.5, 130.5, 131.5, 132.5, 133.5, 134.5,
       135.5, 136.5, 137.5, 138.5, 139.5, 140.5, 141.5, 142.5, 143.5,
       144.5, 145.5, 146.5, 147.5, 148.5, 149.5, 150.5, 151.5, 152.5,
       153.5, 154.5, 155.5, 156.5, 157.5, 158.5, 159.5, 160.5, 161.5,
       162.5, 163.5, 164.5, 165.5, 166.5, 167.5, 168.5, 169.5, 170.5,
       171.5, 172.5, 173.5, 174.5, 175.5, 176.5, 177.5, 178.5, 179.5,
       180.5, 181.5, 182.5, 183.5, 184.5, 185.5, 186.5, 187.5, 188.5,
       189.5, 190.5, 191.5, 192.5, 193.5, 194.5]) * 10
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
def getIrr(dateStart, dateEnd, source):
    """
    Given a starting date and an ending date, automatically download irradiance data from LISIRD for a specific source,
    including FISM2 (daily or stan bands) or SEE (Level 3 daily).
    :param dateStart: str
        The starting date for the data in YYYY-MM-DD format.
    :param dateEnd: str
        The ending date for the data in YYYY-MM-DD format.
    :param source: str
        The type of data to be obtained. Valid inputs are:
        - FISM2 (for daily averages of FISM2 data)
        - FISM2S (for daily averages of FISM2 standard bands, according to Solomon and Qian 2005)
        - SEE (for Level 3 daily averages of TIMED/SEE data)
    :return times: ndarray
        Datetime values for each spectrum.
    :return wavelengths: ndarray
        Wavelength bins (bin boundaries) for the spectral data.
    :return irradiance: ndarray
        A 2D array where each row is a spectrum at a particular time, and the columns are wavelength bands.
    """
    # Converting the input time strings to datetimes:
    dateStartDatetime = datetime.strptime(dateStart, "%Y-%m-%d")
    dateEndDatetime = datetime.strptime(dateEnd, "%Y-%m-%d")

    # Check if the user has asked for a source that can be obtained:
    validSources = ['FISM2', 'FISM2S', 'SEE']
    if source not in validSources:
        raise ValueError("Variable 'source' must be either 'FISM2', 'FISM2S', or 'SEE'.")

    # Download the most recent file for the corresponding source and read it in:
    if source == 'FISM2':
        url = 'https://lasp.colorado.edu/eve/data_access/eve_data/fism/daily_hr_data/daily_data.nc'
        fname = 'FISM2_daily_data.nc'
        toolbox.urlObtain(url, fism2_spectra_folder + fname)
        datetimes, wavelengths, irradiance, uncertainties = obtainFism2(fism2_spectra_folder + fname)
    elif source == 'FISM2S':
        url = 'https://lasp.colorado.edu/eve/data_access/eve_data/fism/daily_bands/daily_bands.nc'
        fname = 'FISM2_daily_bands.nc'
        toolbox.urlObtain(url, fism2_spectra_folder + fname)
        datetimes, wavelengths, irradiance, uncertainties = obtainFism2(fism2_spectra_folder + fname, bands=True)
    else:
        url = 'https://lasp.colorado.edu/data/timed_see/level3/latest_see_L3_merged.ncdf'
        fname = 'TIMED_SEE_Level_3.nc'
        toolbox.urlObtain(url, TIMED_spectra_folder + fname)
        datetimes, wavelengths, irradiance, uncertainties = obtainSEE(TIMED_spectra_folder + fname)

    # Subset the data according to user demands:
    validInds = np.where((datetimes >= dateStartDatetime) & (datetimes <= dateEndDatetime))[0]
    times = datetimes[validInds]
    irradiance = irradiance[validInds, :]

    # Return the resulting data:
    return times, wavelengths, irradiance

def obtainFism1(fismFiles, euv_bins, saveLoc=None):
    """
    Given muliple FISM1 .dat files, get the information from each band, using code developed by Dr. Aaron Ridley.
    :param fismfiles: arraylile
        A list or array of .dat files containing FISM1 data (the full VUV spectrum from .1nm to 195nm at 1 nm
        resolution.
    :param euv_bins: dict
        EUV bins with which to rebin the FISM1 data. Obtained from fism2_process.rebin_fism.
    :param saveLoc: str
        Optional argument that controls where pickle files are saved.
    :return irrTimes: ndarray
        A 1d array of datetimes corresponding to each set of irradiance values.
    :return irrArray: ndarray
        An ndarray containing all of the individual 59 irradiance values in each band from all .dat files.
    """
    myIrrPickleFile = 'myIrrFISM1.pkl'
    myTimePickleFile = 'myTimesFISM1.pkl'
    override = True # If true, forcefully read in the CSV data irrespective if it is previously existing.
    irrTimes = [] #np.zeros(len(fismFiles))
    numBins = len(euv_bins['long'])
    if saveLoc != None:
        searchString = saveLoc + myIrrPickleFile
    else:
        searchString = myIrrPickleFile
    if os.path.isfile(searchString) == False or override == True:
        irrArray = []
        for i in tqdm(range(len(fismFiles))): # Loop through the files.
            with open(fism1_spectra_folder+fismFiles[i]) as fismFileInfo:
                fismFileData = fismFileInfo.readlines()
                currentIrrArray = np.zeros((len(fismFileData[1:]), numBins))
                j = -1
                for line in fismFileData: # Loop over the lines in the file.
                    if j >= 0:
                        fismLineData = line.split()
                        # The first six elements of fismLineData are part of the time stamp:
                        irrTimes.append(datetime(int(fismLineData[0]), int(fismLineData[1]), int(fismLineData[2]), int(fismLineData[3])))
                        # The remaining elements of fismLineData are the irradiance in 59 wavelengths (the order must be flipped):
                        currentIrrArray[j, :] = np.flip(np.asarray(fismLineData[6:]))
                    j += 1
                irrArray.append(currentIrrArray)
        finalIrrArray = np.concatenate(irrArray)
        # Sort the data and order it properly:
        sort_indices = np.argsort(irrTimes)
        irrTimes = np.asarray(irrTimes)[sort_indices]
        for k in range(numBins):
            finalIrrArray[:, k] = finalIrrArray[:, k][sort_indices]
        # Write to a pickle (since the data takes a while to gather and read in):
        if saveLoc != None:
            myTimePkl = open(myTimePickleFile, 'wb')
            myIrrPkl = open(myIrrPickleFile, 'wb')
        else:
            myTimePkl = open(saveLoc + myTimePickleFile, 'wb')
            myIrrPkl = open(saveLoc + myIrrPickleFile, 'wb')
        pickle.dump(np.asarray(irrTimes), myTimePkl)
        pickle.dump(irrArray, myIrrPkl)
    else:
        myTimePkl = open(saveLoc+myTimePickleFile, 'rb')
        irrTimes = pickle.load(myTimePkl)
        myIrrPkl = open(saveLoc+myIrrPickleFile, 'rb')
        finalIrrArray = pickle.load(myIrrPkl)
    return irrTimes, finalIrrArray

def obtainFism2(myFism2File, bands=False):
    """
    Load in spectrum data from a FISM2 file.
    :param myFism2File: str
        The location of the NETCDF4 file.
    :param bands: bool
        If True, loads in the data segmented into the Solomon and Qian 2005 standard bands.
    :return datetimes: ndarray
        An array of datetimes for each TIMED/SEE spectra.
    :return wavelengths: ndarray
        A one-dimensional array of wavelengths at which there are irradiance values.
    :return irradiances: ndarray
        A two-dimensional array of irradiance values at each time.
    :return uncertainties: ndarray
        A two-dimensional array of irradiance uncertainty values at each time.
    """
    fism2Data = Dataset(myFism2File)
    wavelengths = np.asarray(fism2Data.variables['wavelength'])
    if bands==True: # STANDARD BANDS
        flux = np.asarray(fism2Data.variables['ssi']) # photons/cm2/second
        # bandwidths = np.asarray(fism2Data.variables['band_width'])
        pFlux = flux * 1.0e4 # photons/m2/second
        # Convert fluxes to irradiances:
        irr = np.zeros_like(flux)
        for i in range(flux.shape[1]):
            irr[:, i] = spectralAnalysis.spectralIrradiance(pFlux[:, i], wavelengths[i]*10.) # W/m^2
        irradiance = [flux, irr]
        uncertainties = np.asarray(fism2Data.variables['band_width']) # TODO: Replace with an estimation of uncertainty
    else: # NATIVE DATA
        irradiance = np.asarray(fism2Data.variables['irradiance']) # W/m^2/nm
        uncertainties = np.asarray(fism2Data.variables['uncertainty'])
    dates = fism2Data.variables['date']
    datetimes = []
    for i in range(len(dates)):
        year = dates[i][:4]
        day = dates[i][4:]
        currentDatetime = datetime(int(year), 1, 1) + timedelta(int(day) - 1) + timedelta(hours=12)
        datetimes.append(currentDatetime)
    datetimes = np.asarray(datetimes)
    return datetimes, wavelengths, irradiance, uncertainties

def obtainSEE(seeFile):
    """
    Given a TIMED/SEE NETCDF4 file, load in and return the timestamps, wavelengths, irradiances, and uncertainties.
    :param seeFile: str
        The NETCDF4 file containing TIMED/SEE data.
    :return datetimes: ndarray
        An array of datetimes for each TIMED/SEE spectra.
    :return wavelengths: ndarray
        A one-dimensional array of wavelengths at which there are irradiance values.
    :return irradiances: ndarray
        A two-dimensional array of irradiance values at each time.
    :return uncertainties: ndarray
        A two-dimensional array of irradiance uncertainty values at each time.
    """
    seeData = Dataset(seeFile)
    dates = np.squeeze(seeData.variables['DATE'])
    wavelengths = np.squeeze(seeData.variables['SP_WAVE'])
    irradiances = np.squeeze(seeData.variables['SP_FLUX'])
    uncertainties = np.squeeze(seeData.variables['SP_ERR_TOT'])
    precision = np.squeeze(seeData.variables['SP_ERR_MEAS'])
    datetimes = []
    for i in range(len(dates)):
        year = str(dates[i])[:4]
        day = str(dates[i])[4:]
        currentDatetime = datetime(int(year), 1, 1) + timedelta(int(day) - 1) + timedelta(hours=12)
        datetimes.append(currentDatetime)
    datetimes = np.asarray(datetimes)
    return datetimes, wavelengths, irradiances, uncertainties

def obtainNRLSSIS2(filename):
    """
    Given a NRLSSI2 NETCDF4 file, load in and return the timestamps, wavelengths, irradiances, and uncertainties.
    :param filename: str
        The NETCDF4 file containing TIMED/SEE data.
    :return datetimes: ndarray
        An array of datetimes for each TIMED/SEE spectra.
    :return wavelengths: ndarray
        A one-dimensional array of wavelengths at which there are irradiance values.
    :return bandwidths: ndarray
        The width of the corresponding band for each irradiance measurement.
    :return irradiance: ndarray
        A two-dimensional array of irradiance values at each time.
    :return uncertainties: ndarray
        A two-dimensional array of irradiance uncertainty values at each time.
    """
    NRLData = Dataset(filename)
    dates = np.squeeze(NRLData.variables['time'])
    irradiances = np.squeeze(NRLData.variables['SSI'])
    wavelengths = np.squeeze(NRLData.variables['wavelength'])
    bandwidths = np.squeeze(NRLData.variables['Wavelength_Band_Width'])
    uncertainties = np.squeeze(NRLData.variables['SSI_UNC'])
    startingEpoch = datetime(1610, 1, 1)
    datetimes = []
    for i in range(len(dates)):
        currentDatetime = startingEpoch + timedelta(days=int(dates[i]), hours=12) # Move the observations to coincide with noon of each day.
        datetimes.append(currentDatetime)
    datetimes = np.asarray(datetimes)
    return datetimes, wavelengths, bandwidths, irradiances, uncertainties

def obtainSDO(seeFile):
    """
    Given an SDO NETCDF4 file, load in and return the timestamps, wavelengths, irradiances, and uncertainties.
    :param seeFile: str
        The NETCDF4 file containing SDO/EVE data.
    :return datetimes: ndarray
        An array of datetimes for each SDO/EVE spectra.
    :return wavelengths: ndarray
        A one-dimensional array of wavelengths at which there are irradiance values.
    :return irradiances: ndarray
        A two-dimensional array of irradiance values at each time.
    :return uncertainties: ndarray
        A two-dimensional array of irradiance uncertainty values at each time.
    """
    sdoData = Dataset(seeFile)
    dates = np.squeeze(sdoData.variables['MERGEDDATA.YYYYDOY'])
    wavelengths = np.squeeze(sdoData.variables['SPECTRUMMETA.WAVELENGTH'])
    irradiances = np.squeeze(sdoData.variables['MERGEDDATA.SP_IRRADIANCE'])
    uncertainties = np.squeeze(sdoData.variables['MERGEDDATA.SP_STDEV'])
    precision = np.squeeze(sdoData.variables['MERGEDDATA.SP_PRECISION'])
    datetimes = []
    for i in range(len(dates)):
        year = str(dates[i])[:4]
        day = str(dates[i])[4:]
        currentDatetime = datetime(int(year), 1, 1) + timedelta(int(day) - 1) + timedelta(hours=12)
        datetimes.append(currentDatetime)
    datetimes = np.asarray(datetimes)
    return datetimes, wavelengths, irradiances, uncertainties
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Execution
# if __name__=="__main__":
    # euv_data_59 = read_euv_csv_file(euv_folder + 'euv_59.csv', band=False)
    # mids = 0.5*(euv_data_59['short'] + euv_data_59['long'])
    #
    # # Load in the FISM1 data in multiple files, rebin it into 59 GITM wavelength bins, and combine it into a single pickle file:
    # # myFism1Files = os.listdir(fism1_spectra_folder)
    # # myIrrTimesFISM1, myIrrDataAllFISM1 = obtainFism1(myFism1Files, euv_data_59, saveLoc=fism1_spectra_folder)
    #
    # # Load in the FISM2 data, rebin it into 59 GITM wavelength bins, and save it to a pickle file:
    # fism2file = '../empiricalModels/irradiances/FISM2/daily_data_1947-2023.nc'
    # myIrrTimesFISM2, myFISM2Wavelengths, myIrrDataAllFISM2, myIrrUncAllFISM2 = obtainFism2(fism2file)
    # # Rebinning FISM2 into the SEE Bands:
    # rebinnedFISM2Wavelengths_like_SEE, rebinnedFISM2Irr_like_SEE = toolbox.rebin(wavelengths=myFISM2Wavelengths,
    #                                                                      data=myIrrDataAllFISM2,
    #                                                                      limits=[np.min(SEEBands)/10., np.max(SEEBands)/10.], resolution=1.)
    # # Then rebinning it into the GITM Bands:
    # rebinnedFISM2Wavelengths, rebinnedFISM2Irr = toolbox.rebin(wavelengths=myFISM2Wavelengths, data=myIrrDataAllFISM2,
    #                                                    resolution=euv_data_59, zero=False) #limits=[np.min(SEEBands)/10., np.max(SEEBands)/10.])
    #
    # # Load in SDO/EVE data:
    # # sdoFile = '../measurements/SDO_EVE_Level_3/latest_EVE_L3_merged_1nm_2010-2023.ncdf'
    # # myIrrTimesSDO, mySDOWavelengths, myIrrDataAllSDO, myIrrUncAllSDO = obtainSDO(sdoFile)
    # # rebinnedSDOWavelengths, rebinnedSDOIrr = rebin(wavelengths=mySDOWavelengths, data=myIrrDataAllSDO[20166, :],
    # #                                                    resolution=1.,
    # #                                                    limits=[np.min(SEEBands) / 10., np.max(SEEBands) / 10.])
    # # FISM2 20066 = TIMED/SEE 0
    #
    # # Load in TIMED/SEE data and rebin it into 59 GITM wavelength bins:
    # seeFile = '../measurements/TIMED_SEE_Level_3/latest_see_L3_merged_2002-2023.ncdf'
    # myIrrTimesSEE, mySEEWavelengths, myIrrDataAllSEE, myIrrUncAllSEE = obtainSEE(seeFile)
    # myCleanIrrDataAllSEE = myIrrDataAllSEE.copy()
    # myCleanIrrDataAllSEE[myCleanIrrDataAllSEE <= 0] = np.nan
    # rebinnedSEEWavelengths, rebinnedSEEIrr = toolbox.rebin(wavelengths=mySEEWavelengths, data=myCleanIrrDataAllSEE,
    #                                                    resolution=euv_data_59, factor=5, zero=False)
    #
    # # Load in NRLSSI2 data, rebin it into 59 GITM wavelength bins, and save it to pickle file:
    # # NRLFile = '../empiricalModels/irradiances/NRLSSI2/ssi_v02r01_daily_s18820101_e20221231_c20230123.nc'
    # # datetimesNRL, wavelengthsNRL, bandwidthsNRL, irradiancesNRL, uncertaintiesNRL = obtainNRLSSIS2(NRLFile)
    #
    # # VISUALIZE FISM2 and SEE:
    # import matplotlib.pyplot as plt
    # plt.figure()
    # # plt.plot(myFISM2Wavelengths, myIrrDataAllFISM2[20166, :], label='FISM2 (raw)')
    # # plt.plot(rebinnedFISM2Wavelengths_like_SEE, rebinnedFISM2Irr_like_SEE[20166, :], label='FISM2 (rebinned like SEE)')
    # plt.plot(rebinnedFISM2Wavelengths, rebinnedFISM2Irr[20166, :], marker='o', label='FISM2 (rebinned)')
    # # plt.plot(mySEEWavelengths, myIrrDataAllSEE[100, :], label='TIMED/SEE (raw)')
    # plt.plot(rebinnedSEEWavelengths, rebinnedSEEIrr[100, :], marker='o', label='TIMED/SEE (rebinned)')
    # # Plot the bin boundaries:
    # singularInds = np.where(euv_data_59['short'] == euv_data_59['long'])[0]
    # widthInds = np.where(euv_data_59['short'] != euv_data_59['long'])[0]
    # # for line in euv_data_59['short'][singularInds]/10.:
    # #     plt.axvline(x=line, color='b')
    # # for line2 in euv_data_59['short'][widthInds]/10.:
    # #     plt.axvline(x=line2, color='k')
    # # for line3 in euv_data_59['long'][widthInds]/10.:
    # #     plt.axvline(x=line3, color='k')
    # plt.legend()
    # plt.ylim([-0.0017, 0.0157])
    # # Semilogy:
    # plt.figure()
    # ax = plt.gca()
    # plt.scatter(myFISM2Wavelengths, myIrrDataAllFISM2[0, :], label='FISM2 (raw)')
    # plt.scatter(rebinnedFISM2Wavelengths, rebinnedFISM2Irr[20166, :], label='FISM2 (rebinned)')
    # plt.scatter(mySEEWavelengths, myIrrDataAllSEE[100, :], label='TIMED/SEE (raw)')
    # plt.scatter(rebinnedSEEWavelengths, rebinnedSEEIrr[100, :], label='TIMED/SEE (rebinned)')
    # ax.set_yscale('log')
    # plt.legend()
    #
    # # Time series:
    # for i in range(rebinnedFISM2Irr.shape[1]):
    #     # i = 2
    #     plt.figure();
    #     plt.plot(myIrrTimesFISM2, rebinnedFISM2Irr[:, i], label='FISM2')
    #     plt.plot(myIrrTimesSEE, rebinnedSEEIrr[:, i], label='TIMED/SEE')
    #     plt.legend()
    #     plt.title(str(mids[i])+' Angstroms')
    #
    # sys.exit(0)
# -----------------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    datetimes, wavelengths, irradiances, uncertainties = obtainSEE(TIMED_spectra_folder+'latest_see_L3_merged_2002-2023.ncdf')
    datetimesF, wavelengthsF, irradiancesF, uncertaintiesF = obtainFism2(fism2_spectra_folder+'daily_data_1947-2023.nc')

    import matplotlib.pylab as pylab
    params = {'legend.fontsize': 'x-large',
              'figure.figsize': (14, 6),
              'axes.labelsize': 'x-large',
              'axes.titlesize': 'X-Large',
              'xtick.labelsize': 'X-large',
              'ytick.labelsize': 'X-large'}
    sampleVals1 = irradiancesF[21056, :]
    idx = np.flatnonzero(~np.isnan(sampleVals1))
    sampleVals1[idx[sampleVals1[idx] < 0]] = np.nan
    sampleVals2 = irradiances[1000, :]
    idx = np.flatnonzero(~np.isnan(sampleVals2))
    sampleVals2[idx[sampleVals2[idx] < 0]] = np.nan
    pylab.rcParams.update(params)
    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(wavelengthsF*10, sampleVals1, color='b', label='FISM2')
    plt.plot(wavelengths*10, sampleVals2, color='r', label='TIMED/SEE')
    plt.xlabel('Wavelength (Angstroms)')
    plt.ylabel('Irradiance (W/m$^2$/nm)')
    plt.title('Solar EUV Spectra on '+str(datetimes[1000].date()))
    plt.legend(loc='best')
    plt.yscale('log')

    # TODO: Fix plot axis labels (fontsize), add a title, etc.

    ellipsis