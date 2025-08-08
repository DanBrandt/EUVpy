# This module contains functions used for processing/loading in solar indices.

#-----------------------------------------------------------------------------------------------------------------------
# Top-level Imports
import os
from datetime import datetime, timedelta
import pandas as pd
import numpy as np
import sys, csaps
import matplotlib
matplotlib.use('Qt5Agg')
import urllib.request
from pathlib import Path
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Local imports:
from EUVpy.tools.toolbox import uniformSample, imputeData, rollingAverage, rollingStd, savePickle
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Folder for downloading F10.7 data:
F107Folder = '../solarIndices/F107/OMNIWeb/'
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Functions:
def readF107(filename):
    """
    Load in F10.7 data from the files obtained from the FTP server provided by the Government of Canada:
    https://www.spaceweather.gc.ca/forecast-prevision/solar-solaire/solarflux/sx-5-en.php
    Requires ONLY a filename containing F10.7 data.
    Note that F10.7 values are daily values measured at LOCAL NOON of each day.
    Note that the OBSERVED values DO NOT correspond to 1 AU, and they may vary, while the ADJUSTED values are calibrated
    to correspond to 1 AU.

    Parameters
    ----------
    filename : str
        A string containing the location of the txt file with the F10.7 data.

    Returns
    -------
    times : numpy.ndarray
        A 1D array of datetimes corresponding to each F10.7 measurement.
    f107 : numpy.ndarray
        A 1D array of F10.7 values.
    """
    times = []
    f107 = []
    i = 0
    # Loop through the file to collect the data:
    with open(filename, 'r') as myFile:
        for line in myFile:
            # Note that the line structure is different depending on what file is downloaded:
            if filename[-11:-4] != 'current':
                if i > 2:
                    if filename[-8:-4] == '2007':
                        elements = line.split()
                        try:
                            f107.append(float(elements[7])) # Index 6 = Observed; Index 7 = Adjusted
                        except:
                            f107.append(np.nan)
                        currentTime = pd.to_datetime(float(elements[0]), origin='julian', unit='D')
                    else:
                        elements = line.split(',')
                        try:
                            f107.append(float(elements[6])) # Index 5 = Observed; Index 6 = Adjusted
                        except:
                            elements = line.split(' ')
                            f107.append(float(elements[-1]))
                        currentTime = pd.to_datetime(float(line.split(',')[0]), origin='julian', unit='D')
                    times.append(currentTime.to_pydatetime())
            else:
                if i > 3:
                    elements = line.split()
                    f107.append(float(elements[8])) # Index 8 = Adjusted
                    currentTime = pd.to_datetime(float(elements[0]), origin='julian', unit='D')
                    times.append(currentTime.to_pydatetime())
            i += 1

    return np.asarray(times), np.asarray(f107)

def readCLS(filename):
    """
    Load in flare-corrected, Sun-Earth distance adjusted flux values recorded by the Collecte Localisation Satellites
    (CLS).

    Parameters
    ----------
    filename : str
        The location of the data file.

    Returns
    -------
    times : list
        The datetimes for each data value.
    data : numpy.ndarray
        The solar flux data for F30, F15, F10.7, F8, and F3.2.
    """
    times = []
    precisionVals = []
    with open(filename, 'r') as myFile:
        allLines = myFile.readlines()
        data = np.zeros((len(allLines)-25, 5))
        i = 0
        j = 0
        for line in allLines:
            if i >= 25:
                elements = line.split()
                try:
                    data[j, :] = np.array([float(elements[5]), float(elements[9]), float(elements[13]), float(elements[17]), float(elements[21])])
                    times.append( datetime(int(elements[0]), int(elements[1]), int(elements[2]), 12) )
                    precisionVals.append( [float(elements[6]), float(elements[10]), float(elements[14]), float(elements[18]), float(elements[22])] )
                except:
                    raise Exception
                j += 1
            i += 1
    # Print the precision:
    # print('Mean precision values...')
    # print('F30: '+str(np.nanmean([element[0] for element in precisionVals]))+' sfu') # 6
    # print('F15: ' + str(np.nanmean([element[1] for element in precisionVals])) + ' sfu') # 8
    # print('F10.7: ' + str(np.nanmean([element[2] for element in precisionVals])) + ' sfu') # 13
    # print('F8: ' + str(np.nanmean([element[3] for element in precisionVals])) + ' sfu') # 12
    # print('F3.2: ' + str(np.nanmean([element[4] for element in precisionVals])) + ' sfu') # 11
    return times, data

def getCLSF107(dateStart, dateEnd, truncate=True, rewrite=True):
    """
    Obtains Sun-Earth distance adjusted, flare-corrected F10.7 data from Collecte Localisation Satellites. Downloads the
    most recent measurements to a file. Reads the file and extracts the F10.7 values between two dates. Note that if the
    ending date is less than or equal to the last date in the version of the file that has already been downloaded, the
    file IS NOT re-downloaded, but simply parsed. Otherwise, the file is redownloaded.

    Parameters
    ----------
    dateStart : str
        The starting date in YYYY-MM-DD format.
    dateEnd : str
        The ending date in YYYY-MM-DD format.
    truncate : bool
        Controls whether to truncate the data to exclude the most recent 81 days. Defaults is True.
    rewrite : bool
        Controls whether an existing CLS file is rewritten. Defualt is True.

    Returns
    -------
    times : list
        The datetimes for each data value.
    F107 : arraylike
        Solar flux at 10.7 cm.
    F107A : arraylike
        81-day averaged solar flux at 10.7 cm, centered on the current day.
    F107B : arraylike
        54-day averaged solar flux at 10.7 cm, averaged in a backwards-looking window.
    """
    dateTimeStart = datetime.strptime(dateStart, '%Y-%m-%d')
    dateTimeEnd = datetime.strptime(dateEnd, '%Y-%m-%d')
    euvpy_app_folder = Path.home().joinpath(".euvpy")
    if not euvpy_app_folder.exists():
        euvpy_app_folder.mkdir(exist_ok=True)
    # 1: Check if there is ALREADY a F10.7 file present:
    #fname = '../solarIndices/F107/radio_flux_adjusted_observation.txt'
    fname = euvpy_app_folder.joinpath("radio_flux_adjusted_observation.txt")
    if fname.exists() and rewrite == False:
        # Read in the file:
        times, data = readCLS(fname)
        # Check if the ending date exceeds the ending date in the file. If so, redownloading the file:
        if times[-1] > dateTimeEnd:
            out = urllib.request.urlretrieve(
                'ftp://ftpsedr.cls.fr/pub/previsol/solarflux/observation/radio_flux_adjusted_observation.txt', fname)
        times, data = readCLS(fname)
    else:
        # Download the file:
        out = urllib.request.urlretrieve('ftp://ftpsedr.cls.fr/pub/previsol/solarflux/observation/radio_flux_adjusted_observation.txt', fname)
        times, data = readCLS(fname)

    # Compute the 81-day (centered) averaged F10.7 and 54-day averaged (:
    F107 = data[:, 2]
    F107A = rollingAverage(F107, window_length=81, impute_edges=True)
    F107B = rollingAverage(F107, window_length=54, impute_edges=True, center=False)
    print(f'\n\nlen(F107) = {len(F107)}, len(F107B) = {len(F107B)}, len(F107B) = {len(F107B)}')
    # Extract the values in the desired time range:
    goodInds = np.where((np.asarray(times) >= dateTimeStart) & (np.asarray(times) <= dateTimeEnd))[0]
    # Truncation:
    if truncate and len(goodInds) >= 2*81:
        goodInds = goodInds[:-81]
    return np.asarray(times)[goodInds], np.asarray(F107)[goodInds], np.asarray(F107A)[goodInds], np.asarray(F107B)[goodInds]

def getF107(dateStart, dateEnd):
    """
    Given two dates (a start date and an ending date), automatically download F10.7 data from NASA OMNIWeb.

    Parameters
    ----------
    dateStart : str
        The starting date, in YYYY-MM-DD format.
    dateEnd : str
        The ending date, in YYYY-MM-DD format.

    Returns
    -------
    dataFile : str
        The downloaded OMNI data file.
    times : numpy.ndarray
        An array of datetimes for the F10.7 values.
    f107 : numpy.ndarray
        The F10.7 values for the time desired.
    f107A : str
        The 81-day averaged F10.7 values (centered on the current day) for the time desired.
    """
    dateStartStr = dateStart.replace('-','')
    dateEndStr = dateEnd.replace('-','')
    fileStr = F107Folder+'f107_'+dateStartStr+'_'+dateEndStr+'.txt'
    if os.path.isfile(fileStr) == False:
        # Set up the command string:
        cmd = 'wget --post-data "activity=retrieve&res=hour&spacecraft=omni2&start_date='+dateStartStr+'&end_date='+\
              dateEndStr+'&vars=50&scale=Linear&ymin=&ymax=&view=0&charsize=&xstyle=0&ystyle=0&symbol=0&symsize=&' \
              'linestyle=solid&table=0&imagex=640&imagey=480&color=&back=" ' \
              'https://omniweb.gsfc.nasa.gov/cgi/nx1.cgi -O '+fileStr
        # Execute the command:
        os.system(cmd)
    #dataFile = fileStr

    # Obtain the data:
    times = []
    f107 = []
    i = 0
    with open(fileStr, 'r') as f107File:
        f107Data = f107File.readlines()
        for line in f107Data:
            if i >= 8 and i < len(f107Data) - 15:
                lineData = line.split()
                times.append(datetime(int(lineData[0]), 1, 1) + timedelta(days=int(lineData[1]) - 1) + timedelta(
                hours=int(lineData[2])))
                f107.append(float(lineData[-1]))
            i += 1
    # Convert to arrays:
    times = np.asarray(times)
    f107 = np.asarray(f107)
    # Clean the data:
    cleanedF107 = cleanF107(times, f107, bad_value=999.9)
    filteredF107 = F107filter(times, cleanedF107)
    f107A = rollingAverage(filteredF107, window_length=81)
    # Return the data:
    return times, filteredF107, f107A

def readOMNI(dataFile, headerFile):
    """
    Read in file from an OMNI data and parse it according to the specified format.

    Parameters
    ----------
    dataFile : str
        The name of the file where the OMNI data is stored.
    headerFile : str
        The name of the file containing header information for the OMNI data.

    Returns
    -------
    omniTimes : numpy.ndarray
        A 1D array of datetimes for the omni data.
    omniLabels: numpy.ndarray
        A 1D array of strings of each of the variables in the OMNI data.
    omniDataArray: numpy.ndarray
        A 2D array with all the OMNI data. The shape is nxm, where n is the number of time samples (each hour) and
        m is the number of variables collected at each time sample.
    """
    # First, open the header file and obtain the variable names:
    with open(headerFile) as omniHeaderFile:
        omniHeaderInfo = omniHeaderFile.readlines()
        # Ignore the initial lines in the file and just obtain the relevant info:
        omniVariables = omniHeaderInfo[4:]
        # Parse the variables in each line and collect them into a list of indices and a list of variable names:
        omniIndices = []
        omniVariableNames = []
        for element in omniVariables:
            segmentedElement = element.split()
            omniIndices.append(int(int(segmentedElement[0])-1))
            if segmentedElement[1] == 'Scalar' or segmentedElement[1] == 'Vector':
                omniVariableNames.append('Avg '+segmentedElement[1]+' B')
            else:
                omniVariableNames.append(segmentedElement[1].replace('_', ' ').replace('-', ' ').replace(',', ''))
        numVars = len(omniVariableNames)
        omniLabels = np.asarray(omniVariableNames[3:])

    # Second, open the data file and read in the data:
    with open(dataFile) as omniDataFile:
        omniData = omniDataFile.readlines()
        # Go over each line and collect the variables corresponding to each element in the line. For time values, collect them into datetime objects:
        omniDataArray = np.zeros((len(omniData), numVars-3))
        omniTimes = []
        for row in range(omniDataArray.shape[0]):
            omniLine = omniData[row].split()
            omniTimes.append(datetime(int(omniLine[0]), 1, 1) + timedelta(days=int(omniLine[1]) - 1) + timedelta(
                hours=int(omniLine[2])))
            for column in range(omniDataArray.shape[1]+3):
                if column > 2:
                    omniDataArray[row, column-3] = float(omniLine[column])
        omniTimes = np.asarray(omniTimes)
    return omniTimes, omniLabels, omniDataArray

def cleanF107(index_times, index_values, bad_value=999.9):
    """
    Given time stamps and observations of F10.7, determine the location of bad values (using the argument 'bad_value'),
    and replace them with imputed data generated with CSAPS.

    Parameters
    ----------
    index_times : arraylike
        An array or list of datetimes for each F10.7 value.
    index_values : arraylike
        An array or list of F10.7 values with the same shape as index_times.
    bad_value : float or int
        A single value corresponding to the bad values to be removed and imputed over. Default is 999.9.

    Returns
    -------
    clean_values : arraylike
        The gap-filled/imputed data.
    """
    # Check for bad values:
    badVals = np.where(index_values >= bad_value)[0]
    if len(badVals) >= 1:
        clean_values = index_values.copy()
        good_data_inds = np.logical_not(index_values >= bad_value)
        bad_data_inds = np.logical_not(good_data_inds)

        full_times_seconds = np.asarray([(x - index_times[0]).total_seconds() for x in index_times])

        good_data = clean_values[good_data_inds]
        good_times = full_times_seconds[good_data_inds]
        bad_times = full_times_seconds[bad_data_inds]

        yi = csaps.csaps(good_times, good_data, bad_times, smooth=0)

        subsetData = np.zeros_like(clean_values)
        subsetData[good_data_inds] = clean_values[good_data_inds]
        subsetData[bad_data_inds] = yi

        subsetInds = np.full(clean_values.shape, True)  # No subsetting since bad values were filled in

        # Since subsetting took place, change the actual values in the original data:
        clean_values[subsetInds] = subsetData

        # Plot the results for a sanity check:
        # import matplotlib.pyplot as plt
        # plt.figure()
        # plt.plot(index_times, index_values, label='Original')
        # plt.plot(index_times[subsetInds], clean_values[subsetInds], label='New')
        # plt.legend(loc='best')
        return clean_values
    else:
        return index_values

def F107filter(index_times, index_values, window_length=81, n=2):
    """
    Filter F10.7 values like so:
    1. Compute the running average and standard deviation.
    2. Identify the locations of values OUTSIDE the running average +/- n*stddev.
    3. For each datapoint identified in (2), replace it with the running average +/- n*stddev, depending on the sign.

    Parameters
    ----------
    index_times : arraylike
        An array or list of datetime values for each F10.7 value.
    index_values : arraylike
        An array or list of F10.7 values with the same shape as index_times.
    window_length : int
        The length of the running window over which to compute the running average and standard deviation. Default is 81.
    n : int
        The factor by which the standard deviation will be multiplied to determine which values should be replaced.
        Default is 2.

    Returns
    -------
     filteredF107 : arraylike
        The resulting filtered F10.7.
    """
    filteredF107 = index_values.copy()
    rollingMeanF107 = rollingAverage(index_values, window_length=window_length)
    rollingStdF107 = rollingStd(index_values, window_length=window_length)
    upperStdBoundary = np.add(rollingMeanF107, n*rollingStdF107)
    lowerStdBoundary = np.subtract(rollingMeanF107, n*rollingStdF107)
    # Identify locations where F10.7 EXCEED the stddev threshold:
    badLocsUpper = np.where(index_values > upperStdBoundary)[0]
    badTimesUpper = index_times[badLocsUpper]
    # Identify locations where F10.7 FALLS BELOW the stddev threshold:
    badLocsLower = np.where(index_values < lowerStdBoundary)[0]
    badTimesLower = index_times[badLocsLower]
    # Plotting the detected anomalous data:
    # import matplotlib.pyplot as plt
    # plt.figure()
    # plt.plot(index_times, index_values, label='Original')
    # plt.fill_between(index_times, upperStdBoundary, lowerStdBoundary, label='Standard Deviation Boundary', alpha=0.5)
    # plt.plot(index_times, rollingMeanF107, label='Rolling Mean')
    # # Plot the bad locs for EXCESS error:
    # for i in range(len(badTimesUpper)):
    #     plt.axvline(x=badTimesUpper[i], color='k', linestyle='-')
    # # Plot the bad locs for DEFECT error:
    # for i in range(len(badTimesLower)):
    #     plt.axvline(x=badTimesLower[i], color='k', linestyle='--')
    # plt.legend(loc='best')
    # IMPUTATION: Replace excesses with n*stddev
    filteredF107[badLocsUpper] = upperStdBoundary[badLocsUpper]
    # IMPUTATION: Replace defecsts with n*stddev
    filteredF107[badLocsLower] = lowerStdBoundary[badLocsLower]
    # Plot the resulting filtered data, along with the original data:
    # plt.figure()
    # plt.plot(index_times, index_values, label='Original')
    # plt.plot(index_times, filteredF107, label='Filtered')
    # plt.suptitle('Penticton F10.7 and Filtered Penticton F10.7 (81-day 2$\sigma$ criterion)')
    # plt.xlabel('Time')
    # plt.ylabel('F10.7 (sfu)')
    # plt.legend(loc='best')
    return filteredF107
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# # Execution
# if __name__=="__main__":
#     # PENTICTON DATA
#     saveLoc = '../solarIndices/F107/Penticton/'
#     fname = '../solarIndices/F107/Penticton/F107_1947_1996.txt'
#     fname1 = '../solarIndices/F107/Penticton/F107_1996_2007.txt'
#     fname2 = '../solarIndices/F107/Penticton/F107_current.txt'
#     times, f107 = readF107(fname)
#     times1, f1071 = readF107(fname1)
#     times2, f1072 = readF107(fname2)
#     # -------------------------------
#     # Combine all of the outputs from the files above and view the results as a sanity check:
#     allTimes = np.concatenate((times, times1, times2))
#     allF107 = np.concatenate((f107, f1071, f1072))
#     # Resample the times:
#     uniformTimes, uniformF107 = uniformSample(allTimes, allF107, cadence=24)
#     # Clean the data (through either gapification or imputation):
#     cleanedTimes, cleanedF107 = imputeData(uniformTimes, uniformF107, method='interp', bad_values=0)
#     # cleanedF107 = gapify(uniformF107, bad_value=0)
#     # Filter the data:
#     filteredF107 = F107filter(cleanedTimes, cleanedF107)
#     # Compute the centered rolling 81-day average of F10.7:
#     averagedF107 = rollingAverage(filteredF107, window_length=81)
#     # Compute the 54-day average of F10.7 with a backwards window:
#     backwardsAveragedF107 = rollingAverage(filteredF107, window_length=54, center=False)
#     # Plot as a sanity-check:
#     # import matplotlib; matplotlib.use('Qt5Agg'); import matplotlib.pyplot as plt; plt.figure(); plt.plot(cleanedTimes, cleanedF107); plt.plot(cleanedTimes, averagedF107); plt.plot(cleanedTimes, backwardsAveragedF107); plt.show()
#     # Save the data to pickle files:
#     savePickle(cleanedTimes, saveLoc+'F107times.pkl')
#     savePickle(cleanedF107, saveLoc+'F107vals.pkl')
#     savePickle(backwardsAveragedF107, saveLoc+'F107averageVals.pkl') # averagedF107
#     # -------------------------------
#     # Do all of the above with NASA OMNIWEB data:
#     omniSaveloc = '../solarIndices/F107/OMNIWeb/'
#     omniDataFile = '../solarIndices/F107/OMNIWeb/omni2_daily_r1GBiifQTW.lst'
#     omniHeaderFile = '../solarIndices/F107/OMNIWeb/omni2_daily_r1GBiifQTW.fmt'
#     omniTimes, omniLabels, omniData = readOMNI(omniDataFile, omniHeaderFile)
#     omniF107 = np.squeeze(omniData)
#     cleanOmniF107 = cleanF107(omniTimes, omniF107, bad_value=999.9)
#     filteredOMNIF107 = F107filter(omniTimes, cleanOmniF107)
#     averagedOmniData = rollingAverage(filteredOMNIF107, window_length=81)
#     # Save the data to pickle files:
#     savePickle(omniTimes, omniSaveloc + 'OMNIF107times.pkl')
#     savePickle(filteredOMNIF107, omniSaveloc + 'OMNIF107vals.pkl')
#     savePickle(averagedOmniData, omniSaveloc + 'OMNIF107averageVals.pkl')
#     # -------------------------------
#     sys.exit(0)