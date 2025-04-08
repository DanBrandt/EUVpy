# This module contains various helper functions that can be leveraged as generic tools.

#-----------------------------------------------------------------------------------------------------------------------
# Top-level Imports
import pandas as pd
import numpy as np
import math, pickle, os
from numpy import ma
from sklearn.impute import SimpleImputer
from datetime import datetime, timedelta
from scipy.interpolate import CubicSpline
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from sklearn.metrics import mean_squared_error
from tqdm import tqdm
import scipy.integrate as integ
from scipy import interpolate
from urllib.request import urlretrieve
from math import log10, floor
from scipy import stats
import seaborn as sns
# import scipy.optimize as opt
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Local imports
from EUVpy.empiricalModels.models.SOLOMON import solomon


#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Functions
def openDir(directory):
    """
    Create a directory, unless it already exists.
    :param directory:
        A string that is the directory to be created, with respect to the current location.
    """
    if os.path.isdir(directory) == False:
        os.makedirs(directory)
    else:
        print('Directory already exists: '+directory)

def urlObtain(url, fname):
    """
    Download a file from a URL. If the file already exists, don't download it; just print the location of the existing
    file.
    :param url: str
        The URL where the file is to be downloaded from.
    :param fname: str
        The filename which the downloaded file will be saved to.
    """
    if os.path.isfile(fname) == True:
        print('File already exists (loading in data): '+fname)
    else:
        urlretrieve(url, fname)

def savePickle(data, pickleFilename):
    """
    Given some data (a list, dict, or array), save it is a pickle file with a user-supplied name.
    :param: data
        A variable referring to data to be saved as a pickle.
    :param: pickleFilename, str
        A string with which to name the pickle file to be saved.
    """
    with open(pickleFilename, 'wb') as pickleFile:
        pickle.dump(data, pickleFile, protocol=pickle.HIGHEST_PROTOCOL)

def loadPickle(pickleFilename):
    """
    Given the name of a (pre-existing) pickle file, load its contents.
    :param: pickleFilename, str
        A string with the location/name of the filename.
    :return: var
        The loaded data.
    """
    with open(pickleFilename, 'rb') as pickleFile:
        var = pickle.load(pickleFile)
    return var

def firstNonNan(listfloats):
    """
    Find the index of the first non-NaN value in a given sequence.
    Source: https://stackoverflow.com/questions/22129495/return-first-non-nan-value-in-python-list
    :param listfloats:
    :return item:
        The element that is the first non-NaN value.
    :return idx:
        The index corresponding to the first non-NaN value.
    """
    indices = np.indices(np.asarray(listfloats).shape)[0]
    i = 0
    if np.isnan(listfloats[0]) == False:
        return listfloats[0], 0
    else:
        for item in listfloats:
            if math.isnan(item) == False:
                idx = indices[i]
                return item, idx
            i += 1

def find_nearest(array, value):
    """
    Given an array, find the index and value of an item closest to a supplied value.
    :param array: ndarray
        An array of values over which to search.
    :param value: float, int, str, datetime
        A value for which the closest value will be searched for.
    :return idx: int
        The index of the nearest value.
    :return array[idx]: float, int, str, datetime
        The actual value.
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]

def uniformSample(xdata, ydata, cadence):
    """
    Given some data with associated time stamps, resample the data so that the data have the time resolution equal to
    the given cadence.
    :param xdata: arraylike
        A list of time stamps.
    :param ydata: arraylike
        A list of data/values associated with the time stamps.
    :param cadence: int
        An integer denoting (in hours) the desired time cadence of the data.
    :return uniformXdata: ndarray
        A uniformly-sampled array of datetimes.
    :return uniformYdata: ndarray
        A uniformly-sampled array of corresponding data/values.
    """
    # Instantiate an array of timestamps with the desired cadence:
    start = datetime(xdata[0].year, xdata[0].month, xdata[0].day, 12)
    end = datetime(xdata[-1].year, xdata[-1].month, xdata[-1].day, 12)
    uniformXdata = np.arange(start, end, timedelta(hours=cadence)).astype(datetime)
    # Get indices of values nearest to the time array above:
    uniformIndices = []
    for i in range(len(uniformXdata)):
        goodIndex, goodValue = find_nearest(xdata, uniformXdata[i])
        uniformIndices.append( goodIndex )
    uniformIndicesArray = np.asarray(uniformIndices)
    # Extract the good values:
    uniformYdata = ydata[uniformIndicesArray]
    return uniformXdata, uniformYdata

def imputeData(timestamps, values, method='mean', bad_values=np.nan):
    """
    Given timeseries data, impute bad values and returned the cleaned data.
    :param timestamps: arraylike
        A 1D list or array of timestamps.
    :param values: arraylike
        A 1D list or array of timeseries data.
    :param method: str
        A string indicating the method to be used. If 'mean', 'median', 'most_frequent', or 'constant', uses the
        SimpleImputer routine from sklearn. If 'gam', uses a Generalized Additive Model to fill in the data. If 'gpr',
        use Gaussian Process Regression to perform imputation. Otherwise, cubic spline interpolation is used for
        imputation.
    :param bad_values: float, int, or NaN (type)
        The type of the data you wish to impute.
    :return cleanTimes: ndarray
        The timestamps corresponding to the cleaned data. This will differ from 'timestamps' if and only if some data
        are necessarily excluded by the chosen imputation method. For example, the edges of the data are often removed
        when the cubic spline interpolation is used for imputation, to avoid Runge's phenomenon.
    :return cleanData: ndarray
        The cleaned 1D timeseries.
    """
    cleanTimes = timestamps
    # Replace the bad values with NaNs:
    if bad_values != np.nan:
        bad_inds = np.where(values==bad_values)[0]
        values[bad_inds] = np.nan
    # Use the SimpleImputer routine from sklearn:
    if method=='mean' or method=='median' or method=='most_frequent' or method=='constant':
        imp = SimpleImputer(missing_values=np.nan, strategy=method)
        imp.fit(np.asarray(values).reshape(-1, 1))
        cleanData = imp.transform(values.reshape(-1, 1))
    elif method=='gpr':
        # TODO: Fix the Gaussian Process approach below so that sensible values are imputed:
        subset = values # [100:10000]
        DOYvals = np.array([fractionalDOY(element) for element in timestamps])
        DOYsubset = DOYvals # [100:10000]
        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import RBF
        kernel = 1.0 * RBF(length_scale=1, length_scale_bounds=(1e-3, 1e3))
        XPred = np.asarray(DOYsubset).reshape(-1, 1)
        XAxis = np.linspace(0,len(subset)-1,len(subset)).reshape(-1, 1)
        X = XAxis[~np.isnan(subset)]
        y = subset[~np.isnan(subset)].reshape(-1, 1)
        gpr = GaussianProcessRegressor(kernel=kernel, random_state=0).fit(X, y)
        yMeanPred, yStdPred = gpr.predict(XAxis, return_std=True)
        # plt.figure(); plt.plot(XAxis, yMeanPred, linestyle='-.'); plt.plot(XAxis, subset); plt.show()
    elif method=='gam':
        # Use a GAM to parameterize perform gap-filling:
        DOY = np.array([fractionalDOY(element) for element in timestamps])
        from pygam import LinearGAM, s
        X = DOY[~np.isnan(values)]
        y = values[~np.isnan(values)]
        gam = LinearGAM(s(0)).fit(X, y)
        cleanedData = gam.predict(DOY[np.isnan(values)])
        cleanData = values.copy()
        cleanData[np.isnan(values)] = cleanedData
    else:
        # Clip the ends of the data:
        firstNonNanValue, firstNonNanIndex = firstNonNan(values)
        lastNonNanValue, lastNonNanIndex = firstNonNan(np.flip(values, -1))
        lastNonNanIndex = len(values) - lastNonNanIndex
        subset = values[firstNonNanIndex:lastNonNanIndex]
        cleanTimes = timestamps[firstNonNanIndex:lastNonNanIndex]
        # Perform the imputation
        XAxis = np.linspace(0, len(subset) - 1, len(subset)).reshape(-1, 1)
        X = XAxis[~np.isnan(subset)]
        y = subset[~np.isnan(subset)].reshape(-1, 1)
        spl = InterpolatedUnivariateSpline(X, y, k=3)
        cleanData = spl(XAxis)
        # View the results as a sanity check:
        # plt.figure();
        # plt.plot(XAxis, cleanData, linestyle='-.');
        # plt.plot(XAxis, subset);
        # plt.show()

    # plt.figure(); plt.plot(timestamps, values); plt.plot(cleanTimes, cleanData); plt.show()

    return cleanTimes, np.squeeze(cleanData)

def gapify(timeseries, bad_value=999, replace_val=np.nan):
    """
    Take a timeseries and replace the bad values (signified by 'bad_values') with NaNs, unless otherwise desired.
    :param timestamps: arraylike
        A 1D list or array of timestamps.
    :param bad_values: int or float
        A float or int corresponding to the bad values to be 'gapped out'. Default is 999.
    :param replace_val: int, float, NaN, or None
        An int, float, NaN, or Nonetype describing what will replace the bad values.
    :return gappedData:
        The gapified 1D data.
    """
    gappedData = timeseries.copy()
    bad_inds = np.where(gappedData == bad_value)[0]
    gappedData[bad_inds] = replace_val
    return gappedData

def fractionalDOY(myDatetime):
    """
    Convert a generic datetime object into a float corresponding to the fractional day of the year.
    :param: myDatetime: datetime
        A datetime object.
    :return: fracDOY: float
        A float corresponding to the fractional Day of the Year, with the decimal portion included to show contributions
        from the time of day in hours, minutes, and seconds.
    """
    fracDOY = myDatetime.timetuple().tm_yday + myDatetime.hour/24. + myDatetime.minute/3600. + myDatetime.second/86400.
    return fracDOY

# TODO: Fix the function below to handle VERY LARGE data gaps.
def rollingAverage(myData, window_length=1, impute_edges=True, center=True):
    """
    Using pandas, compute a rolling average of over 'data' using a window length of 'windowlength'. Sets the leading and
    trailing windows to the values of the original data.
    :param myData: arraylike
        The data over which to compute the rolling average.
    :param window_length: int
        The size of the window over which to average.
    :param impute_edges: bool
        A boolean determining whether the edges will be interpolated. Default is True.
    :param center: bool
        A boolean determining whether the centered average will be used.
    :return: rolled, arraylike
        The rolling average data.
    """
    myDataframe = pd.DataFrame(data=myData, columns=['Var'])
    myDataframe['Rolling'] = myDataframe['Var'].rolling(window=window_length, center=center).mean()
    firstValidIndex = myDataframe['Rolling'].first_valid_index()
    lastValidIndex = myDataframe['Rolling'].last_valid_index()
    if impute_edges == True:
        # Sample x-axis:
        sampleXaxis = np.linspace(0, window_length, window_length)
        middleIndex = int(0.5*window_length)
        # Use cubic interpolation to fill the gaps on the edges:
        # leadingEdgeStartingVal = myDataframe['Var'][0]
        # goodLeadingVals = myDataframe['Var'][:window_length][myDataframe['Var'][:window_length] > 0]
        # leadingEdgeStartingVal = np.percentile(goodLeadingVals.values, 25)
        leadingEdgeStartingVal = myDataframe['Var'][:window_length].values[0]
        leadingEndingVal = myDataframe['Rolling'][firstValidIndex]
        leadingEdgeMiddleVal = np.mean([leadingEdgeStartingVal, leadingEndingVal])
        leadingSpline = CubicSpline([sampleXaxis[0], sampleXaxis[middleIndex], sampleXaxis[-1]],
                                    [leadingEdgeStartingVal, leadingEdgeMiddleVal, leadingEndingVal])
        leadingImputedValues = leadingSpline(sampleXaxis)
        # plt.figure(); plt.plot(sampleXaxis, myDataframe['Var'][:window_length].values); plt.plot(sampleXaxis, leadingImputedValues); plt.show()

        trailingEdgeStartingVal = myDataframe['Rolling'][lastValidIndex]
        trailingEndingVal = myDataframe['Var'].values[-1]
        trailingEdgeMiddleVal = np.mean([trailingEdgeStartingVal, trailingEndingVal])
        trailingSpline = CubicSpline([sampleXaxis[0], sampleXaxis[middleIndex], sampleXaxis[-1]],
                                    [trailingEdgeStartingVal, trailingEdgeMiddleVal, trailingEndingVal])
        trailingImputedValues = trailingSpline(sampleXaxis)
        # plt.figure(); plt.plot(sampleXaxis, myDataframe['Var'][-window_length:].values); plt.plot(sampleXaxis, trailingImputedValues); plt.show()

        myDataframe['Rolling'][:window_length] = leadingImputedValues
        myDataframe['Rolling'][-window_length:] = trailingImputedValues
    else:
        myDataframe['Rolling'][:window_length] = myDataframe['Var'][:window_length]
        myDataframe['Rolling'][-window_length:] = myDataframe['Var'][-window_length:]
    rolled = myDataframe['Rolling'].values
    return rolled

def rollingStd(myData, window_length=2, axis=-1):
    """
    Given some data, compute the rolling standard deviation. If the data is two dimensional, compute the rolling
    standard deviation along a specific axis of the data (specified by the user).
    :param myData: arraylike
        The data over which to compute the rolling average.
    :param window_length: int
        The size of the window over which to average. Default is 2.
    :param axis: int
        For 2D data, the axis along which to compute the rolling standard deviation. Defaults is -1.
    :return stdData: ndarray
        The rolling standard deviation values of the data.
    """
    # Define the generic rolling std function to be used repeatedly:
    def stdRoller(data, windowLength):
        myDataframe = pd.DataFrame(data=data, columns=['Var'])
        myDataframe['Rolling'] = myDataframe['Var'].rolling(window=windowLength, center=True).std()
        # Set the leading and trailing values equal to the mean standard deviation:
        meanStd = np.nanmean(myDataframe['Rolling'].values)
        infillData = np.full_like(myDataframe['Rolling'][:windowLength].values, fill_value=meanStd)
        myDataframe['Rolling'][:windowLength] = infillData
        myDataframe['Rolling'][-windowLength:] = infillData
        stdRes = myDataframe['Rolling'].values
        return stdRes
    # Perform the computation:
    if len(myData.shape) < 2:
        # Case for unidimensional data:
        stdData = stdRoller(myData, window_length)
    else:
        # Case for 2D data: Iterate through the data along the desired axis:
        stdData = np.zeros_like(myData)
        if axis == 0:
            # First Axis:
            for iRow in range(myData.shape[0]-1):
                stdData[iRow, :] = stdRoller(myData[iRow, :], window_length)
        elif axis == 1 or axis == -1:
            # Second Axis:
            for iCol in range(myData.shape[1]-1):
                stdData[:, iCol] = stdRoller(myData[:, iCol], window_length)
        else:
            raise ValueError('The axis specified exceeds the dimensions of the data.')
    return stdData

def normalize(myData, axis=-1):
    """
    Normalize data with respect to the mean, along the axis specified by the user.
    :param myData: ndarray
        A 1d or 2d array of data.
    :param axis: int
        The axis along which to perform normalization with respect to the mean. Default is -1.
    :return normedData: ndarray
        The normalized data.
    """
    def normFunc(data):
        return (data - np.nanmean(data)) / np.nanstd(data)
    if len(myData.shape) < 2:
        normedData = normFunc(myData)
    else:
        normedData = np.zeros_like(myData)
        if axis == 0:
            # First dimension:
            for iRow in range(myData.shape[0]-1):
                normedData[iRow, :] = normFunc(myData[iRow, :])
        elif axis==1 or axis==-1:
            # Second dimension:
            for iCol in range(myData.shape[1]-1):
                normedData[:, iCol] = normFunc(myData[:, iCol])
        else:
            raise ValueError('The axis specified exceeds the dimensions of the data.')
    return normedData

def corrCol(myData, otherData, saveLoc=None):
    """
    Given two sets of data, find the correlation between them. If the first set of data
    is multidimensional, correlate it with the second set by column. This function requires that
    'otherData' is 1D and equivalent to the lengths of the columns of 'myData'. Otherwise, if 'myData'
    is 1D, it must be the same length as 'otherData'.
    :param myData: ndarray
        The dependent variable data.
    :param otherData: ndarray
        The independent variable data.
    :param saveLoc: str
        A location where to save figures of the correlations. Optional argument. Default is None.
    :return fitParams: ndarray
        A poly1d object for the fit of the data. The last element is Pearson's R.
    """
    # Loop through the columns and perform the correlation:
    fitParams = []
    sortInds = np.argsort(otherData)
    sortedOtherData = otherData[sortInds]
    referenceData = np.linspace(np.nanmin(sortedOtherData), np.nanmax(sortedOtherData), num=100)
    for i in range(myData.shape[1]):
        currentData = myData[:, i][sortInds]
        # Find the best model fit, up to a polynomial of order 10:
        p, rss, order = bestPolyfit(sortedOtherData, currentData, 1)
        Rval = pearsonr(sortedOtherData, currentData)
        # View the data:
        plt.figure()
        plt.scatter(sortedOtherData, currentData, color='b')
        plt.plot(referenceData, p(referenceData), 'r-', label='Model Fit: Order '+str(order)+' (RSS= '+str(np.round(rss, 2))+')')
        text = "Pearson's R: "+str(np.round(Rval[0], 2))
        plt.plot([], [], ' ', label=text)
        plt.ylabel('$\sigma_{r_{\Theta_{\lambda}}}$')
        plt.xlabel('F10.7 (sfu)')
        plt.title('$\sigma_{r_{\Theta_{\lambda}}}$ vs. F10.7 (Band '+str(i+1)+')')
        plt.legend(loc='best')
        fitParams.append([p, Rval[0]])
        if saveLoc != None:
            plt.savefig(saveLoc+'neuvacResidStdCorrelation_Band'+str(i+1)+'.png', dpi=300)
            print('Plot saved to '+saveLoc+'neuvacResidStdCorrelation_Band'+str(i+1)+'.png')
    return fitParams

def bestPolyfit(xdata, ydata, maxOrder=5, func=None, **kwargs):
    """
    Given independent variable and dependent variable 1D data fit the data with polynomial functions of orders up to a
    user-defined limit set to 'maxOrder'. For the best-fitting model, return the parameters for the polynomial (as a
    poly1d object) along with the associated Residual Sum of Squares.
    :param xdata: ndarray
        1D independent variable data.
    :param ydata: ndarray
        1D dependent variable data.
    :param maxOrder: int
        A number below which (inclusive) to consider models to fit to the data. Default is 5.
    :param func: str
        Specifies what function should be used for fitting. Valid strings are: 'exp', 'linear', 'log', and 'cubic'. If
        this argument is passed, maxOrder is ignored.
    :return modelRes: list
        A list where the first element are the model parameters, the second is the Residual Sum of Squares, and the
        third is the model order.
    """
    # Example functions:
    if func is not None:
        if func == 'exp':
            def myFunc(x, a, b, c, d):
                return a*np.exp(-b*x + c) + d
        elif func == 'cubic':
            def myFunc(x, a, b, c, d):
                return a*x**3 + b*x**2 + c*x + d
        elif func == 'linear':
            def myFunc(x, a, b):
                return a*x + b
        elif func == 'log':
            def myFunc(x, a, b):
                a*np.log(x) + b
        elif func == 'quadratic':
            def myFunc(x, a, b, c):
                return a*x**2 + b*x + c
        else:
            raise ValueError('Invalid argument supplied for argument "func".')

    models = []
    rss_vals = []
    orders = []
    # Fit models of various orders:
    if func is not None:
        if func == 'exp':
            p0 = [2e-3, 2e-1, 3, 8e-7]
            newP0 = [2e-3, 0.1, np.nanmax(ydata), np.nanmin(ydata)]
        else:
            p0 = None
        popt, pcov = curve_fit(myFunc, xdata, ydata, p0=p0)
        # If convergence is not achieved, try again with new initializing parameters:
        if np.where(np.isinf(pcov))[0].shape[0] > 0:
            popt, pcov = curve_fit(myFunc, xdata, ydata, p0=newP0)
        p = popt
        # Plotting for a sanity check:
        sampleXdata = np.linspace(np.nanmin(xdata), np.nanmax(xdata), num=100)
        # plt.figure(); plt.scatter(xdata, ydata, color='b', label='Data'); plt.plot(sampleXdata, myFunc(sampleXdata, *popt), 'r-', label='Fit')
        # plt.plot(sampleXdata, myFunc(sampleXdata, *newP0), color='orange', linestyle='--'); plt.legend(loc='best')
        # Compute the R-squared (https://stackoverflow.com/questions/19189362/getting-the-r-squared-value-using-curve-fit):
        residuals = ydata - myFunc(xdata, *popt)
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((ydata - np.mean(ydata)) ** 2)
        r_squared = 1 - (ss_res / ss_tot)
        rss_vals.append(r_squared)
        models.append([myFunc, p])
        modelRes = [[myFunc, p], r_squared, func]
    else:
        for i in range(maxOrder):
            coeffs = np.polyfit(xdata, ydata, i+1)
            p = np.poly1d(coeffs)
            rss_vals.append(np.sum((ydata - p(xdata)) ** 2))
            models.append(p)
            orders.append(i+1)
        # Determine which model worked the best:
        locBest = np.nanargmin(np.asarray(rss_vals))
        modelRes = [models[locBest], rss_vals[locBest], orders[locBest]]
    return modelRes

def binRMSE(xdata, ydataEst, ydataTrue, step=10, saveLoc=None, titleStr=None, normalize=False):
    """
    Given some 1D independent variable data, some 1D estimates of dependent variable data, 1D true values of dependent
    variable data, and a step size, divide the xdata into bins of width equal to the step size and compute the RMSE
    error in each bin. Then compute the correlation between the RMSE and the binned xdata. Automatically saves
    a figure for the results at a user-defined location.
    :param xdata: ndarray
        1D independent variable data.
    :param ydataEst: ndarray
        1D estimates of dependent variable data.
    :param ydataEst: ndarray
        1D actual values of dependent variable data.
    :param step: int
        The bin width for the independent variable data. Default is 10.
    :param saveLoc: str
        A string for the location where the figure should be saved.
    :param titleStr: str
        A string for the title of the figure to be generated. Assumes that a single number representing a figure number
        or wavelength is given.
    :param normalize: Bool
        Controls whether the RMSE is normalized or not. Default is False.
    :return binCenters: list
        The bin centers for the dependent variable data.
    :return RMSE: list
        The RMSE values for each bin.
    """
    # Create the bins:
    start = round_mult(np.nanmin(xdata), step, direction='down')
    stop = round_mult(np.nanmax(xdata), step, direction='up')
    bins = np.arange(start, stop, step=step)
    binCenters = np.asarray([(a + b) / 2 for a, b in zip(bins[::2], bins[1::2])])
    # Loop the bins and compute the RMSE values:
    RMSE = []
    i = 0
    for element in binCenters:
        # Isolate the data that correspond to the given bin:
        goodInds = np.where((xdata >= bins[2*i-1]) & (xdata <= bins[2*i]))[0]
        if len(goodInds) > 0:
            goodEstData = ydataEst[goodInds]
            goodTrueData = ydataTrue[goodInds]
            rmse = mean_squared_error(goodTrueData[~np.isnan(goodTrueData)], goodEstData[~np.isnan(goodTrueData)], squared=False)
            if normalize == True:
                rmse = rmse / (np.nanmax(goodTrueData) - np.nanmin(goodTrueData))
        else:
            # If there is nothing in the bin, just record NaN for that bin:
            rmse = np.nan
        RMSE.append(rmse)
        i += 1
    # Plot the results for a sanity check:
    plt.figure()
    plt.plot(binCenters, RMSE, 'bo-')
    plt.xlabel('F10.7 (sfu)')
    plt.ylabel('RMSE (W/m$^2$/nm)')
    if titleStr != None:
        plt.suptitle('RMSE vs. F10.7: '+titleStr+' Angstroms')
    if saveLoc != None:
        plt.savefig(saveLoc+'RootMeanSquareDeviationVsF107_'+titleStr.replace('.', '_')+'.png', dpi=300)
    return binCenters, RMSE

def binCorrelation(xdata, ydataEst, ydataTrue, step=10, saveLoc=None, titleStr=None, root=False, normalize=False):
    """
    Given some 1D independent variable data, some 1D estimates of dependent variable data, 1D true values of dependent
    variable data, and a step size, divide the xdata into bins of width equal to the step size and compute the squared
    difference in each bin. Then compute the correlation between the the squared difference and the binned xdata.
    Automatically saves a figure for the results at a user-defined location.
    :param xdata: ndarray
        1D independent variable data.
    :param ydataEst: ndarray
        1D estimates of dependent variable data.
    :param ydataEst: ndarray
        1D actual values of dependent variable data.
    :param step: int
        The bin width for the independent variable data. Default is 10.
    :param saveLoc: str
        A string for the location where the figure should be saved.
    :param titleStr: str
        A string for the title of the figure to be generated. Assumes that a single number representing a figure number
        or wavelength is given.
    :param root: Bool
        Controls whether the square root of the SQDF is taken. Default is False.
    :param normalize: Bool
        Controls whether the SQDF is normalized (by dividing by the true data and multiplying by 100). Default is False.
        If True, ignores the value of 'root'.
    :return binCenters: list
        The bin centers for the dependent variable data.
    :return SQDF: list
        The Squared Difference values for each bin.
    """
    # Create the bins:
    start = round_mult(np.nanmin(xdata), step, direction='down')
    stop = round_mult(np.nanmax(xdata), step, direction='up')
    bins = np.arange(start, stop, step=step)
    binCenters = np.asarray([(a + b) / 2 for a, b in zip(bins[::2], bins[1::2])])
    # Loop the bins and compute the RMSE values:
    SQDF = []
    i = 0
    for element in binCenters:
        # Isolate the data that correspond to the given bin:
        goodInds = np.where((xdata >= bins[2*i-1]) & (xdata <= bins[2*i]))[0]
        if len(goodInds) > 0:
            goodEstData = ydataEst[goodInds]
            goodTrueData = ydataTrue[goodInds]
            N = len(goodTrueData[~np.isnan(goodTrueData)])
            sqdf = np.mean(
                np.square(np.subtract(goodEstData[~np.isnan(goodTrueData)], goodTrueData[~np.isnan(goodTrueData)])) / N)
            ylabelString = 'Squared Differences (W/m$^2$/nm)'
            realTitleString = 'Squared Differences'
            if normalize == False:
                if root == True:
                    sqdf = np.sqrt(sqdf)
                    ylabelString = 'RMSE (W/m$^2$/nm)'
                    realTitleString = 'RMSE'
            else:
                sqdf = np.mean(np.divide(
                    np.square(np.subtract(goodEstData[~np.isnan(goodTrueData)], goodTrueData[~np.isnan(goodTrueData)])),
                    goodTrueData[~np.isnan(goodTrueData)]) * 100)
                ylabelString = 'Normalized RMSE (W/m$^2$/nm)'
                realTitleString = 'Normalized RMSE'
        else:
            # If there is nothing in the bin, just record NaN for that bin:
            sqdf = np.nan
        SQDF.append(sqdf)
        i += 1
    # Plot the results for a sanity check:
    titleFontSize = 20
    fontSize = 18
    labelSize = 16
    plt.figure(figsize=(12, 8))
    plt.plot(binCenters, SQDF, 'mo-')
    plt.xlabel('F10.7 (sfu)', fontsize=fontSize)
    plt.ylabel(ylabelString, fontsize=fontSize)
    plt.tick_params(axis='both', labelsize=labelSize)
    if titleStr != None:
        plt.suptitle(realTitleString+' vs. F10.7: '+titleStr+' Angstroms', fontsize=titleFontSize)
    if saveLoc != None:
        plt.savefig(saveLoc+realTitleString.replace(' ', '_')+'vsF107_'+titleStr.replace('.', '_')+'.png', dpi=300)
    return binCenters, SQDF

def round_mult(num, divisor, direction='down'):
    """
    Round a number to the nearest integer multiple of a given divisor.
    :param num: int or float
        The number to round down.
    :param divisor: int or float
        The number which the result must be an integer multiple of.
    :param direction: str
        Either 'down' or 'up'. Specifies whether rounding should be done up or down/
    :return rounded: float
        The resulting number.
    """
    if direction == 'down':
        return num - (num%divisor)
    else:
        return divisor*(round(num/divisor))

def stringList(myList):
    """
    Convert a list of numbers into a single string.
    :param myList: list
        A list of numbers.
    :return myStr: str
        The string.
    """
    strList = [str(element) for element in myList]
    myStr = " ".join(strList)
    return myStr

def rebin(wavelengths, data, resolution, limits=None, factor=None, zero=True, unc=False):
    """
    Rebin 2D [irradiance/flux] data, where each row is a complete spectrum, and each column is a specific wavelength
    bin. Do the rebinning according to a desired resolution, and restrict the returned information to wavelength
    boundaries set by the user.
    :param wavelengths: arraylike
        A list of wavelengths at which each irradiance value is taken.
    :param data: ndarray
        A 2D array of spectrum data (irradiances or fluxes).
    :param resolution: float or dict
        The desired spectral resolution of the rebinned data. Units must be in nanometers. If argument is a dict, the
        key 'short' should contain wavelength boundaries for the beginning of bins, and the key 'long' should contain
        wavelength boundaries for the ending of bins. These boundaries should be in Angstroms.
    :param limits: list
        A list of two elements, where the first is a lower limit for wavelengths and the second is an upper limit for
        wavelengths, both in units of nanometers.
    :param factor: int
        A factor by which to upsample the data before performing rebinning. Tends to make the binning much more
        accurate, due to the use of integrals. Default is None.
    :param zero: bool
        Controls whether singular lines are set to a value of zero after they are extracted. Default is True.
    :param unc: bool
        Indicates whether or not uncertainties are being rebinned. If so, they will be combined in each bin using an
        interpolation method.
    :return newWaves: ndarray
        The centers of the wavelength bins corresponding to the rebinned data.
    :return newData: ndarray
        The rebinned data as a 2D array, arranged like the input data but at the new wavelength resolution.
    """
    # Get the native wavelength resolution of the input data:
    nativeResolution = np.concatenate((np.diff(wavelengths), np.array([np.diff(wavelengths)[-1]])), axis=0)
    nativeWavelengths = wavelengths.copy()

    # Upsample the wavelengths and irradiance data so that integrals are more accurate:
    if factor is not None:
        if len(data.shape) < 2:
            f = interpolate.interp1d(wavelengths, data, kind='linear')
            wavelengths = np.linspace(wavelengths[0], wavelengths[-1], factor * len(wavelengths))
            data = f(wavelengths)
        else:
            # Replace all NaNs with zeros:
            data[np.isnan(data)] = 0
            obs = np.linspace(0, data.shape[0], data.shape[0])
            f = interpolate.interp2d(wavelengths, obs, data, kind='linear')
            wavelengths = np.linspace(wavelengths[0], wavelengths[-1], factor*len(wavelengths))
            data = f(wavelengths, obs)

    # For uncertainty quantification, interpolate the grid of uncertainty values:
    if unc == True:
        # d_vals = data.flatten()
        # x_vals = np.repeat(wavelengths, int(d_vals.shape[0] / wavelengths.shape[0]))
        # numrows = data.shape[0]
        # y_vals = np.asarray([np.repeat(element, data.shape[1]) for element in np.linspace(0, numrows-1, numrows)]).flatten()
        # rbf_fun = Rbf(x_vals, y_vals, d_vals)
        newWaves = wavelengths * 10
        mids = 0.5 * (resolution['long'] + resolution['short'])
        if len(data.shape) < 2:
            f = interpolate.interp1d(newWaves, data, kind='linear')
            newData = f(mids)
        else:
            data[np.isnan(data)] = 0
            obs = np.linspace(0, data.shape[0], data.shape[0])
            f = interpolate.interp2d(newWaves, obs, data, kind='linear')
            newData = f(mids, obs)
        print()
    else:
        # Loop through each row the new data array, and fill it in with the rebinned data:
        if type(resolution) is not dict:
            # TODO: The code for this if cell is VERY SLOW when processing large arrays. It needs to be sped up.
            # CONDITION FOR A SINGLE WAVELENGTH RESOLUTION.
            # Compute the wavelength centers for the new bins:
            if limits is not None:
                newWaves = np.arange(limits[0], limits[-1] + resolution, resolution)
            else:
                step = np.nanmean(np.diff(wavelengths))
                newWaves = np.arange(wavelengths[0], wavelengths[-1] + step, step)
            # Instantiate the new data array:
            if len(data.shape) < 2:
                newData = np.zeros((1, newWaves.shape[0]))
            else:
                newData = np.zeros((data.shape[0], newWaves.shape[0]))
            # Loop over each spectrum (each individual observation):
            for i in tqdm(range(newData.shape[0])):
                # Loop through the wavelength intervals:
                currentSpectrum = np.zeros_like(newWaves)
                for j in range(newWaves.shape[0]-1):
                    # Get the indices of elements in each wavelength bin - ASSUME that the original data are associated with
                    # bin centers - apply an offset equal to the respective interval of the bin:
                    binInds = np.where((wavelengths > (newWaves[j] - resolution)) & (wavelengths < (newWaves[j] + resolution)))[0]
                    # Gather the elements in the bin and multiply them by the respective native bin width:
                    try:
                        binElements = data[i, binInds]
                    except:
                        binElements = data[binInds]

                    binWidths = nativeResolution[binInds]
                    binProducts = [a*b for a,b in zip(binElements, binWidths)]
                    # Take the sum:
                    binSum = np.sum(binProducts)

                    # Integrate the bin:
                    binVal = integ.trapz(binElements, wavelengths[binInds])
                    # binVal = binSum/(2*resolution)

                    # Add the value to the current spectrum:
                    currentSpectrum[j] = binVal

                # Once a single spectrum is obtained, add it to the new data array:
                newData[i, :] = currentSpectrum
        else:
            # CONDITION FOR A BINNING WITH UNIQUE BIN WIDTHS.
            shorts = resolution['short']/10.
            longs = resolution['long']/10.
            newWaves = 0.5*(shorts + longs)

            # Instantiate the new data array:
            if len(data.shape) < 2:
                newData = np.zeros((1, newWaves.shape[0]))
            else:
                newData = np.zeros((data.shape[0], newWaves.shape[0]))

            # First go through all of the wavelengths that are singular
            myData = data
            for iWave, short in enumerate(shorts):
                long = longs[iWave]
                if (long == short):
                    i = np.argmin(np.abs(wavelengths - short))
                    i2 = np.argmin(np.abs(nativeWavelengths - short))
                    try:
                        newData[:, iWave] = myData[:, i] * (nativeWavelengths[i2 + 1] - nativeWavelengths[i2])
                    except:
                        newData[:, iWave] = myData[i] * (nativeWavelengths[i2 + 1] - nativeWavelengths[i2])
                    if zero == True:
                        # Zero out bin so we don't double count it.
                        try:
                            myData[:, i] = np.zeros_like(myData[:, i])
                        except:
                            myData[i] = 0.0

            # Then go through the ranges
            for iWave, short in enumerate(shorts):
                long = longs[iWave]
                if (long != short):
                    d1 = np.abs(wavelengths - short)
                    iStart = np.argmin(d1)
                    d2 = np.abs(wavelengths - long)
                    iEnd = np.argmin(d2)
                    wave_int = 0.0
                    # For wavelengths at or below 0.2 nm, just compute the sum:
                    if long <= 0.2:
                        for i in range(iStart + 1, iEnd + 1):
                            newData[:, iWave] += myData[:, i] * \
                                                 (wavelengths[i + 1] - wavelengths[i])
                            wave_int += (wavelengths[i + 1] - wavelengths[i])
                    else:
                        try:
                            newData[:, iWave] = integ.trapz(myData[:, iStart:iEnd], wavelengths[iStart:iEnd], axis=1)
                        except:
                            newData[:, iWave] = integ.trapz(myData[iStart:iEnd], wavelengths[iStart:iEnd])

                        # # Plotting for a sanity check:
                        # plt.figure();
                        # plt.plot(wavelengths[iStart:iEnd], myData[iStart:iEnd], marker='o')
                        # plt.scatter(newWaves[iWave], newData[:, iWave])

                    # for i in range(iStart + 1, iEnd + 1):
                    #     binWidths = nativeResolution[iStart + 1:iEnd + 1]
                    #     try:
                    #         binProducts = [a * b for a, b in zip(myData[:, i], binWidths)]
                    #         binSum = np.sum(binProducts)
                    #         newData[:, iWave] = binSum / (iEnd - iStart)
                    #         # newData[:, iWave] += myData[:, i]  / (iEnd - iStart) # * (wavelengths[i + 1] - wavelengths[i])
                    #     except:
                    #         binProducts = [a * b for a, b in zip([myData[i]], binWidths)]
                    #         binSum = np.sum(binProducts)
                    #         newData[:, iWave] = binSum / (iEnd - iStart)
                    #         # newData[:, iWave] += (myData[i] * () )/ (iEnd - iStart)
                    #     wave_int += (wavelengths[i + 1] - wavelengths[i])

    # If only a single spectrum was generated, remove the singleton dimension:
    if newData.shape[0] == 1:
        newData = np.squeeze(newData)

    return newWaves, newData

def newbins(wavelengths, data, bins, zero=False, interpolation=False):
    """
    Rebin data according to a user-defined binning scheme.
    :param wavelengths: arraylike
        The wavelengths of the native data. Units should be in nm.
    :param data: arraylike
        The native data to be rebinned. Assumes each row is an observation and each column is a bin.
    :param bins: dict
        A dictionary containing the left-ward wavelength bin limits (first key) and the right-ward wavelength bin limits
        (right key). Units should be in Angstroms.
    :param zero: bool
        Controls whether a wavelength bin is 'zeroed' out after it is considered. Default is False.
    :param interpolation: bool
        Controls whether interpolation or 4-nm centered windowed averaging is used to handle singular emission
        lines. Default is False, corresponding to the windowed averaging.
    :return newWaves: arraylike
        The wavelength boundaries of the new binning scheme.
    :return newData: arraylike
        The values of the rebinned data.
    """
    # originalData = data
    nativeBinWidth = np.round(np.nanmean(np.diff(wavelengths)), 2)
    lowBins = bins['short']
    highBins = bins['long']
    newData = np.zeros((data.shape[0], len(lowBins)))
    newWaves = 0.5*(lowBins + highBins)
    # Relevant line indices in FISM2:
    fismLines = np.array([256, 284, 303, 304, 368, 465, 554, 584, 609, 629, 703, 765, 770, 787, 976, 1025, 1031, 1215])

    # Loop through the desired wavelengths:
    halfWindow = 2
    singulars = []
    j = 0
    for i in range(len(lowBins)): # 32
        # Go through the singular wavelengths first.
        if lowBins[i] == highBins[i]:
            # idx, val = find_nearest(wavelengths, lowBins[i] / 10.)

            # Simply assign the values at that index to the new bin values:
            # newData[:, i] = data[:, fismLines[j]] * nativeBinWidth

            if interpolation == False:
                # Sum the data in a 0.4-nm window centered on the line:
                newData[:, i] = np.sum(data[:, fismLines[j]-halfWindow:fismLines[j]+halfWindow] * nativeBinWidth, axis=-1 )
            else:
                # Aaronic method:
                # 1: Linearly-interpolate between the adjacent bins to approximate the background value
                interpolatedValues = np.zeros_like(data[:, 0])
                x = [wavelengths[fismLines[j]-1], wavelengths[fismLines[j]+1]]
                for k in range(data[:, fismLines[j]].shape[0]):
                    y = [data[k, fismLines[j]-1], data[k, fismLines[j]+1]]
                    interpVal = np.interp(wavelengths[fismLines[j]], x, y)
                    interpolatedValues[k] = interpVal

                # 2: Take the difference between the native FISM2 bin value and interpolated value
                difference = data[:, fismLines[j]] - interpolatedValues

                # 3: Assign the difference as the line irradiance (increase in brightness above background)
                newData[:, i] = difference * nativeBinWidth

                # 4: Sanity check - view the results:
                idx = k
                plt.figure()
                plt.plot(wavelengths, data[idx, :])
                plt.axvline(x=wavelengths[fismLines[j]], color='k')
                plt.scatter(x=wavelengths[fismLines[j]], y=data[idx, fismLines[j]], color='b')
                plt.scatter(x=wavelengths[fismLines[j]], y=interpolatedValues[idx], color='r')
                plt.scatter(x=wavelengths[fismLines[j]], y=difference[idx], color='m')

            singulars.append(lowBins[i])

            # Performing zeroing, if desired:
            if zero == True:
                # data[:, fismLines[j]] = np.zeros_like(data[:, fismLines[j]])
                if interpolation == False:
                    data[:, fismLines[j] - halfWindow:fismLines[j] + halfWindow] = np.zeros_like(data[:, fismLines[j]-halfWindow:fismLines[j]+halfWindow])
                else:
                    data[:, fismLines[j]-1:fismLines[j] + 1] = np.zeros_like(data[:, fismLines[j]-1:fismLines[j] + 1])
            j += 1

    # plt.figure(); plt.plot(wavelengths*10, originalData[0, :], 'bo-'); plt.plot(wavelengths*10, data[0, :], 'ro-')
    # for j in range(len(singulars)):
    #     plt.axvline(x=singulars[j], color='k')
    # plt.plot(newWaves, newData[0, :], 'go-')

    for i in range(len(lowBins)):
        # Go through the wavelength ranges:
        if lowBins[i] != highBins[i]:
            # binWidth = (highBins[i] - lowBins[i])/10.
            # Isolate the wavelength bands to sum together:
            validInds = np.where((wavelengths >= lowBins[i]/10.) & (wavelengths < highBins[i]/10.))[0]
            # Sum the values in the bins:
            newData[:, i] = np.sum(data[:, validInds] * nativeBinWidth, axis=-1)

    return newWaves, newData

def forecast(y, coefs, trend_coefs, steps, exog=None):
    """
    Produce linear minimum MSE forecast [taken from statsmodels.tsa.vector_ar.var_model.py]

    Parameters
    ----------
    y : ndarray (k_ar x neqs)
    coefs : ndarray (k_ar x neqs x neqs)
    trend_coefs : ndarray (1 x neqs) or (neqs)
    steps : int
    exog : ndarray (trend_coefs.shape[1] x neqs)

    Returns
    -------
    forecasts : ndarray (steps x neqs)

    Notes
    -----
    Lütkepohl p. 37
    """
    p = len(coefs)
    k = len(coefs[0])
    if y.shape[0] < p:
        raise ValueError(
            f"y must by have at least order ({p}) observations. "
            f"Got {y.shape[0]}."
        )
    # initial value
    forcs = np.zeros((steps, k))
    if exog is not None and trend_coefs is not None:
        forcs += np.dot(exog, trend_coefs)
    # to make existing code (with trend_coefs=intercept and without exog) work:
    elif exog is None and trend_coefs is not None:
        forcs += trend_coefs

    # h=0 forecast should be latest observation
    # forcs[0] = y[-1]

    # make indices easier to think about
    for h in range(1, steps + 1):
        # y_t(h) = intercept + sum_1^p A_i y_t_(h-i)
        f = forcs[h - 1]
        for i in range(1, p + 1):
            # slightly hackish
            if h - i <= 0:
                # e.g. when h=1, h-1 = 0, which is y[-1]
                prior_y = y[h - i - 1]
            else:
                # e.g. when h=2, h-1=1, which is forcs[0]
                prior_y = forcs[h - i - 1]

            # i=1 is coefs[0]
            f = f + np.dot(coefs[i - 1], prior_y)

        forcs[h - 1] = f

    return forcs

def forecastInversion(forecastData, history, lastDiff, trend, window):
    """
    Designed to operate on the output of the statsmodels VAR package. Takes a forecast generated by a VAR model, and
    for each sample in the forecast, inverts the forecast to obtain the actual prediction. Essentially, this function
    takes forecasts generated for non-stationary data that has been transformed to stationary data, and inverts it
    in order to give the non-stationary forecat result.
    :param forecastData: ndarray
        An n x m array of forecasts, where n is the number of forecasts and m is the variable.
    :param history: ndarray
        A 2D array of training data which was used to generate the forecastData.
    :param lastDiff, ndarray
        A 1 x m array of the most recent differences for all of the variables.
    :param trend: ndarray
        An 2D array of historical trend data for the non-stationary version of the data.
    :param window: int
        The window over which the historical trend was calculated (assumes some sort of centered average was done).
    :return invertedForecastDataNew: ndarray
        The inverted forecasted data.
    """
    data_and_preds_combined = history
    invertedForecastData = np.zeros_like(forecastData)
    # Loop over the number of samples:
    firstInversions = []
    for i in range(forecastData.shape[0]):
        # Undoing the differencing:
        if i == 0:
            firstInv = np.add(forecastData[i, :], lastDiff)
        elif i == 1:
            currentDiff = firstInversions[-1] - firstInv
            firstInv = np.add(forecastData[i, :], currentDiff)
        else:
            currentDiff = firstInversions[-1] - firstInversions[-2]
            firstInv = np.add(forecastData[i, :], currentDiff)
        firstInversions.append(firstInv)

        # Adding back on the trend
        if i == 0:
            secondInv = np.multiply(firstInv, trend[-1, :])
        else:
            currentTrend = np.mean(data_and_preds_combined[-window:], axis=0)
            secondInv = np.multiply(firstInv, currentTrend)
        # Enforce positive definiteness:
        secondInv[secondInv < 0] = 0
        invertedForecastData[i, :] = secondInv
        data_and_preds_combined = np.vstack([data_and_preds_combined, invertedForecastData])

    base = history[-1, :]
    offset = np.subtract(invertedForecastData[0, :], base)
    invertedForecastDataNew = invertedForecastData - offset

    # firstX = np.linspace(0, 26617, 26618)
    # secondX = np.linspace(26617, 26621, 5)
    # plt.figure()
    # plt.plot(firstX, history[:, 5])
    # plt.plot(secondX, invertedForecastDataNew[:, 5])

    return invertedForecastDataNew

def find_exp(number) -> int:
    base10 = log10(abs(number))
    return 1*(10**(floor(base10)))

def get_cc(array1, array2, normalize=True):
    """
    Compute the cross-correlation of two 1D arrays of the same length.
    :param array1: ndarray
        A 1D array of length n.
    :param array2: ndarray
        A 1D array of length n.
    :return c: float
        The normalized correlation of the two arrays.
    """
    if normalize:
        a = (array1 - np.mean(array1)) / (np.std(array1) * len(array1))
        b = (array2 - np.mean(array2)) / (np.std(array2))
        c = np.correlate(a, b)
    else:
        c = np.correlate(array1, array2)
    return c

def mycorrelate2d(df, normalized=False):
    """
    Compute the correlation matrix from 2D data, where each row is cross correlated with the others.
    This function handles NaN values by ignoring them.
    :param df: ndarray
        A 2D array of dimensions n x m.
    :param normalized: bool
        Determines whether the resulting correlation matrix is normalized. Default is False.
    :returns ccm: ndarray
        The [normalized] cross-correlation matrix.
    Source: https://stackoverflow.com/questions/54292947/basics-of-normalizing-cross-correlation-with-a-view-to-comparing-signals
    """
    # Initialize cross correlation matrix with zeros
    ccm = np.zeros((df.shape[1], df.shape[1]))
    # Fill in each entry of the matrix one-by-one:
    for i in range(df.shape[1]):
        outer_row = df[:, i]
        for j in range(df.shape[1]):
            inner_row = df[:, j]
            goodInds = np.logical_and(~np.isnan(outer_row), ~np.isnan(inner_row))
            if (not normalized):
                x = np.correlate(outer_row[goodInds], inner_row[goodInds])
            else:
                x = get_cc(outer_row[goodInds], inner_row[goodInds])
                # a = (inner_row - np.mean(inner_row)) / (np.std(inner_row) * len(inner_row))
                # b = (outer_row - np.mean(outer_row)) / (np.std(outer_row) )
                # x = np.correlate(a, b)
            ccm[i, j] = x
    return ccm

def linear(x, a, b):
    """
    A generic linear function for fitting lines with.
    """
    return a*x + b

def squareDiff(x, y):
    """
    Compute the difference squared between two arrays or lists (of the same length).
    """
    res = []
    for i in range(len(x)):
        res.append((y[i] - x[i])**2)
    return np.asarray(res)

def mape(x, y):
    """
    Compute the mean absolute percentage error between two arrays.
    :param x: arraylike
        The modeled values.
    :param y: arraylike
        The actual (true) values.
    :return mape:
        The mean absolute percentage error.
    """
    n = len(x)
    quantity = np.abs(np.divide(np.subtract(np.asarray(y), np.asarray(x)), np.asarray(y)))
    mape = (1./n) * np.sum(quantity)
    return mape

def plotHist(data, bins, color, saveLoc=None, labels=None, logScale=None, density=True):
    """
    Given data and bins, plot a histogram and fit a distribution to it.
    :param data: arraylike
        1D data to make a histogram with.
    :param bins: arraylike
        The bins with which to bin the data.
    :param color: str
        The color of the histogram data that will be in the plot.
    :param saveLoc: str
        A string containing the location with which to save the file.
    :param labels: list
        A 3-element list containing string for the xlabel, ylabel, and the title.
    :param logScale: bool
        If 'x', scales the x-axis on a log scale. If 'y', scales the y-axis on a log scale. If 'both', scales both.
    :param density: bool
        If True, computes the density curve of the histogram. If False, simply fits a skew normal distribution to the
        data.
    """
    validLocs = ~np.isnan(data)
    xVals = np.linspace(bins[0], bins[-1], 500)
    cleanData = data[validLocs]
    prunedData = cleanData[np.where((cleanData >= bins[0]) & (cleanData <= bins[-1]))[0]]
    if density == False:
        fig, ax = plt.subplots()
        ax.hist(data, bins=bins, density=True, color=color, label='Data')
        a, loc, scale = stats.skewnorm.fit(prunedData)
        p = stats.skewnorm.pdf(xVals, a, loc, scale)
        ax.plot(xVals, p, 'k', linewidth=2, label=r'Fit: $\alpha$='+str(np.round(a, 2))+r', $\xi$='+str(np.round(loc, 2))+r', $\omega$='+str(np.round(scale, 2)))
    else:
        fig, ax = plt.subplots()
        sns.histplot(prunedData, kde=True, bins=bins, color=color, line_kws={'lw': 5}).lines[0].set_color('black')
    # Text and labeling:
    wavelength_str = labels[2][-18:] # Valid only for wavelengths under 1000 Angstroms
    textstr = r''+wavelength_str+'\n'+'$\mu$: '+str(np.round(np.mean(prunedData), 2))+'\n $\sigma$: '+str(np.round(np.std(prunedData), 2))
    plt.text(x=0.9, y=0.8, s=textstr, ha='center', va='center', fontsize=50, transform = ax.transAxes)
    plt.xlabel(labels[0], fontsize=37)
    plt.ylabel(labels[1], fontsize=37)
    # plt.title(labels[2], fontsize=35)
    # plt.legend(loc='best')
    plt.tick_params(axis='both', labelsize=35)
    plt.xlim([-25., 25.])
    if logScale == 'x':
        plt.xscale('log')
    elif logScale == 'y':
        plt.yscale('log')
    elif logScale == 'both':
        plt.xscale('log')
        plt.yscale('log')
    else:
        pass
    plt.axvline(x=0, color='gray', linestyle='--', alpha=0.8, linewidth=5)
    plt.tight_layout()
    if density == False:
        print('Skewnormal Parameters: Alpha=' + str(a) + ', Loc=' + str(loc) + ', Scale=' + str(scale)+', Kurtosis='+
              str(stats.kurtosis(data[validLocs]))+', Skew='+str(stats.skew(data[validLocs])))
    if saveLoc:
        plt.savefig(saveLoc, dpi=300)
        print('Saved figure to '+saveLoc)
    return fig

def percDev(x, y):
    """
    Compute the percentage deviation (a.k.a. percentage error).
    :param x: float or arraylike
        The observed value(s).
    :param y: float or arraylike
        The true value(s).
    """
    return np.divide(np.subtract(x, y), y) * 100

def readFISM2(fname):
    """
    Given a FISM2 netcdf file, read it in and output a .dat file with the same format as the .dat files here:
    https://github.com/aaronjridley/GITM/tree/master/srcData/FISM
    """
    from netCDF4 import Dataset
    fism2Data = Dataset(fname)

    irradiance = np.asarray(fism2Data.variables['irradiance'])  # W/m^2/nm
    dates = fism2Data.variables['date']
    datetimes = []
    for i in range(len(dates)):
        year = dates[i][:4]
        day = dates[i][4:]
        currentDatetime = datetime(int(year), 1, 1) + timedelta(int(day) - 1) + timedelta(hours=12)
        datetimes.append(currentDatetime)

    def numStr(num):
        if int(num) < 10:
            return ' '+str(int(num))
        else:
            return str(int(num))

    i = 0
    with open('fism2irr_daily.dat', 'w') as fism2File:
        fism2File.write('#START\n')
        dateStr = str(datetimes[0])[:4]+' '+numStr(str(datetimes[0])[5:7])+' '+numStr(str(datetimes[0])[8:10])
        if i == 0:
            timeStr = '  0  0  0 '
        else:
            timeStr = ' 12  0  0 '
        dateTimeStr = dateStr + timeStr
        strList = [str(element) for element in irradiance[i, :]]
        irrStr = ' '.join(strList)
        fileStr = dateTimeStr + irrStr + '\n'
        fism2File.write(fileStr)
        i += 1
    print('Wrote file to: '+os.getcwd()+'/fism2irr_daily.dat')

def readCLS(filename):
    """
    Load in flare-corrected, Sun-Earth distance adjusted flux values recorded by the Collecte Localisation Satellites
    (CLS).
    :param filename: str
        The location of the data file.
    :return times: list
        The datetimes for each data value.
    :return data: ndarray
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
                data[j, :] = np.array([float(elements[5]), float(elements[9]), float(elements[13]), float(elements[17]), float(elements[21])])
                times.append( datetime(int(elements[0]), int(elements[1]), int(elements[2]), 12) )
                precisionVals.append( [float(elements[6]), float(elements[10]), float(elements[14]), float(elements[18]), float(elements[22])] )
                j += 1
            i += 1
    # Print the precision:
    print('Mean precision values...')
    print('F30: '+str(np.nanmean([element[0] for element in precisionVals]))+' sfu') # 6
    print('F15: ' + str(np.nanmean([element[1] for element in precisionVals])) + ' sfu') # 8
    print('F10.7: ' + str(np.nanmean([element[2] for element in precisionVals])) + ' sfu') # 13
    print('F8: ' + str(np.nanmean([element[3] for element in precisionVals])) + ' sfu') # 12
    print('F3.2: ' + str(np.nanmean([element[4] for element in precisionVals])) + ' sfu') # 11
    return times, data

def band_info(bin_boundaries_left, bin_boundaries_right, solomon=False):
    """
    Given a two array-likes of wavelength bin boundaries, return the indices of those corresponding only to wavelength
    ranges, rather than individual lines, and also return the wavelength boundaries themselves.
    :param bin_boundaries_left: array-like
        The left-sided wavelength band boundaries. May either be in nm or Angstroms.
    :param bin_boundaries_right: array-like
        The right-sided wavelength band boundaries. May either be in nm or Angstroms. Units should match those of
        bin_boundaries_left.
    :param solomon: boolean
        If True, does the same thing, but accounts for overlapping bands in the solomon bins. Default is False.
    :return sortInds: array-like
        The indices of the wavelength ranges only.
    :return xPosSorted: array-like
        The values of the wavelength ranges themselves, for the purposes of making spectral stair plots.
    """
    if not solomon:
        # Individual wavelength lines should have the same value for the left and right boundaries, so find the indices where
        # that is NOT true:
        bandInds = np.where(np.asarray(bin_boundaries_left) != np.asarray(bin_boundaries_right))[0]
        # Append the right-most boundary to ensure that the resulting boundary values are 'closed' on both sides.
        xPos = np.append(np.asarray(bin_boundaries_left)[bandInds],
                         np.asarray(bin_boundaries_right)[-1])
        # Sort the indices:
        sortInds = np.argsort(xPos)
        # Return the sorted boundary values for ease of plotting:
        xPosSorted = xPos[sortInds]
        return sortInds[:-1], xPosSorted
    else:
        xPosSolomonInitial = 0.5 * (bin_boundaries_left + bin_boundaries_right)
        xPosSolomon = np.append(xPosSolomonInitial, 1130.)
        bandWidths = np.append(bin_boundaries_right, 1130.) - np.append(bin_boundaries_left, bin_boundaries_right[-1])
        # Constrain the bins; i.e. for the overlapping bins, we will want to sum everything in there together:
        xPosSolomonNew = []
        j = 0
        for i in range(len(xPosSolomon) - 1):
            if xPosSolomon[i] != xPosSolomon[i + 1]:
                xPosSolomonNew.append(xPosSolomon[i])
                j += 1
        xPosSolomon = np.asarray(xPosSolomonNew)
        sortIndsSolomon = np.argsort(xPosSolomon)
        # The edges of the bins are what we want:
        xPosSortedSolomon = np.array([0.5, 4, 8, 18, 32, 70, 155, 224, 290, 320, 540, 650, 798, 913, 975, 987,
                                      1027])
        return sortIndsSolomon[:-1], xPosSortedSolomon

def solomonRebin(solomonIrr):
    """
    Given some irradiance in the Solomon Bins, do a simple rebinning procedure to simply sum the irradiance in the
    overlapping bins, so that the resulting spectra can be plotted rather easily.
    :param solomonIrr: numpy array
        An array of irradiance in the Solomon Bins. Rows are observations and columns are wavelength bins.
    :return solomonIrr_for_plotting: numpy array
        The rebinned irradiance.
    """
    solomonTable = solomon.solomonBands
    # sortIndsSolomon, xPosSortedSolomon = band_info(solomonTable['short'], solomonTable['long'], solomon=True)
    xPosSolomonInitial = 0.5 * (solomonTable['short'] + solomonTable['long'])
    xPosSolomon = np.append(xPosSolomonInitial, 1130.)
    goodInds = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 16, 20, 21])  # badInds: 12, 14, 15, 17, 18
    solomonIrr_for_plotting = solomonIrr[:, goodInds]
    j = 0
    for i in range(len(xPosSolomon)-1):
        if xPosSolomon[i] != xPosSolomon[i+1]:
            # Do nothing
            j += 1
            pass
        else:
            solomonIrr_for_plotting[:, j] += solomonIrr[:, i]
    return solomonIrr_for_plotting
########################################################################################################################
# TODO: Python versions of Liying Qian's Wavelength Rebinning Code from here: https://download.hao.ucar.edu/pub/lqian/tlsm/idl/



#-----------------------------------------------------------------------------------------------------------------------