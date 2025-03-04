#!/usr/bin/env python

import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as dates
import argparse
import os
import re

# ----------------------------------------------------------------------
# Function to parse input arguments
# ----------------------------------------------------------------------

def parse_args_fism():

    parser = argparse.ArgumentParser(description = 'Download and process FISM data')
    parser.add_argument('-start', \
                        help = 'start date as YYYYMMDD', default = '0')
    parser.add_argument('-end', \
                        help = 'end date as YYYYMMDD', default = '0')
    parser.add_argument('-euvfile', \
                        help='EUV file that provides wavelength bins',
                        default = 'euv_59.csv')
    parser.add_argument('-fismfile', \
                        help='FISM file to convert',
                        default = 'none')
    parser.add_argument('-gitm', \
                        help='Write out GITM-style file',
                        action='store_true')
    parser.add_argument('-flare', \
                        help='Download and process FISM flare data',
                        action='store_true')

    args = parser.parse_args()

    return args

# ----------------------------------------------------------------------
# Take string time YYYYMMDD.HHMM and convert to datetime
# ----------------------------------------------------------------------

def convert_ymdhm_to_dt(ymdhm):

    year = ymdhm[0:4]
    mo = ymdhm[4:6]
    da = ymdhm[6:8]

    if (len(ymdhm) >= 11):
        hr = ymdhm[9:11]
        if (len(ymdhm) >= 13):
            mi = ymdhm[11:13]
        else:
            mi = '00'
    else:
        hr = '00'
        mi = '00'
        
    time = dt.datetime(int(year), int(mo), int(da), int(hr), int(mi), 0)
    return time

# ----------------------------------------------------------------------
# Function to download FISM2 data.
# Need call to be in form:
# https://lasp.colorado.edu/lisird/latis/dap/fism_daily_hr.csv?&time>=2020-01-01T00:00:00.000Z&time<=2020-01-31T00:00:00.000Z
# ----------------------------------------------------------------------

def download_fism2(start, end, isFlare):

    sStart = start.isoformat()+'.000Z'
    sEnd = end.isoformat()+'.000Z'

    if (isFlare):
        site = 'https://lasp.colorado.edu/lisird/latis/dap/fism_flare_hr.csv'
    else:
        site = 'https://lasp.colorado.edu/lisird/latis/dap/fism_daily_hr.csv'

    url = site + '?&time>=' + sStart + '?&time<=' + sEnd
    
    ymdS = start.strftime('%Y%m%d')
    ymdE = end.strftime('%Y%m%d')
    filename = '.fism_raw_' + ymdS + '_to_' + ymdE + '.txt'

    command = 'curl "' + url + '" > ' + filename
    print("Running Command : ", command)
    os.system(command)
    
    return filename
    

#------------------------------------------------------------------------------
# Convert yearday (YYYYDDD) to date time
#   - center the day at 12 UT.
#------------------------------------------------------------------------------

def convert_time(yearday):
    year = int(yearday/1000.0)
    base = dt.datetime(year,1,1,12)
    doy = yearday - year * 1000.0
    time = base + dt.timedelta(days = doy-1)
    return time

#------------------------------------------------------------------------------
# Convert seconds since 1970 to date time
#   - center the day at 12 UT.
#------------------------------------------------------------------------------

def convert_time_seconds(seconds):
    base = dt.datetime(1970,1,1)
    time = base + dt.timedelta(seconds = seconds)
    return time

#------------------------------------------------------------------------------
# read euv file - this determines the wavelength bins
#------------------------------------------------------------------------------

def read_euv_csv_file(file, band=False):

    fpin = open(file, 'r')

    iFound = 0
    for line in fpin:
        aline = line.split(',')
        s = aline[-1].strip().split('.')[0]
        if (aline[0].strip() == "Short"):
            if (s.isnumeric()):
                short = np.asarray(aline[5:], dtype = float)
            else:
                short = np.asarray(aline[5:-1], dtype = float)
            iFound += 1
        if (aline[0].strip() == "Long"):
            if (s.isnumeric()):
                long = np.asarray(aline[5:], dtype = float)
            else:
                long = np.asarray(aline[5:-1], dtype = float)
            iFound += 1
        if (iFound == 2):
            break
    wavelengths = {'short': short,
                   'long': long}
    if band==True:
        diffs = long - short
        bandInds = np.where( diffs != 0)[0]
        wavelengths = {'short': short[bandInds],
                       'long': long[bandInds]}
    return wavelengths
    
#------------------------------------------------------------------------------
# read fism2 csv file
#------------------------------------------------------------------------------

def read_fism_csv_file(file):

    fpin = open(file, 'r')

    header = fpin.readline()
    vars = header.split(',')

    isSeconds = False

    m = re.match(r'.*seconds.*',vars[0])
    if m:
        isSeconds = True

    # Read in a 1d time array and 1d wavelength array
    # read in all of the irradiance and uncertainty data as
    # a large 1d array, then reshape it into a time vs wavelength 2d array
    # structure of the file is all wavelengths for one time, then
    # a new time and all wavelengths for that time.
    
    iRow = 0
    oldtime = 0.0
    nWaves = 0
    nTimes = 0
    times = []
    wavelengths = []
    irradiance = []
    uncertainty = []
    for line in fpin:
        aline = line.split(',')
        yearday = float(aline[0])
        irradiance.append(float(aline[2]))
        uncertainty.append(float(aline[3]))
        if (iRow == 0):
            oldtime = yearday

        if (yearday > oldtime):
            nWaves = 1
            if isSeconds:
                times.append(convert_time_seconds(oldtime))
            else:
                times.append(convert_time(oldtime))
            oldtime = yearday
            nTimes = nTimes + 1
        else:
            nWaves = nWaves + 1

        if (nTimes == 0):
            wavelengths.append(float(aline[1]))
            
        iRow = iRow+1

    if isSeconds:
        times.append(convert_time_seconds(oldtime))
    else:
        times.append(convert_time(oldtime))
    nTimes = nTimes + 1
    
    irradiance = np.array(irradiance).reshape((nTimes, nWaves))
    uncertainty = np.array(uncertainty).reshape((nTimes, nWaves))
    
    fism_data = {'vars': vars,
                 'time': times,
                 'nTimes': nTimes,
                 'wave': np.array(wavelengths)*10.0,  # convert to Angstroms
                 'irr': irradiance,
                 'unc': uncertainty}
            
    return fism_data

#------------------------------------------------------------------------------
# rebin FISM data into new wavelength bins
#  some bins are single wavelength, and some span lots of wavelenghts
#------------------------------------------------------------------------------

def rebin_fism(fism_waves, fism_vals, wavelengths):

    shorts = wavelengths['short']
    longs = wavelengths['long']
    nWaves = len(shorts)
    nvals = fism_vals.shape[0]
    new_irr = np.zeros((nvals, nWaves)) # np.zeros(nWaves)
    ave_wav = np.zeros(nWaves) # np.zeros(nWaves)

    # first go through all of the wavelengths that are singular
    for iWave, short in enumerate(shorts):
        long = longs[iWave]
        if (long == short):
            d = np.abs(fism_waves - short)
            i = np.argmin(d)
            new_irr[:, iWave] = fism_vals[:, i] * \
                    ((fism_waves[i+1] - fism_waves[i])/10.0)
            # zero out bin so we don't double count it.
            fism_vals[i] = 0.0

    # then go through the ranges
    for iWave, short in enumerate(shorts):
        long = longs[iWave]
        ave_wav[iWave] = (short + long)/2.0
        if (long != short):
            d = np.abs(fism_waves - short)
            iStart = np.argmin(d)
            d = np.abs(fism_waves - long)
            iEnd = np.argmin(d)
            wave_int = 0.0
            for i in range(iStart+1, iEnd+1):
                new_irr[:, iWave] += fism_vals[:, i] * \
                    ((fism_waves[i+1] - fism_waves[i])/10.0)
                wave_int += (fism_waves[i+1] - fism_waves[i])/10.0
    return new_irr, ave_wav

#------------------------------------------------------------------------------
# isolate FISM data into new wavelength bins
# choose those closest to the center wavelength for all bins
#------------------------------------------------------------------------------
# from scipy.interpolate import UnivariateSpline
def isolate_fism(fism_waves, fism_vals, wavelengths):

    shorts = wavelengths['short']
    longs = wavelengths['long']
    mids = (shorts + longs)*0.5
    nWaves = len(shorts)
    new_irr = np.zeros((fism_vals.shape[0], nWaves)) # np.zeros(nWaves)
    ave_wav = np.zeros(nWaves)

    # first go through all of the wavelengths that are singular
    for iWave, short in enumerate(shorts):
        long = longs[iWave]
        if (long == short):
            d = np.abs(fism_waves - short)
            i = np.argmin(d)
            new_irr[:, iWave] = fism_vals[:, i] * ((fism_waves[i + 1] - fism_waves[i]) / 10.0)
            # zero out bin so we don't double count it.
            fism_vals[i] = 0.0

    # then go through the ranges
    for iWave, short in enumerate(shorts):
        mid = mids[iWave]
        long = longs[iWave]
        ave_wav[iWave] = (short + long)/2.0
        if (long != short):
            d1 = np.abs(fism_waves - short)
            iStart = np.argmin(d1)
            d2 = np.abs(fism_waves - long)
            iEnd = np.argmin(d2)
            wave_int = 0.0
            for i in range(iStart + 1, iEnd + 1):
                new_irr[:, iWave] += fism_vals[:, i] * \
                                  ((fism_waves[i + 1] - fism_waves[i]) / 10.0)
                wave_int += (fism_waves[i + 1] - fism_waves[i]) / 10.0

            # Add up all of the irradiances in a bin, but average the irradiance value over that bin:
            #for i in range(iStart + 1, iEnd + 1):
            # irrMaxes[iWave] = max(fism_vals[iStart+1:iEnd+1])
            # diffs = (fism_waves[iStart+1:iEnd+1]-fism_waves[iStart:iEnd])/10.
            # try:
            #     # new_irr[iWave] = np.sum(fism_vals[iStart+1:iEnd+1])/(iEnd-iStart) # Standard Average (wrt number of elements)
            #     # new_irr[iWave] = np.sum(fism_vals[iStart+1:iEnd+1]*fism_waves[iStart+1:iEnd+1*diffs]) / ((fism_waves[iEnd+1] + fism_waves[iStart+1])/2.) # Wavelength-specific average
            #     # expVal = 2
            #     myWeights = fism_vals[iStart+1:iEnd+1] / sum(fism_vals[iStart+1:iEnd+1]) # fism_vals[iStart+1:iEnd+1] / np.mean(fism_vals[iStart+1:iEnd+1]) ** expVal # Weighted sum
            #     # new_irr[iWave] = np.average(fism_vals[iStart+1:iEnd+1], weights=myWeights)
            #     # Form a univariate spline and take its central value:
            #     # print(len(fism_waves[iStart+1:iEnd+1]))
            #     # spl = UnivariateSpline(fism_waves[iStart+1:iEnd+1], fism_vals[iStart+1:iEnd+1], k=3)
            #     # new_irr[iWave] = np.average(spl(fism_waves[iStart+1:iEnd+1]), weights=fism_vals[iStart+1:iEnd+1] / (np.sum(fism_vals[iStart+1:iEnd+1]))) #fism_vals[iStart+1:iEnd+1])
            #     # Integrals:
            #     # new_irr[iWave] = np.sum(fism_vals[iStart+1:iEnd+1] * (fism_waves[iEnd+1] - fism_waves[iStart+1])/10.) / len(fism_vals[iStart+1:iEnd+1])
            # except:
            #     new_irr[iWave] = np.mean(fism_vals[iStart+1:iEnd+1])
        # d = np.abs(fism_waves - mid)
        # i = np.argmin(d)
        # choices[iWave] = fism_waves[i]
        # new_irr[iWave] = fism_vals[i]

    # plt.figure(figsize=(16, 10))
    # # plt.plot(shorts, 'ro-')
    # # plt.plot(longs, 'bo-')
    # # plt.plot(np.linspace(0, len(longs), len(fism_waves)), fism_waves, 'go-')
    # # plt.plot(choices, 'mo-')
    # # plt.plot(ave_wav, 'co-')
    # # plt.semilogy(fism_waves, fism_vals, 'ko', label='FISM2')
    # plt.plot(ave_wav, new_irr, 'yo', label='FISM2Ave')
    # # euvacSingulars = np.array([  75.,  125.,  175.,  225.,  275.,  325.,  375.,  425.,  475.,
    # #     525.,  575.,  625.,  675.,  725.,  775.,  825.,  875.,  925.,
    # #     975., 1025.])
    # euvacAll = np.array([75.,  125.,  175.,  225.,  256.3,  284.15,  275.,
    #     303.31,  303.78,  325.,  368.07,  375.,  425.,  465.22,
    #     475.,  525.,  554.37,  584.33,  575.,  609.76,  629.73,
    #     625.,  675.,  703.31,  725.,  765.15,  770.41,  789.36,
    #     775.,  825.,  875.,  925.,  977.02,  975., 1025.72,
    #    1031.91, 1025.])
    # euvacIrr = np.array([5.69353262e-05, 1.32411915e-05, 9.37927458e-05, 4.36505616e-05,
    #    3.46218785e-04, 0.00000000e+00, 1.75707120e-05, 3.87532953e-04,
    #    4.35535292e-03, 9.03853673e-06, 3.26719434e-04, 2.06070896e-06,
    #    3.30406238e-06, 1.14166610e-04, 1.88167308e-06, 3.10849797e-06,
    #    2.49217276e-04, 8.68566411e-05, 2.37267582e-06, 1.43577129e-04,
    #    4.84189165e-04, 1.90682223e-06, 1.29357633e-06, 9.84585057e-05,
    #    7.34438803e-07, 4.23650873e-05, 5.80973952e-05, 1.70635152e-04,
    #    3.69285746e-06, 7.43286848e-06, 1.51104669e-05, 1.22163808e-05,
    #    8.57872454e-04, 5.73383084e-06, 6.41234926e-04, 3.82012071e-04,
    #    9.12641890e-06])
    # # np.array([7.65540784e-05, 1.63809850e-05, 1.38699401e-04, 7.64534817e-05,
    # #    3.73511747e-05, 1.71982336e-05, 5.80616255e-06, 4.12128830e-06,
    # #    3.36710758e-06, 4.03152129e-06, 2.65070695e-06, 2.69721293e-06,
    # #    1.47153941e-06, 8.47519491e-07, 4.26354587e-06, 8.59421745e-06,
    # #    1.79185283e-05, 1.41947633e-05, 6.55172863e-06, 1.04153719e-05])
    # ridleyIrr = np.array([2.01523195e-05, 8.00565189e-06, 1.61690278e-04, 4.27308566e-05,
    #    7.49159642e-04, 1.97787269e-03, 6.21483716e-05, 3.08510416e-03,
    #    4.75250291e-03, 1.28075616e-04, 6.24942071e-04, 2.04546434e-05,
    #    9.22183305e-06, 2.36149361e-04, 1.13273939e-05, 8.74336221e-06,
    #    2.39752226e-04, 3.30081441e-04, 1.48277447e-05, 1.29773281e-04,
    #    2.87132439e-04, 1.61944521e-05, 3.96214298e-06, 1.25996697e-04,
    #    4.87309962e-06, 1.31103060e-04, 1.19993795e-04, 2.07571994e-04,
    #    1.07477687e-05, 1.22880354e-05, 2.13009427e-05, 1.89073677e-05,
    #    1.25476626e-03, 1.83125851e-05, 6.64609359e-04, 7.61778963e-04,
    #    2.66737074e-05])
    # # np.array([2.88915208e-05, 9.80923420e-06, 2.26321683e-04, 6.06076516e-05,
    # #    9.35931824e-05, 1.80927841e-04, 3.00750583e-05, 1.22766422e-05,
    # #    1.40985717e-05, 1.20683072e-05, 1.93956245e-05, 2.01386524e-05,
    # #    4.50497544e-06, 5.66952742e-06, 1.18097192e-05, 1.51887053e-05,
    # #    2.83487155e-05, 2.50380541e-05, 2.53424192e-05, 3.62951817e-05])
    # plt.plot(euvacAll, euvacIrr, 'ro', label='EUVAC')
    # plt.plot(euvacAll, ridleyIrr, 'co', label='RIDLEY')
    # plt.title('2016-1-19')
    # plt.legend(loc='best')
    # plt.show()

    return new_irr, ave_wav

#------------------------------------------------------------------------------
# main code:
#------------------------------------------------------------------------------
if __name__ == '__main__':
    args = parse_args_fism()

    if (args.fismfile == 'none'):

        if (args.start == '0'):
            print('Need to specify -start time! Use -h for help!')
            exit()

        start = convert_ymdhm_to_dt(args.start)

        if (args.flare):
            print('*** Downloading Flare data - Can only get 24 hours of data! ***')
            end = start + dt.timedelta(days = 1)
        else:
            if (args.end == '0'):
                print('Need to specify -end time! Use -h for help!')
                exit()
            end = convert_ymdhm_to_dt(args.end)

        fism_file = download_fism2(start, end, args.flare)

    else:
        fism_file = args.fismfile

    fism = read_fism_csv_file(fism_file)
    wavelengths = read_euv_csv_file(args.euvfile)

    nWaves = len(wavelengths['short'])

    filetime = fism["time"][0].strftime('fism%Y%m')
    filestart = filetime+'_nWaves_%03d' % nWaves

    fig = plt.figure(figsize = (10,10))
    ax = fig.add_subplot()

    if (args.gitm):
        fileout = filestart + '_gitm.dat'
    else:
        fileout = filestart + '.dat'

    fp = open(fileout, 'wb')

    if (args.gitm):
        fp.write('#START\n'.encode())

    else:

        shortline = ' 0000 00 00 00 00 00 '
        for short in wavelengths['short']:
            shortline = shortline + "%8.1f" % short
        shortline = shortline + '\n'
        fp.write(shortline.encode())

        longline = ' 0000 00 00 00 00 00 '
        for long in wavelengths['long']:
            longline = longline + "%8.1f" % long
        longline = longline + '\n'
        fp.write(longline.encode())

    for iTime, time in enumerate(fism['time']):
        new_irr, ave_wav = rebin_fism(fism['wave'], fism['irr'][iTime], wavelengths)

        if (args.gitm):
            ave_wav = np.flip(ave_wav)
            new_irr = np.flip(new_irr)
        ax.scatter(ave_wav, new_irr)

        sTime = time.strftime(' %Y %m %d %H %M %S')
        sData = ' '
        for irr in new_irr:
            sData = sData + "%15.8e" % irr
        line = sTime + sData + '\n'
        fp.write(line.encode())

    fp.close()

    ax.set_xlabel('Wavelength (A)')
    ax.set_ylabel(fism['vars'][2])

    plotfile = filestart + '.png'
    print('writing : ',plotfile)
    fig.savefig(plotfile)
    plt.close()
