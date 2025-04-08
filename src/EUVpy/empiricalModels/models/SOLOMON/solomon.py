# Code for computing solar irradiance according to Solomon and Qian 2005.
# Reference: Solomon, S. C. and Qian, L. (2005) Solar extreme-ultraviolet irradiance for general circulation models,
# Journal of Geophysical Research: Space Physics, 110, A10, 10.1029/2005JA011160

#-----------------------------------------------------------------------------------------------------------------------
# Top-level imports:
import numpy as np
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Local imports:
from EUVpy.tools.spectralAnalysis import spectralIrradiance
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Global Variable:
# solomonTable = np.array([
#         [1, 0.5, 4, 5.010e1, 0, 2.948e2, 5.010e1, 6.240e-1, 3.188e4, 7.847e5],
#         [2, 4, 8, 1.0e4, 0, 7.6e3, 1.0e4, 3.710e-1, 3.643e4, 8.968e5],
#         [3, 8, 18, 2.0e-6, 0, 4.6e5, 2.0e6, 2.0e-1, 5.485e6, 1.046e8],
#         [4, 18, 32, 7.600e6, 7.470e5, 9.220e5, 2.850e7, 6.247e-2, 6.317e6, 9.234e7],
#         [5, 32, 70, 1.659e8, 6.623e7, 4.293e6, 5.326e8, 1.343e-2, 3.710e8, 1.475e9],
#         [6, 70, 155, 4.012e8, 1.662e8, 5.678e6, 1.270e9, 9.182e-3, 1.023e9, 3.752e9],
#         [7, 155, 224, 2.078e9, 1.510e8, 6.273e7, 5.612e9, 1.433e-2, 2.953e9, 1.144e10],
#         [8, 224, 290, 1.724e9, 3.310e8, 9.834e7, 4.342e9, 2.575e-2, 4.927e9, 1.436e10],
#         [9, 290, 320, 6.793e9, 2.220e9, 4.286e7, 8.380e9, 7.059e-3, 6.942e9, 1.234e10],
#         [10, 320, 540, 2.750e9, 5.469e8, 1.080e8, 2.861e9, 1.458e-2, 6.486e9, 1.591e10],
#         [11, 540, 650, 5.035e9, 2.969e9, 1.590e7, 4.830e9, 5.857e-3, 3.499e9, 6.213e9],
#         [12, 650, 798, 1.562e9, 6.938e8, 8.208e6, 1.459e9, 5.719e-3, 1.869e9, 2.631e9],
#         [13, 650, 798, 1.264e9, 6.690e8, 5.445e5, 1.142e9, 3.680e-3, 1.136e9, 1.540e+09],
#         [14, 798, 913, 3.011e9, 3.011e9, 0, 2.364e9, 5.310e-3, 3.494e9, 5.868e9],
#         [15, 798, 913, 4.661e9, 4.213e9, 0, 3.655e9, 5.261e-3, 5.138e9, 8.562e9],
#         [16, 798, 913, 1.020e9, 1.020e9, 0, 8.448e8, 5.437e-3, 1.306e9, 2.157e9],
#         [17, 913, 975, 5.441e8, 4.187e8, 0, 3.818e8, 4.915e-3, 8.343e8, 1.373e9],
#         [18, 913, 975, 1.483e9, 1.307e9, 0, 1.028e9, 4.955e-3, 1.866e9, 2.862e9],
#         [19, 913, 975, 8.642e8, 8.440e8, 0, 7.156e8, 4.422e-3, 6.840e8, 1.111e9],
#         [20, 975, 987, 6.056e9, 3.671e9, 0, 4.482e9, 3.950e-3, 4.139e9, 6.801e9],
#         [21, 987, 1027, 5.569e9, 4.984e9, 0, 4.419e9, 5.021e-3, 6.274e9, 1.019e10],
#         [22, 1027, 1050, 6.309e9, 5.796e9, 0, 4.235e9, 4.825e-3, 4.389e9, 7.153e9]
#     ])

solomonTable = np.array([
    [1, 0.05, 0.4, 5.010e+01, 0.000e+00, 2.948e+02, 5.010e+01, 6.240e-01],
    [2, 0.4, 0.8, 1.000e+04, 0.000e+00, 7.600e+03, 1.000e+04, 3.710e-01],
    [3, 0.8, 1.8, 2.000e+06, 0.000e+00, 4.600e+05, 2.000e+06, 2.000e-01],
    [4, 1.8, 3.2, 7.600e+06, 7.470e+05, 9.220e+05, 2.850e+07, 6.247e-02],
    [5, 3.2, 7.0, 1.659e+08, 6.623e+07, 4.293e+06, 5.326e+08, 1.343e-02],
    [6, 7.0, 15.5, 4.012e+08, 1.662e+08, 5.678e+06, 1.270e+09, 9.182e-03],
    [7, 15.5, 22.4, 2.078e+09, 1.510e+08, 6.273e+07, 5.612e+09, 1.433e-02],
    [8, 22.4, 29.0, 1.724e+09, 3.310e+08, 9.834e+07, 4.342e+09, 2.575e-02],
    [9, 29.0, 32.0, 6.793e+09, 2.220e+09, 4.286e+07, 8.380e+09, 7.059e-03],
    [10, 32.0, 54.0, 2.750e+09, 5.469e+08, 1.080e+08, 2.861e+09, 1.458e-02],
    [11, 54.0, 65.0, 5.035e+09, 2.969e+09, 1.590e+07, 4.830e+09, 5.857e-03],
    [12, 65.0, 79.8, 1.562e+09, 6.938e+08, 8.208e+06, 1.459e+09, 5.719e-03],
    [13, 65.0, 79.8, 1.264e+09, 6.690e+08, 5.445e+05, 1.142e+09, 3.680e-03],
    [14, 79.8, 91.3, 3.011e+09, 3.011e+09, 0.000e+00, 2.364e+09, 5.310e-03],
    [15, 79.8, 91.3, 4.661e+09, 4.213e+09, 0.000e+00, 3.655e+09, 5.261e-03],
    [16, 79.8, 91.3, 1.020e+09, 1.020e+09, 0.000e+00, 8.448e+08, 5.437e-03],
    [17, 91.3, 97.5, 5.441e+08, 4.187e+08, 0.000e+00, 3.818e+08, 4.915e-03],
    [18, 91.3, 97.5, 1.483e+09, 1.307e+09, 0.000e+00, 1.028e+09, 4.955e-03],
    [19, 91.3, 97.5, 8.642e+08, 8.440e+08, 0.000e+00, 7.156e+08, 4.422e-03 ],
    [20, 97.5, 98.7, 6.056e+09, 3.671e+09, 0.000e+00, 4.482e+09, 3.950e-03],
    [21, 98.7, 102.7, 5.569e+09, 4.984e+09, 0.000e+00, 4.419e+09, 5.021e-03],
    [22, 102.7, 105.0, 6.309e+09, 5.796e+09, 0.000e+00, 4.235e+09, 4.825e-03]
    ])

solomonBands = {
    'short': solomonTable[:, 1]*10,
    'long': solomonTable[:, 2]*10
}
solomonBandWidths = np.array([3.5, 0.4, 1., 1.4, 3.8, 8.5, 6.9, 6.6, 3., 22., 11., 14.8, 14.8, 11.5, 11.5, 11.5, 6.2,
                              6.2, 6.2, 1.2,  4., 2.3 , 16.], dtype=np.float32)/10.
# The tables above are in units of Angstroms - source is Table A1 from Solomon and Qian 2005.
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
def solomon(F107, F107A, model='HFG'):
    """
    Compute the solar EUV irradiance in 23 standard bands.
    :param F107: ndarray
        Values of the F10.7 solar flux.
    :param F107A: ndarray
        Values of the 81-day averaged solar flux, centered on the present day.
    :param model: str
        Either 'HFG' or 'EUVAC'. Controls whether or not the empirical EUV data returned corresponds to the HFG model or
        the EUVAC model.
    :return solomonFlux: ndarray
        Values of the solar radiant flux in 23 distinct wavelength bands. Units of photon/m^2/s
    :return solomonIrr: ndarray
        Values of the solar EUV irradiance in 23 distinct wavelength bands. Units of W/m^2
    """
    # Instantiate the output data:
    if type(F107) == list or type(F107) != np.ndarray:
        F107 = np.array([F107])
        F107A = np.array([F107A])
        solomonFlux = np.zeros((1, solomonTable.shape[0]))
        solomonIrr = np.zeros((1, solomonTable.shape[0]))
    else:
        solomonFlux = np.zeros((len(F107), solomonTable.shape[0]))
        solomonIrr = np.zeros((len(F107), solomonTable.shape[0]))

    # Calculate the model coefficients (for every F107, F107A pair):
    r1 = 0.0138*(F107 - 71.5) + 0.005*(F107 - F107A + 3.9)
    r2 = 0.5943*(F107 - 71.5) + 0.381*(F107 - F107A + 3.9)

    # Compute P:
    P = 0.5*(F107  + F107A)

    # Loop across every F107, F107A pair:
    # for i in range(len(F107)):
    # Loop across every bin and fill in the data:
    for j in range(solomonIrr.shape[1]):
        waves = solomonTable[j, 1]
        wavel = solomonTable[j, 2]
        # dwav = wavel - waves
        mid = 0.5*(waves + wavel)
        if model == 'EUVAC':
            if j <= 3:
                prod = np.abs(1. + solomonTable[j, 7]*(P - 80.))
            else:
                prod = (1. + solomonTable[j, 7]*(P - 80.))
            flux = (solomonTable[j, 6] * prod) * 1.0e4 # convert from ph cm^-2 s^-1 to photon m^-2 s^-1
        else:
            flux = (solomonTable[j, 3] + r1*solomonTable[j, 4] + r2*solomonTable[j, 5]) * 1.0e4 # Units of photon m^-2 s^-1
        irradiance = spectralIrradiance(flux, mid*10)
        solomonFlux[:, j] = np.squeeze(flux)
        solomonIrr[:, j] = np.squeeze(irradiance)

    # mins = [np.nanmin(element) for element in solomonIrr.T]
    # np.where(np.asarray(mins) < 0)
    # 0, 1, 2, 3, 7 (EUVAC)
    # 0, 1, 2, 3, 4, 7, 9 (HFG)

    return solomonFlux, solomonIrr
#-----------------------------------------------------------------------------------------------------------------------